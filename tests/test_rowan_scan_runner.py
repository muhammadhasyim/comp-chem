"""Tests for the Rowan distance scan runner."""

import pytest
from unittest.mock import MagicMock, patch

from zn2_adsorption.neb_calculator import NebCalculator
from zn2_adsorption.rowan_scan_runner import (
    ROWAN_AVAILABLE,
    _structure_to_xyz_string,
    run_distance_scan_rowan,
)

pytest.importorskip("pymatgen")


def test_structure_to_xyz_string():
    """Test pymatgen structure to XYZ string conversion."""
    calc = NebCalculator()
    base = calc.prepare_scan_base(
        num_carboxyl=0,
        num_hydroxyl=0,
        surface_size=(2, 2),
    )
    struct = calc.surface_builder.add_zn2_ion(
        base["surface_structure"].copy(),
        distance=5.0,
    )
    xyz = _structure_to_xyz_string(struct)
    lines = xyz.strip().split("\n")
    assert lines[0].strip() == str(len(struct))
    assert "energy:" in lines[1]
    # Line 2+: Symbol x y z
    for i, line in enumerate(lines[2:], start=0):
        parts = line.split()
        assert len(parts) == 4
        assert parts[0] in ("C", "H", "Zn")


@pytest.mark.skipif(not ROWAN_AVAILABLE, reason="Rowan dependencies not installed")
def test_run_distance_scan_rowan_workflow_structure():
    """Test that run_distance_scan_rowan constructs correct workflow with constraints."""
    calc = NebCalculator()
    base = calc.prepare_scan_base(
        num_carboxyl=0,
        num_hydroxyl=0,
        surface_size=(2, 2),
    )
    distances = [5.0, 4.0]

    mock_workflow = MagicMock()
    mock_workflow.uuid = "test-uuid-123"
    mock_workflow.data = {
        "scan_points": [
            {"energy": -100.0, "coordinate": 5.0},
            {"energy": -101.0, "coordinate": 4.0},
        ],
    }

    with patch.dict("os.environ", {"ROWAN_API_KEY": "test-key"}):
        with patch(
            "zn2_adsorption.rowan_scan_runner.rowan.submit_workflow",
            return_value=mock_workflow,
        ) as mock_submit:
            with patch.object(mock_workflow, "wait_for_result"):
                with patch.object(mock_workflow, "fetch_latest"):
                    result = run_distance_scan_rowan(
                        surface_structure=base["surface_structure"],
                        n_graphene=base["n_graphene"],
                        fg_bond_constraints=base["fg_bond_constraints"],
                        distances=distances,
                        neb_calculator=calc,
                        charge=2,
                        multiplicity=1,
                        verbose=False,
                    )

    # Verify submit_workflow was called with scan workflow
    call_args = mock_submit.call_args
    assert call_args.kwargs.get("workflow_type") == "scan"
    workflow_data = call_args.kwargs.get("workflow_data", {})
    assert "scan_settings" in workflow_data
    assert "calc_settings" in workflow_data
    # scan_settings is a list (stjames serialization)
    scan_cfg = workflow_data["scan_settings"][0] if isinstance(workflow_data["scan_settings"], list) else workflow_data["scan_settings"]
    assert scan_cfg["type"] == "bond"
    assert len(scan_cfg["atoms"]) == 2
    opt_settings = workflow_data["calc_settings"].get("opt_settings", {})
    constraints = opt_settings.get("constraints", [])
    # At least freeze_atoms for surface
    freeze = [c for c in constraints if c.get("constraint_type") == "freeze_atoms"]
    assert len(freeze) >= 1
    assert len(freeze[0].get("atoms", [])) == base["n_graphene"]

    # Verify result shape
    assert result["distances"] == [5.0, 4.0]
    assert len(result["energies"]) == 2
    assert result["workflow_uuid"] == "test-uuid-123"
    assert "labs.rowansci.com" in result["workflow_url"]


@pytest.mark.skipif(not ROWAN_AVAILABLE, reason="Rowan dependencies not installed")
def test_run_distance_scan_rowan_freezes_functional_groups():
    """Test that freeze_functional_groups=True includes FG atoms in freeze constraint."""
    from zn2_adsorption.rowan_scan_runner import _build_scan_workflow

    import random
    random.seed(42)
    calc = NebCalculator()
    base = calc.prepare_scan_base(
        num_carboxyl=1,
        num_hydroxyl=1,
        surface_size=(3, 3),
    )
    struct = calc.surface_builder.add_zn2_ion(
        base["surface_structure"].copy(),
        distance=5.0,
    )
    struct = calc._ensure_atom_ordering(struct)
    n_graphene = base["n_graphene"]
    zn_idx = len(struct) - 1
    n_fg = zn_idx - n_graphene  # 1 carboxylate (3 atoms) + 1 hydroxyl (2 atoms) = 5

    # Charge = 2 - n_carboxyl (carboxylate is -1 each)
    charge = 2 - 1  # 1 carboxylate
    wf_frozen, _ = _build_scan_workflow(
        struct, n_graphene, 0, [5.0], calc, charge, 1,
        freeze_surface=True,
        freeze_functional_groups=True,
        fg_bond_constraints=base["fg_bond_constraints"],
    )
    wf_bonded, _ = _build_scan_workflow(
        struct, n_graphene, 0, [5.0], calc, charge, 1,
        freeze_surface=True,
        freeze_functional_groups=False,
        fg_bond_constraints=base["fg_bond_constraints"],
    )

    # Serialize to get constraint atoms (stjames uses 1-based)
    data_frozen = wf_frozen.model_dump(mode="json")
    data_bonded = wf_bonded.model_dump(mode="json")
    constraints_frozen = data_frozen.get("calc_settings", {}).get("opt_settings", {}).get("constraints", [])
    constraints_bonded = data_bonded.get("calc_settings", {}).get("opt_settings", {}).get("constraints", [])

    freeze_frozen = next((c for c in constraints_frozen if c.get("constraint_type") == "freeze_atoms"), None)
    freeze_bonded = next((c for c in constraints_bonded if c.get("constraint_type") == "freeze_atoms"), None)
    bond_constraints = [c for c in constraints_bonded if c.get("constraint_type") == "bond"]
    assert freeze_frozen is not None
    assert freeze_bonded is not None

    # With freeze_functional_groups=True: freeze surface + FGs = n_graphene + n_fg
    assert len(freeze_frozen.get("atoms", [])) == n_graphene + n_fg
    # With freeze_functional_groups=False: freeze surface only; C-FG bonds constrained
    assert len(freeze_bonded.get("atoms", [])) == n_graphene
    assert len(bond_constraints) == 2  # 1 carboxylate + 1 hydroxyl C-FG bonds


def test_run_distance_scan_rowan_requires_api_key():
    """Test that run_distance_scan_rowan raises without API key."""
    if not ROWAN_AVAILABLE:
        pytest.skip("Rowan dependencies not installed")

    import zn2_adsorption.rowan_scan_runner as mod

    calc = NebCalculator()
    base = calc.prepare_scan_base(
        num_carboxyl=0,
        num_hydroxyl=0,
        surface_size=(2, 2),
    )

    with patch.dict("os.environ", {"ROWAN_API_KEY": ""}):
        with patch.object(mod.rowan, "api_key", None, create=True):
            with pytest.raises(RuntimeError, match="API key"):
                run_distance_scan_rowan(
                    surface_structure=base["surface_structure"],
                    n_graphene=base["n_graphene"],
                    fg_bond_constraints=base["fg_bond_constraints"],
                    distances=[5.0],
                    neb_calculator=calc,
                    verbose=False,
                )


def test_run_distance_scan_rowan_import_error_when_unavailable():
    """Test helpful error when Rowan deps not installed."""
    if ROWAN_AVAILABLE:
        pytest.skip("Rowan is installed, cannot test ImportError path")

    calc = NebCalculator()
    base = calc.prepare_scan_base(
        num_carboxyl=0,
        num_hydroxyl=0,
        surface_size=(2, 2),
    )

    with pytest.raises(ImportError, match="rowan"):
        run_distance_scan_rowan(
            surface_structure=base["surface_structure"],
            n_graphene=base["n_graphene"],
            fg_bond_constraints=base["fg_bond_constraints"],
            distances=[5.0],
            neb_calculator=calc,
            verbose=False,
        )
