"""Tests for the distance scan runner."""

import pytest

try:
    from openbabel import openbabel
    OPENBABEL_AVAILABLE = True
except ImportError:
    OPENBABEL_AVAILABLE = False

from zn2_adsorption.neb_calculator import NebCalculator
from zn2_adsorption.scan_runner import run_distance_scan_uff

pytest.importorskip("ase")


@pytest.mark.skipif(not OPENBABEL_AVAILABLE, reason="Open Babel not installed")
def test_run_distance_scan_uff_basic():
    """Test distance scan produces distances and energies."""
    calc = NebCalculator()
    base = calc.prepare_scan_base(
        num_carboxyl=0,
        num_hydroxyl=0,
        surface_size=(2, 2),
    )
    distances = [5.0, 4.0, 3.0]

    result = run_distance_scan_uff(
        surface_structure=base["surface_structure"],
        n_graphene=base["n_graphene"],
        fg_bond_constraints=base["fg_bond_constraints"],
        distances=distances,
        surface_builder=calc.surface_builder,
        neb_calculator=calc,
        output_dir=None,
        verbose=False,
        save_geometries=False,
    )

    assert result["distances"] == distances
    assert len(result["energies"]) == 3
    assert all(e == e for e in result["energies"])
    assert len(result["geometries"]) == 3


@pytest.mark.skipif(not OPENBABEL_AVAILABLE, reason="Open Babel not installed")
def test_run_distance_scan_uff_output_dir(tmp_path):
    """Test scan writes XYZ files when output_dir is set."""
    calc = NebCalculator()
    base = calc.prepare_scan_base(
        num_carboxyl=0,
        num_hydroxyl=0,
        surface_size=(2, 2),
    )
    distances = [5.0, 4.0]

    run_distance_scan_uff(
        surface_structure=base["surface_structure"],
        n_graphene=base["n_graphene"],
        fg_bond_constraints=base["fg_bond_constraints"],
        distances=distances,
        surface_builder=calc.surface_builder,
        neb_calculator=calc,
        output_dir=tmp_path,
        verbose=False,
        save_geometries=True,
    )

    assert (tmp_path / "scan_000.xyz").exists()
    assert (tmp_path / "scan_001.xyz").exists()
