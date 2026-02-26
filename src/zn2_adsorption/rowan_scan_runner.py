"""
Rowan Zn-surface distance scan runner.

Submits 6 parallel scans -- one for each carbon in the hexagonal ring nearest
to the Zn starting position. Each scan constrains perpendicular approach via
two 90-degree angle constraints (DOF analysis in _build_scan_workflow docstring).
Graphene is frozen; C-FG bond lengths are constrained so FGs stay attached but can flex/rotate.

The neural network potential (NNP) model is configurable via the ``method``
parameter (default: UMA Small OMol25).
"""

from __future__ import annotations

import json
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np

try:
    import rowan
    import stjames
    from stjames import Molecule, Method, Settings
    from stjames.constraint import Constraint, ConstraintType
    from stjames.opt_settings import OptimizationSettings
    from stjames.workflows.scan import ScanSettings, ScanWorkflow
    ROWAN_AVAILABLE = True
except ImportError:
    ROWAN_AVAILABLE = False
    Molecule = None  # type: ignore[misc, assignment]
    Method = None
    Settings = None
    Constraint = None
    ConstraintType = None
    OptimizationSettings = None
    ScanSettings = None
    ScanWorkflow = None

HARTREE_TO_EV = 27.211386245988


def _structure_to_xyz_string(structure) -> str:
    """Convert pymatgen Structure to XYZ string (no PBC; for molecular calc)."""
    symbols = [site.specie.symbol for site in structure]
    positions = structure.cart_coords
    comment = "energy: 0.0"
    lines = [str(len(symbols)), comment]
    for sym, pos in zip(symbols, positions):
        lines.append(f"  {sym:4s} {pos[0]:18.12f} {pos[1]:18.12f} {pos[2]:18.12f}")
    return "\n".join(lines) + "\n"


def _resolve_method(method):
    """Resolve a method string or Method enum to a stjames Method.

    Accepts a ``stjames.Method`` enum member directly, or a case-insensitive
    string matching a ``Method`` name (e.g. ``"uma_s_omol"``, ``"GFN2_XTB"``).
    """
    if isinstance(method, Method):
        return method
    name = str(method).upper().strip()
    try:
        return Method[name]
    except KeyError:
        valid = ", ".join(m.name.lower() for m in Method)
        raise ValueError(
            f"Unknown method '{method}'. Valid options: {valid}"
        )


def _build_scan_workflow(
    struct_with_zn,
    n_graphene: int,
    target_c_idx: int,
    distances: List[float],
    neb_calculator,
    charge: int,
    multiplicity: int,
    freeze_surface: bool = True,
    freeze_functional_groups: bool = False,
    fg_bond_constraints: Optional[List[Tuple[int, int, float]]] = None,
    method=None,
) -> Tuple["ScanWorkflow", "Molecule"]:
    """
    Build a ScanWorkflow for Zn approaching a specific graphene carbon.

    Perpendicular approach is enforced via two angle constraints
    (Ci-C_target-Zn = 90 deg for two non-collinear graphene neighbors Ci).
    The scan bond fixes |Zn-C_target| = d (sphere, 2 angular DOF). Each 90 deg
    angle removes 1 DOF by confining Zn to the plane perpendicular to the
    C_target->Ci vector. Two such planes intersect along the surface normal,
    leaving Zn at exactly one point above C_target for each d.

    Args:
        struct_with_zn: pymatgen Structure with Zn as last atom.
        n_graphene: Number of graphene carbon atoms.
        target_c_idx: 0-based index of the target graphene carbon.
        distances: Scan distances (Angstrom).
        neb_calculator: NebCalculator for neighbor lookups.
        charge: Total system charge.
        multiplicity: Spin multiplicity.
        freeze_surface: If True, freeze all graphene carbons. If False,
            surface is free to relax; only angle constraints remain.
        freeze_functional_groups: If True, freeze all FG atoms (no motion).
            If False (default), FGs can flex/rotate; C-FG bond lengths are
            constrained via fg_bond_constraints.
        fg_bond_constraints: List of (anchor_c_idx, fg_atom_idx, bond_length).
            Constrains C-FG attachment bonds so they cannot stretch; FGs can
            still rotate/flex. Required when FGs present and not frozen.
        method: stjames.Method enum or string name. Default: UMA_S_OMOL.

    Returns:
        (ScanWorkflow, stjames Molecule).
    """
    zn_idx = len(struct_with_zn) - 1

    xyz_str = _structure_to_xyz_string(struct_with_zn)
    mol = Molecule.from_xyz(xyz_str, charge=charge, multiplicity=multiplicity)

    scan_settings = ScanSettings(
        type="bond",
        atoms=[zn_idx + 1, target_c_idx + 1],
        start=float(distances[0]),
        stop=float(distances[-1]),
        num=len(distances),
    )

    constraints: List[Constraint] = []
    # Atoms to freeze: graphene (optional) + functional groups (optional)
    freeze_indices: List[int] = []
    if freeze_surface:
        freeze_indices.extend(range(n_graphene))
    if freeze_functional_groups:
        # FG atoms are between n_graphene and zn_idx (exclusive)
        freeze_indices.extend(range(n_graphene, zn_idx))
    if freeze_indices:
        constraints.append(
            Constraint(
                constraint_type=ConstraintType.FREEZE_ATOMS,
                atoms=[i + 1 for i in freeze_indices],
            )
        )

    # C-FG bond constraints: keep attachment fixed, allow FG to flex/rotate
    if not freeze_functional_groups and fg_bond_constraints:
        for anchor_c, fg_atom, bond_len in fg_bond_constraints:
            if 0 <= anchor_c < zn_idx and 0 <= fg_atom < zn_idx:
                constraints.append(
                    Constraint(
                        constraint_type=ConstraintType.BOND,
                        atoms=[anchor_c + 1, fg_atom + 1],
                        value=float(bond_len),
                    )
                )

    neighbors = neb_calculator.surface_builder._get_neighbor_indices(
        struct_with_zn, {target_c_idx}, n_graphene
    )
    neighbor_list = sorted(neighbors)
    if len(neighbor_list) < 2:
        raise ValueError(
            f"Graphene C at index {target_c_idx} has {len(neighbor_list)} neighbor(s); "
            f"need at least 2 for perpendicular approach constraints."
        )
    for nbr_idx in neighbor_list[:2]:
        constraints.append(
            Constraint(
                constraint_type=ConstraintType.ANGLE,
                atoms=[nbr_idx + 1, target_c_idx + 1, zn_idx + 1],
                value=90.0,
            )
        )

    resolved = _resolve_method(method) if method is not None else Method.UMA_S_OMOL

    opt_settings = OptimizationSettings(constraints=constraints)
    calc_settings = Settings(
        method=resolved,
        tasks=["optimize"],
        corrections=[],
        mode="auto",
        opt_settings=opt_settings,
    )

    workflow = ScanWorkflow(
        initial_molecule=mol,
        scan_settings=scan_settings,
        calc_settings=calc_settings,
        calc_engine=resolved.default_engine(),
        wavefront_propagation=True,
    )
    return workflow, mol


def _submit_and_wait(
    workflow: "ScanWorkflow",
    mol: "Molecule",
    name: str,
    distances: List[float],
    verbose: bool,
) -> Dict[str, Any]:
    """Submit a single scan workflow, wait for result, return parsed output."""
    result = rowan.submit_workflow(
        workflow_type="scan",
        workflow_data=workflow.model_dump(mode="json"),
        initial_molecule=mol.model_dump(mode="json"),
        name=name,
    )
    workflow_url = f"https://labs.rowansci.com/scan/{result.uuid}"
    if verbose:
        print(f"  Submitted {name}: {workflow_url}")

    result.wait_for_result()
    result.fetch_latest(in_place=True)

    data = getattr(result, "data", None) or {}
    energies_hartree: List[float] = []
    distances_out: List[float] = list(distances)

    # Rowan result format varies; inspect and parse defensively
    if verbose and data:
        print(f"  Result keys: {list(data.keys())}")

    if "scan_points" in data:
        pts = data["scan_points"]
        if isinstance(pts, list) and pts:
            if isinstance(pts[0], dict):
                energies_hartree = [float(p.get("energy", p.get("E", 0))) for p in pts]
                if "coordinate" in pts[0] or "distance" in pts[0]:
                    distances_out = [
                        float(p.get("coordinate", p.get("distance", d)))
                        for p, d in zip(pts, distances)
                    ]
            elif isinstance(pts[0], (int, float)):
                energies_hartree = [float(e) for e in pts]
    elif "energies" in data:
        raw = data["energies"]
        if isinstance(raw, list):
            energies_hartree = [float(e) for e in raw]
    elif "energy" in data:
        energies_hartree = [float(data["energy"])] if len(distances) == 1 else []

    # Try alternative result structures from Rowan
    if not energies_hartree:
        for key in ("results", "calculations", "points"):
            items = data.get(key, [])
            if isinstance(items, list) and items:
                if isinstance(items[0], dict) and "energy" in items[0]:
                    energies_hartree = [float(p["energy"]) for p in items]
                    break

    energies_ev = [e * HARTREE_TO_EV for e in energies_hartree] if energies_hartree else []

    return {
        "distances": distances_out,
        "energies": energies_ev,
        "energies_hartree": energies_hartree,
        "workflow_uuid": str(result.uuid),
        "workflow_url": workflow_url,
    }


def run_distance_scan_rowan(
    surface_structure,
    n_graphene: int,
    fg_bond_constraints: List[Tuple[int, int, float]],
    distances: List[float],
    neb_calculator,
    charge: int = 2,
    multiplicity: int = 1,
    output_dir: Optional[Union[str, Path]] = None,
    name: str = "Zn-surface scan",
    verbose: bool = True,
    target_c_idx: Optional[int] = None,
    freeze_surface: bool = True,
    freeze_functional_groups: bool = False,
    method=None,
) -> Dict[str, Any]:
    """
    Run Zn-surface distance scan via Rowan for a single target carbon.

    If target_c_idx is None, the nearest graphene carbon to Zn is used.

    Args:
        surface_structure: pymatgen Structure with FGs but no Zn.
        n_graphene: Number of graphene carbon atoms.
        fg_bond_constraints: Unused; constraints use carbon surface only.
        distances: List of Zn-surface distances (Angstrom) to scan.
        neb_calculator: NebCalculator instance.
        charge: Total charge (default 2 for Zn2+).
        multiplicity: Unused; computed from electron count.
        output_dir: Directory for JSON output (optional).
        name: Workflow name in Rowan.
        verbose: If True, print progress.
        target_c_idx: 0-based index of the target graphene carbon.
            If None, nearest graphene carbon to Zn is used.
        freeze_surface: If True, freeze graphene. If False, surface relaxes.
        freeze_functional_groups: If True, freeze FG atoms entirely. If False
            (default), C-FG bonds are constrained; FGs can flex/rotate.
        method: stjames.Method enum or string name. Default: UMA_S_OMOL.

    Returns:
        Dict with distances, energies, workflow_uuid, workflow_url.
    """
    if not ROWAN_AVAILABLE:
        raise ImportError(
            "Rowan dependencies required. Install with: pip install zn2-adsorption[rowan]"
        )

    api_key = getattr(rowan, "api_key", None) or __import__("os").environ.get("ROWAN_API_KEY")
    if not api_key:
        raise RuntimeError("Rowan API key required. Set ROWAN_API_KEY or rowan.api_key.")

    # If no target specified, find nearest graphene C to center
    if target_c_idx is None:
        tmp = neb_calculator.surface_builder.add_zn2_ion(
            surface_structure.copy(), distance=distances[0]
        )
        tmp = neb_calculator._ensure_atom_ordering(tmp)
        zn_pos = tmp.cart_coords[-1]
        target_c_idx, _ = neb_calculator.surface_builder.find_nearest_surface_atom(
            tmp, zn_pos, graphene_only=True, n_graphene=n_graphene
        )

    # Place Zn directly above target carbon so initial geometry matches constraints
    struct_with_zn = neb_calculator.surface_builder.add_zn2_ion(
        surface_structure.copy(), distance=distances[0], above_atom_idx=target_c_idx,
    )
    struct_with_zn = neb_calculator._ensure_atom_ordering(struct_with_zn)

    ase_atoms = neb_calculator._pmg_to_ase(struct_with_zn)
    multiplicity = neb_calculator._calculate_multiplicity(ase_atoms, charge)

    resolved = _resolve_method(method) if method is not None else Method.UMA_S_OMOL

    if verbose:
        print(f"Submitting scan to Rowan ({resolved.name}), target C{target_c_idx + 1}...")
        print(f"  Points: {len(distances)}, range: {distances[0]:.2f}–{distances[-1]:.2f} Å")

    workflow, mol = _build_scan_workflow(
        struct_with_zn, n_graphene, target_c_idx,
        distances, neb_calculator, charge, multiplicity,
        freeze_surface=freeze_surface,
        freeze_functional_groups=freeze_functional_groups,
        fg_bond_constraints=fg_bond_constraints,
        method=resolved,
    )

    output = _submit_and_wait(workflow, mol, name, distances, verbose)

    if output_dir:
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        summary = {
            "path_method": "scan",
            "calculator": "rowan",
            "target_c_idx": target_c_idx,
            "distances": output["distances"],
            "energies": output["energies"],
            "workflow_url": output["workflow_url"],
            "workflow_uuid": output["workflow_uuid"],
        }
        with open(output_path / "scan_result.json", "w") as f:
            json.dump(summary, f, indent=2)

    return output


def run_hexagon_scans_rowan(
    surface_structure,
    n_graphene: int,
    fg_bond_constraints: List[Tuple[int, int, float]],
    distances: List[float],
    neb_calculator,
    charge: int = 2,
    multiplicity: int = 1,
    output_dir: Optional[Union[str, Path]] = None,
    verbose: bool = True,
    freeze_surface: bool = True,
    freeze_functional_groups: bool = False,
    method=None,
) -> Dict[str, Any]:
    """
    Run a single Zn-surface distance scan targeting the central carbon site
    that is away from functional groups (not an FG anchor or neighbor).

    Zn is placed above the chosen carbon. The target is the graphene carbon
    closest to the surface center that is not an FG attachment site and not
    adjacent to one. Uses perpendicular approach constraints (2 angle constraints).

    Args:
        surface_structure: pymatgen Structure with FGs but no Zn.
        n_graphene: Number of graphene carbon atoms.
        fg_bond_constraints: Unused.
        distances: List of Zn-surface distances (Angstrom).
        neb_calculator: NebCalculator instance.
        charge: Total charge (default 2).
        multiplicity: Unused; computed from electron count.
        output_dir: Directory for output (optional).
        verbose: If True, print progress.
        freeze_surface: If True, freeze graphene. If False, surface relaxes.
        freeze_functional_groups: If True, freeze FG atoms entirely. If False
            (default), C-FG bonds are constrained; FGs can flex/rotate.
        method: stjames.Method enum or string name. Default: UMA_S_OMOL.

    Returns:
        Dict with:
            - sites: list of one per-site result dict (target_c_idx, distances,
              energies, workflow_uuid, workflow_url)
            - hexagon_indices: [target_c_idx] for API compatibility
    """
    if not ROWAN_AVAILABLE:
        raise ImportError(
            "Rowan dependencies required. Install with: pip install zn2-adsorption[rowan]"
        )

    api_key = getattr(rowan, "api_key", None) or __import__("os").environ.get("ROWAN_API_KEY")
    if not api_key:
        raise RuntimeError("Rowan API key required. Set ROWAN_API_KEY or rowan.api_key.")

    # Find central carbon away from functional groups
    tmp_struct = neb_calculator.surface_builder.add_zn2_ion(
        surface_structure.copy(), distance=distances[0]
    )
    tmp_struct = neb_calculator._ensure_atom_ordering(tmp_struct)

    ase_atoms = neb_calculator._pmg_to_ase(tmp_struct)
    multiplicity = neb_calculator._calculate_multiplicity(ase_atoms, charge)

    c_coords = np.array([
        s.coords for s in tmp_struct if s.specie.symbol == "C"
    ])
    center_xy = tuple(np.mean(c_coords[:, :2], axis=0))
    fg_anchors = [i for i, _, _ in fg_bond_constraints] if fg_bond_constraints else []

    target_c_idx = neb_calculator.surface_builder.find_central_carbon_away_from_fg(
        tmp_struct, center_xy, n_graphene, fg_anchor_indices=fg_anchors,
    )
    hex_indices = [target_c_idx]

    resolved = _resolve_method(method) if method is not None else Method.UMA_S_OMOL

    if verbose:
        surface_mode = "surface frozen" if freeze_surface else "surface relaxed"
        print(f"Target carbon: C{target_c_idx + 1} (central, away from FGs)")
        print(f"Submitting scan to Rowan ({resolved.name}, {surface_mode})...")
        print(f"  Points: {len(distances)}, range: {distances[0]:.2f}–{distances[-1]:.2f} Å")

    workflows: List[Tuple[int, "ScanWorkflow", "Molecule"]] = []
    for c_idx in hex_indices:
        site_struct = neb_calculator.surface_builder.add_zn2_ion(
            surface_structure.copy(), distance=distances[0], above_atom_idx=c_idx,
        )
        site_struct = neb_calculator._ensure_atom_ordering(site_struct)
        wf, mol = _build_scan_workflow(
            site_struct, n_graphene, c_idx,
            distances, neb_calculator, charge, multiplicity,
            freeze_surface=freeze_surface,
            freeze_functional_groups=freeze_functional_groups,
            fg_bond_constraints=fg_bond_constraints,
            method=resolved,
        )
        workflows.append((c_idx, wf, mol))

    def _run_one(c_idx: int, wf: "ScanWorkflow", mol: "Molecule") -> Dict[str, Any]:
        site_name = f"Zn→C{c_idx + 1} scan"
        result = _submit_and_wait(wf, mol, site_name, distances, verbose)
        result["target_c_idx"] = c_idx
        return result

    site_results: List[Dict[str, Any]] = []
    with ThreadPoolExecutor(max_workers=6) as pool:
        futures = {
            pool.submit(_run_one, c_idx, wf, mol): c_idx
            for c_idx, wf, mol in workflows
        }
        for future in as_completed(futures):
            c_idx = futures[future]
            try:
                res = future.result()
                site_results.append(res)
                if verbose:
                    n_e = len(res.get("energies", []))
                    print(f"  C{c_idx + 1}: {n_e} energies, {res['workflow_url']}")
            except Exception as exc:
                if verbose:
                    print(f"  C{c_idx + 1}: FAILED - {exc}")
                site_results.append({
                    "target_c_idx": c_idx,
                    "distances": list(distances),
                    "energies": [],
                    "energies_hartree": [],
                    "workflow_uuid": None,
                    "workflow_url": None,
                    "error": str(exc),
                })

    site_results.sort(key=lambda r: r["target_c_idx"])

    output: Dict[str, Any] = {
        "hexagon_indices": hex_indices,
        "sites": site_results,
    }

    if output_dir:
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        with open(output_path / "hexagon_scan_results.json", "w") as f:
            json.dump(output, f, indent=2, default=str)

    return output
