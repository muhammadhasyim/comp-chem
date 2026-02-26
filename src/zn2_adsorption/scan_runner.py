"""
Zn-surface distance scan runner using UFF.

Performs constrained geometry optimization at each distance. Surface and Zn
are frozen. C-FG bond lengths are constrained so FGs stay attached but can
flex/rotate. Use --freeze-functional-groups to fully freeze FGs (no motion).
Output is E(d) potential curve.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import numpy as np

try:
    from ase import Atoms
    from ase.io import write
    ASE_AVAILABLE = True
except ImportError:
    ASE_AVAILABLE = False
    Atoms = None  # type: ignore[misc, assignment]

from zn2_adsorption.uff_optimizer import optimize_with_uff
from zn2_adsorption.uff_calculator import UFFCalculator


def run_distance_scan_uff(
    surface_structure,
    n_graphene: int,
    fg_bond_constraints: List[Tuple[int, int, float]],
    distances: List[float],
    surface_builder,
    neb_calculator,
    output_dir: Optional[Union[str, Path]] = None,
    verbose: bool = True,
    save_geometries: bool = True,
    use_continuation: bool = True,
    freeze_functional_groups: bool = False,
) -> Dict:
    """
    Run Zn-surface distance scan with UFF optimization at each distance.

    At each distance d: place Zn at d, freeze surface and Zn (and optionally
    FGs), record energy.

    Args:
        surface_structure: pymatgen Structure with FGs but no Zn.
        n_graphene: Number of graphene carbon atoms.
        fg_bond_constraints: List of (i, j, d) for FG-surface bond constraints.
        distances: List of Zn-surface distances (Angstrom) to scan.
        surface_builder: FunctionalizedGrapheneBuilder instance for add_zn2_ion.
        neb_calculator: NebCalculator instance for _pmg_to_ase, _ensure_atom_ordering.
        output_dir: Directory to write scan_*.xyz files.
        verbose: If True, print progress per distance.
        save_geometries: If True and output_dir set, save optimized XYZ per point.
        use_continuation: If True, use optimized geometry from previous distance as
            starting guess for next (only Zn position updated). Faster and more stable.
        freeze_functional_groups: If True, freeze FG atoms entirely. If False
            (default), C-FG bonds are constrained; FGs can flex/rotate.

    Returns:
        Dict with distances, energies, and optionally geometries.
    """
    if not ASE_AVAILABLE:
        raise ImportError("ASE is required. Install with: pip install ase")

    output_dir = Path(output_dir) if output_dir else None
    if output_dir and save_geometries:
        output_dir.mkdir(parents=True, exist_ok=True)

    energies: List[float] = []
    geometries: List[Atoms] = []

    for i, d in enumerate(distances):
        if i == 0 or not use_continuation or not geometries:
            # Fresh structure: surface + FGs + Zn at distance d
            struct_with_zn = surface_builder.add_zn2_ion(
                surface_structure.copy(),
                distance=d,
            )
            struct_with_zn = neb_calculator._ensure_atom_ordering(struct_with_zn)
            atoms = neb_calculator._pmg_to_ase(struct_with_zn)
        else:
            # Continuation: use previous optimized geometry, move Zn closer
            atoms = geometries[-1].copy()
            avg_pos = np.mean(atoms.positions[:n_graphene], axis=0)
            atoms.positions[-1] = np.array([
                avg_pos[0],
                avg_pos[1],
                avg_pos[2] + d,
            ])

        zn_idx = len(atoms) - 1
        if freeze_functional_groups:
            freeze_indices = list(range(zn_idx + 1))  # surface + FGs + Zn
        else:
            freeze_indices = list(range(n_graphene)) + [zn_idx]  # surface + Zn

        if verbose:
            print(f"  Scan point {i + 1}/{len(distances)}: d = {d:.3f} Angstrom ...", end=" ", flush=True)

        success = optimize_with_uff(
            atoms,
            freeze_indices=freeze_indices,
            bond_constraints=fg_bond_constraints,
        )

        if success:
            # Get energy via UFFCalculator (optimizer doesn't attach one)
            calc = UFFCalculator(
                freeze_indices=freeze_indices,
                bond_constraints=fg_bond_constraints,
            )
            atoms.calc = calc
            energy = atoms.get_potential_energy()
            energies.append(energy)
            geometries.append(atoms.copy())
            if verbose:
                print(f"E = {energy:.4f} eV")
            if output_dir and save_geometries:
                write(output_dir / f"scan_{i:03d}.xyz", atoms, format="xyz")
        else:
            energies.append(np.nan)
            geometries.append(atoms.copy())
            if verbose:
                print("(optimization failed)")
            if output_dir and save_geometries:
                write(output_dir / f"scan_{i:03d}.xyz", atoms, format="xyz")

    return {
        "distances": list(distances),
        "energies": energies,
        "geometries": geometries,
    }
