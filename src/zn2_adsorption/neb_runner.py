"""
NEB (Nudged Elastic Band) runner using UFF via ASE.

Runs NEB optimization in-process with Open Babel UFF as the calculator.
Supports constrained optimization (freeze surface+Zn, FG-surface bonds).
"""

from __future__ import annotations

from pathlib import Path
from typing import List, Optional, Tuple, Dict, Any

import numpy as np

try:
    from ase import Atoms
    from ase.constraints import FixAtoms
    from ase.io import read, write
    from ase.mep.neb import NEB
    from ase.optimize import FIRE
    ASE_AVAILABLE = True
except ImportError:
    ASE_AVAILABLE = False
    Atoms = None  # type: ignore[misc, assignment]
    FixAtoms = None  # type: ignore[misc, assignment]
    NEB = None  # type: ignore[misc, assignment]
    FIRE = None  # type: ignore[misc, assignment]

from zn2_adsorption.uff_calculator import UFFCalculator


# Maximum energy (eV) before NEB is considered diverged; stop early to avoid runaway
DIVERGENCE_ENERGY_THRESHOLD_EV = 1e4


def run_neb_uff(
    initial: Atoms,
    final: Atoms,
    n_images: int = 7,
    freeze_indices: Optional[List[int]] = None,
    fg_bond_constraints: Optional[List[Tuple[int, int, float]]] = None,
    interpolation: str = "idpp",
    fmax: float = 0.05,
    max_steps: int = 200,
    output_dir: Optional[Path] = None,
    verbose: bool = True,
    max_force: float = 100.0,
    maxstep: float = 0.05,
) -> Dict[str, Any]:
    """
    Run NEB optimization using UFF.

    Args:
        initial: Initial structure (reactant).
        final: Final structure (product).
        n_images: Total number of images including endpoints (minimum 3).
        freeze_indices: 0-based atom indices to freeze (e.g. surface + Zn).
        fg_bond_constraints: List of (i, j, d) for FG-surface bond constraints.
        interpolation: 'linear' or 'idpp' for initial path guess.
        fmax: Force convergence criterion (eV/Å).
        max_steps: Maximum optimization steps.
        output_dir: Directory to write optimized images (neb_000.xyz, ...).
        verbose: If True, print optimization progress (step, energy, fmax).
        max_force: Cap on force magnitude (eV/Å) to prevent blow-up from bad geometries.
        maxstep: Max atom displacement per step (Å); smaller = more stable.

    Returns:
        Dict with:
            - images: List of optimized ASE Atoms
            - energies: List of energies (eV) per image
            - ts_index: Index of highest-energy image (transition state estimate)
            - converged: Whether optimization converged
    """
    if not ASE_AVAILABLE:
        raise ImportError("ASE is required. Install with: pip install ase")

    freeze_indices = freeze_indices or []
    fg_bond_constraints = fg_bond_constraints or []

    if len(initial) != len(final):
        raise ValueError(
            f"Initial and final must have same atom count: {len(initial)} vs {len(final)}"
        )
    if initial.get_chemical_symbols() != final.get_chemical_symbols():
        raise ValueError("Initial and final must have identical atom ordering")

    n_images = max(3, n_images)
    n_interior = n_images - 2  # Exclude endpoints

    # Build image list: [initial, interior..., final]
    images = [initial.copy()]
    for _ in range(n_interior):
        img = initial.copy()
        img.set_positions(initial.get_positions())
        images.append(img)
    images.append(final.copy())

    # Interpolate interior images
    neb = NEB(images)
    neb.interpolate(method=interpolation)

    # Apply constraints and calculator to each image (fresh calc per image)
    for img in images:
        if freeze_indices:
            img.set_constraint(FixAtoms(indices=freeze_indices))
        img.calc = UFFCalculator(
            freeze_indices=freeze_indices,
            bond_constraints=fg_bond_constraints,
            max_force=max_force,
        )

    # Optimize with step-by-step divergence check
    import sys
    logfile = sys.stdout if verbose else None
    optimizer = FIRE(
        neb,
        trajectory=None,
        logfile=logfile,
        maxstep=maxstep,
        downhill_check=True,
    )
    converged = False
    for is_converged in optimizer.irun(fmax=fmax, steps=max_steps):
        energies = [img.get_potential_energy() for img in images]
        max_energy = max(energies)
        if max_energy > DIVERGENCE_ENERGY_THRESHOLD_EV:
            if verbose:
                print(
                    f"  NEB diverged (max E = {max_energy:.1e} eV > {DIVERGENCE_ENERGY_THRESHOLD_EV:.0e}); stopping early."
                )
            break
        converged = is_converged
        if is_converged:
            break

    # Extract results
    energies = [img.get_potential_energy() for img in images]
    ts_index = int(np.argmax(energies))

    if output_dir:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        for i, img in enumerate(images):
            write(output_dir / f"neb_{i:03d}.xyz", img, format="xyz")

    return {
        "images": images,
        "energies": energies,
        "ts_index": ts_index,
        "converged": converged,
    }
