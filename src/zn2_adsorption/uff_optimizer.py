"""
UFF (Universal Force Field) geometry optimizer for endpoint structures.

Uses Open Babel's UFF implementation for fast, constrained geometry optimization.
Carbon surface and Zn ion are frozen; only functional groups relax.
Runs in seconds instead of hours compared to DFT (ORCA).
"""

from __future__ import annotations

from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np

try:
    from ase import Atoms
    from ase.io import read, write
    ASE_AVAILABLE = True
except ImportError:
    ASE_AVAILABLE = False
    Atoms = None  # type: ignore[misc, assignment]

try:
    from openbabel import openbabel

    OPENBABEL_AVAILABLE = True
    # Suppress "Failed to kekulize aromatic bonds" and similar warnings for
    # graphene XYZ files (no bond info); calculations still run correctly.
    openbabel.obErrorLog.SetOutputLevel(0)
except ImportError:
    OPENBABEL_AVAILABLE = False
    openbabel = None  # type: ignore[misc, assignment]


def optimize_with_uff(
    atoms: Atoms,
    freeze_indices: Optional[List[int]] = None,
    bond_constraint: Optional[Tuple[int, int, float]] = None,
    bond_constraints: Optional[List[Tuple[int, int, float]]] = None,
    steps: int = 500,
) -> bool:
    """
    Optimize geometry using UFF with optional constraints.

    Args:
        atoms: ASE Atoms object. Positions are updated in-place on success.
        freeze_indices: 0-based atom indices to freeze (e.g. graphene carbons).
        bond_constraint: Optional single (atom1_idx, atom2_idx, distance_angstrom).
        bond_constraints: Optional list of (i, j, d) for multiple bonds (e.g. FG-surface).
        steps: Maximum conjugate gradient steps (default 500).

    Returns:
        True if optimization succeeded, False otherwise.
    """
    if not ASE_AVAILABLE:
        raise ImportError("ASE is required. Install with: pip install ase")
    if not OPENBABEL_AVAILABLE:
        raise ImportError(
            "Open Babel is required for UFF. Install with: pip install openbabel-wheel"
        )

    import tempfile
    with tempfile.NamedTemporaryFile(
        mode="w",
        suffix=".xyz",
        delete=False,
    ) as f:
        tmp_xyz = f.name
    try:
        write(tmp_xyz, atoms, format="xyz")
        obMol = openbabel.OBMol()
        conv = openbabel.OBConversion()
        conv.SetInFormat("xyz")
        if not conv.ReadFile(obMol, tmp_xyz):
            return False

        constraints = openbabel.OBFFConstraints()
        if freeze_indices:
            for idx in freeze_indices:
                if 0 <= idx < obMol.NumAtoms():
                    constraints.AddAtomConstraint(idx + 1)
        if bond_constraint is not None:
            i, j, dist = bond_constraint
            if 0 <= i < obMol.NumAtoms() and 0 <= j < obMol.NumAtoms():
                constraints.AddDistanceConstraint(i + 1, j + 1, dist)
        if bond_constraints:
            for i, j, dist in bond_constraints:
                if 0 <= i < obMol.NumAtoms() and 0 <= j < obMol.NumAtoms():
                    constraints.AddDistanceConstraint(i + 1, j + 1, dist)

        forcefield = openbabel.OBForceField.FindForceField("UFF")
        if not forcefield:
            return False
        if not forcefield.Setup(obMol, constraints):
            return False
        forcefield.SetConstraints(constraints)

        forcefield.ConjugateGradients(steps)
        forcefield.GetCoordinates(obMol)

        conv.SetOutFormat("xyz")
        conv.WriteFile(obMol, tmp_xyz)
        opt_atoms = read(tmp_xyz, format="xyz")
        atoms.positions[:] = opt_atoms.get_positions()
        return True
    finally:
        Path(tmp_xyz).unlink(missing_ok=True)


def get_endpoint_constraints(
    atoms: Atoms,
    n_graphene: int,
) -> Tuple[List[int], None]:
    """
    Compute freeze indices for endpoint optimization (surface + Zn frozen).

    Only functional groups are allowed to relax; carbon surface and Zn ion
    stay fixed to avoid changing Zn height or surface geometry.

    Args:
        atoms: ASE Atoms (Zn must be last).
        n_graphene: Number of graphene carbon atoms (indices 0..n_graphene-1).

    Returns:
        (freeze_indices, None) where freeze_indices includes graphene + Zn.
    """
    zn_idx = len(atoms) - 1
    freeze_indices = list(range(n_graphene)) + [zn_idx]
    return freeze_indices, None
