"""
ASE Calculator wrapping Open Babel UFF for NEB and other workflows.

Provides get_potential_energy() and get_forces() for ASE optimizers.
Supports freeze constraints and FG-surface bond distance constraints.
"""

from __future__ import annotations

from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np

try:
    from ase import Atoms
    from ase.calculators.calculator import Calculator
    from ase.io import read, write
    ASE_AVAILABLE = True
except ImportError:
    ASE_AVAILABLE = False
    Atoms = None  # type: ignore[misc, assignment]
    Calculator = None  # type: ignore[misc, assignment]

try:
    from openbabel import openbabel

    OPENBABEL_AVAILABLE = True
    # Suppress "Failed to kekulize aromatic bonds" and similar warnings for
    # graphene XYZ files (no bond info); calculations still run correctly.
    openbabel.obErrorLog.SetOutputLevel(0)
except ImportError:
    OPENBABEL_AVAILABLE = False
    openbabel = None  # type: ignore[misc, assignment]

# Open Babel UFF uses kcal/mol and Angstrom. ASE uses eV and Angstrom.
KCAL_MOL_TO_EV = 0.043364  # 1 kcal/mol ≈ 0.043364 eV
# Forces: kcal/mol/Angstrom -> eV/Angstrom (same factor)


class UFFCalculator(Calculator):
    """
    ASE Calculator wrapping Open Babel UFF force field.

    Supports freeze constraints (fixed atoms) and distance constraints
    (e.g. FG-surface bonds). Charge/multiplicity are not used by UFF
    but can be set for compatibility.
    """

    implemented_properties = ["energy", "forces"]

    def __init__(
        self,
        freeze_indices: Optional[List[int]] = None,
        bond_constraint: Optional[Tuple[int, int, float]] = None,
        bond_constraints: Optional[List[Tuple[int, int, float]]] = None,
        max_force: Optional[float] = None,
        **kwargs,
    ):
        """
        Args:
            freeze_indices: 0-based atom indices to freeze.
            bond_constraint: Optional single (i, j, distance_angstrom).
            bond_constraints: Optional list of (i, j, d) for multiple bonds.
            max_force: Optional cap on force magnitude per atom (eV/Å). When set,
                forces exceeding this are scaled down to prevent NEB blow-up from
                bad intermediate geometries. Default 100.0 for NEB robustness.
            **kwargs: Passed to Calculator base class.
        """
        super().__init__(**kwargs)
        self.freeze_indices = freeze_indices or []
        self.bond_constraint = bond_constraint
        self.bond_constraints = bond_constraints or []
        if bond_constraint is not None:
            self.bond_constraints = list(self.bond_constraints) + [bond_constraint]
        self.max_force = max_force

    def calculate(
        self,
        atoms: Optional[Atoms] = None,
        properties: List[str] = None,
        system_changes: List[str] = None,
    ):
        if not ASE_AVAILABLE:
            raise ImportError("ASE is required. Install with: pip install ase")
        if not OPENBABEL_AVAILABLE:
            raise ImportError(
                "Open Babel is required for UFF. Install with: pip install openbabel-wheel"
            )

        if atoms is None:
            atoms = self.atoms
        if atoms is None:
            raise ValueError("No atoms supplied")

        if properties is None:
            properties = ["energy", "forces"]

        import tempfile
        with tempfile.NamedTemporaryFile(
            mode="w",
            suffix=".xyz",
            delete=False,
        ) as f:
            tmp_xyz = f.name
        try:
            write(tmp_xyz, atoms, format="xyz")
            ob_mol = openbabel.OBMol()
            conv = openbabel.OBConversion()
            conv.SetInFormat("xyz")
            if not conv.ReadFile(ob_mol, tmp_xyz):
                raise RuntimeError("Open Babel failed to read XYZ")

            constraints = openbabel.OBFFConstraints()
            for idx in self.freeze_indices:
                if 0 <= idx < ob_mol.NumAtoms():
                    constraints.AddAtomConstraint(idx + 1)
            for i, j, dist in self.bond_constraints:
                if 0 <= i < ob_mol.NumAtoms() and 0 <= j < ob_mol.NumAtoms():
                    constraints.AddDistanceConstraint(i + 1, j + 1, dist)

            ff = openbabel.OBForceField.FindForceField("UFF")
            if not ff:
                raise RuntimeError("Open Babel UFF force field not found")
            if not ff.Setup(ob_mol, constraints):
                raise RuntimeError("Open Babel UFF setup failed")
            ff.SetConstraints(constraints)

            energy_kcal = ff.Energy()

            n_atoms = ob_mol.NumAtoms()
            forces = np.zeros((n_atoms, 3))
            for i in range(1, n_atoms + 1):
                atom = ob_mol.GetAtom(i)
                g = ff.GetGradient(atom)
                if g:
                    # Force = -gradient; convert kcal/(mol*A) -> eV/A
                    forces[i - 1, 0] = -g.GetX() * KCAL_MOL_TO_EV
                    forces[i - 1, 1] = -g.GetY() * KCAL_MOL_TO_EV
                    forces[i - 1, 2] = -g.GetZ() * KCAL_MOL_TO_EV

            energy_eV = energy_kcal * KCAL_MOL_TO_EV
            # Clamp forces to prevent NEB blow-up when intermediate geometries
            # produce huge UFF repulsions (atoms too close).
            if self.max_force is not None and self.max_force > 0:
                for i in range(n_atoms):
                    f_norm = np.linalg.norm(forces[i])
                    if f_norm > self.max_force:
                        forces[i] = forces[i] * (self.max_force / f_norm)
            # Cap energy for robustness (avoids inf propagation)
            ENERGY_CAP_EV = 1e5
            if abs(energy_eV) > ENERGY_CAP_EV:
                energy_eV = np.sign(energy_eV) * ENERGY_CAP_EV
            self.results = {
                "energy": energy_eV,
                "forces": forces,
            }
        finally:
            Path(tmp_xyz).unlink(missing_ok=True)
