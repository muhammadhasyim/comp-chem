"""Tests for UFF ASE Calculator."""

import pytest

try:
    from ase import Atoms
    from zn2_adsorption.uff_calculator import UFFCalculator, KCAL_MOL_TO_EV
    ASE_AVAILABLE = True
except ImportError:
    ASE_AVAILABLE = False
    pytestmark = pytest.mark.skip("ASE not available")

try:
    from openbabel import openbabel
    OPENBABEL_AVAILABLE = True
except ImportError:
    OPENBABEL_AVAILABLE = False


@pytest.mark.skipif(not OPENBABEL_AVAILABLE, reason="Open Babel not installed")
def test_uff_calculator_energy_forces():
    """Test UFFCalculator returns energy and forces."""
    atoms = Atoms("H2", positions=[[0, 0, 0], [0, 0, 0.74]])
    atoms.calc = UFFCalculator()
    e = atoms.get_potential_energy()
    f = atoms.get_forces()
    assert e > 0
    assert f.shape == (2, 3)
    assert abs(f.sum()) < 1e-6  # Net force near zero


@pytest.mark.skipif(not OPENBABEL_AVAILABLE, reason="Open Babel not installed")
def test_uff_calculator_constrained():
    """Test UFFCalculator with freeze and bond constraints."""
    atoms = Atoms(
        symbols=["C", "C", "Zn"],
        positions=[[0, 0, 0], [1.4, 0, 0], [0.7, 0, 3.0]],
    )
    atoms.calc = UFFCalculator(
        freeze_indices=[0, 1],
        bond_constraints=[(2, 0, 3.0)],
    )
    e = atoms.get_potential_energy()
    f = atoms.get_forces()
    assert e is not None
    assert f.shape == (3, 3)
