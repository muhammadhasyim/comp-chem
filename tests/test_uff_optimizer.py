"""Tests for UFF optimizer module."""

import pytest

try:
    from ase import Atoms
    from zn2_adsorption.uff_optimizer import optimize_with_uff, get_endpoint_constraints
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
def test_optimize_with_uff_simple():
    """Test UFF optimization on simple system."""
    atoms = Atoms(
        symbols=["H", "H"],
        positions=[[0.0, 0.0, 0.0], [0.0, 0.0, 1.0]],
    )
    result = optimize_with_uff(atoms, freeze_indices=None, bond_constraint=None)
    assert result is True
    assert atoms.get_positions().shape == (2, 3)


@pytest.mark.skipif(not OPENBABEL_AVAILABLE, reason="Open Babel not installed")
def test_optimize_with_uff_constrained():
    """Test UFF optimization with atom freeze and distance constraint."""
    atoms = Atoms(
        symbols=["C", "C", "Zn"],
        positions=[[0.0, 0.0, 0.0], [1.4, 0.0, 0.0], [0.7, 0.0, 3.0]],
    )
    result = optimize_with_uff(
        atoms,
        freeze_indices=[0, 1],
        bond_constraint=(2, 0, 3.0),
    )
    assert result is True


def test_get_endpoint_constraints():
    """Test constraint extraction: surface + Zn frozen, no bond constraint."""
    atoms = Atoms(
        symbols=["C", "C", "O", "Zn"],
        positions=[[0, 0, 0], [1.4, 0, 0], [0.7, 1.2, 0], [0.7, 0.5, 3.0]],
    )
    freeze, bc = get_endpoint_constraints(atoms, n_graphene=2)
    assert freeze == [0, 1, 3]
    assert bc is None
