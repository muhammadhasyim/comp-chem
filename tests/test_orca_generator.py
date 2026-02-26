import pytest
from ase import Atoms
from zn2_adsorption.orca_generator import OrcaInputGenerator

def test_orca_generator_basic():
    atoms = Atoms(symbols=['H', 'H'], positions=[[0, 0, 0], [0, 0, 0.74]])
    gen = OrcaInputGenerator(method="B3LYP", basis_set="def2-SVP")
    input_content = gen.generate_from_ase(atoms, calc_type="opt")
    assert "!B3LYP def2-SVP" in input_content
    assert "OPT" in input_content
    assert "* xyz 0 1" in input_content
    assert "H " in input_content
    assert "0.0000000000" in input_content
    assert "0.7400000000" in input_content

def test_orca_generator_constraints():
    atoms = Atoms(symbols=['Zn', 'C'], positions=[[0, 0, 3.0], [0, 0, 0]])
    gen = OrcaInputGenerator()
    input_content = gen.generate_from_ase(
        atoms, 
        calc_type="opt",
        atom1_idx=0,
        atom2_idx=1,
        distance=3.0
    )
    assert "%geom" in input_content
    assert "Constraints" in input_content
    assert "{ B 0 1 3.000000 C }" in input_content

def test_orca_generator_freeze_atoms():
    """Test ORCA input with freeze_atoms constraint."""
    atoms = Atoms(symbols=['C', 'C', 'H', 'H'], positions=[[0, 0, 0], [1.4, 0, 0], [0, 0, 1], [1.4, 0, 1]])
    gen = OrcaInputGenerator()
    input_content = gen.generate_from_ase(
        atoms,
        calc_type="opt",
        freeze_atoms=[0, 1]
    )
    assert "%geom" in input_content
    assert "Constraints" in input_content
    assert "{ C 0 C }" in input_content
    assert "{ C 1 C }" in input_content


def test_orca_generator_combined_constraints():
    """Test ORCA input with both freeze_atoms and distance constraint."""
    atoms = Atoms(symbols=['C', 'O', 'Zn'], positions=[[0, 0, 0], [0, 0, 1.2], [0, 0, 3.0]])
    gen = OrcaInputGenerator()
    input_content = gen.generate_from_ase(
        atoms,
        calc_type="opt",
        atom1_idx=2,
        atom2_idx=1,
        distance=3.0,
        freeze_atoms=[0]
    )
    assert "%geom" in input_content
    assert "Constraints" in input_content
    assert "{ C 0 C }" in input_content
    assert "{ B 2 1 3.000000 C }" in input_content


def test_orca_generator_solvation():
    atoms = Atoms(symbols=['H', 'H'], positions=[[0, 0, 0], [0, 0, 0.74]])
    gen = OrcaInputGenerator(solvent="water", solvation_model="CPCM")
    input_content = gen.generate_from_ase(atoms, calc_type="sp")
    assert "CPCM(water)" in input_content
    assert "%cpcm" in input_content
    assert "solvent water" in input_content
