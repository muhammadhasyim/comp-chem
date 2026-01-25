import pytest
import numpy as np
from pymatgen.core import Structure
from zn2_adsorption.surface_builder import FunctionalizedGrapheneBuilder

def test_builder_pristine():
    builder = FunctionalizedGrapheneBuilder()
    structure = builder.build_pristine_surface(supercell_size=(2, 2))
    assert isinstance(structure, Structure)
    assert len(structure) == 8 # 2 atoms per unit cell * 2 * 2 supercell
    assert structure.lattice is not None

def test_builder_functionalized():
    builder = FunctionalizedGrapheneBuilder()
    structure = builder.build_pristine_surface(supercell_size=(4, 4))
    functionalized = builder.add_functional_groups(
        structure,
        num_carboxyl=1,
        num_hydroxyl=1
    )
    assert len(functionalized) > len(structure)
    symbols = [s.symbol for s in functionalized.species]
    assert 'O' in symbols
    assert 'H' in symbols

def test_add_zn2():
    builder = FunctionalizedGrapheneBuilder()
    structure = builder.build_pristine_surface(supercell_size=(4, 4))
    zn_structure = builder.add_zn2_ion(structure, distance=3.0)
    assert 'Zn' in [s.symbol for s in zn_structure.species]
    
    # Check distance
    zn_idx = [i for i, s in enumerate(zn_structure.species) if s.symbol == 'Zn'][0]
    zn_pos = zn_structure.cart_coords[zn_idx]
    
    # Find surface z
    c_coords = np.array([s.coords for s in zn_structure if s.specie.symbol == "C"])
    surface_z = np.mean(c_coords[:, 2])
    
    assert abs(zn_pos[2] - surface_z - 3.0) < 0.1
