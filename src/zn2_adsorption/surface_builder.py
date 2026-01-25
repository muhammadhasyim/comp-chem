"""
Build functionalized carbon surfaces for adsorption studies.

This module creates carbon surface models with functional groups suitable
for Zn2+ adsorption energy calculations, using pymatgen for periodic
boundary conditions and adapting GOPY's functionalization logic.
"""

from typing import List, Tuple, Optional, Dict, Union
import numpy as np
import random

try:
    from pymatgen.core import Structure, Lattice, Molecule, Element
    from pymatgen.core.surface import SlabGenerator
    PYMATGEN_AVAILABLE = True
except ImportError:
    PYMATGEN_AVAILABLE = False

class FunctionalizedGrapheneBuilder:
    """
    Build functionalized graphene surfaces with periodic boundary conditions.
    """
    
    def __init__(
        self,
        miller_indices: Tuple[int, int, int] = (0, 0, 1),
        min_slab_size: float = 1.0, # Graphene is 1 layer
        min_vacuum_size: float = 15.0
    ):
        if not PYMATGEN_AVAILABLE:
            raise ImportError("pymatgen is required. Install with: pip install pymatgen")
        self.miller_indices = miller_indices
        self.min_slab_size = min_slab_size
        self.min_vacuum_size = min_vacuum_size
        
        # GOPY typical bond lengths
        self.bond_lengths = {
            'C-C': 1.42,
            'C-OH': 1.49,
            'O-H': 0.97,
            'C-COOH': 1.52,
            'C=O': 1.21,
            'C-O_acid': 1.32,
        }

    def build_pristine_surface(
        self,
        supercell_size: Tuple[int, int] = (4, 4)
    ) -> Structure:
        """
        Build pristine graphene slab with periodic boundary conditions.
        """
        # Graphene unit cell (hexagonal)
        a = 2.46
        lattice = Lattice.hexagonal(a, 10.0) # c=10.0 is dummy
        structure = Structure(lattice, ["C", "C"], [[1/3, 2/3, 0], [2/3, 1/3, 0]])
        
        # Create slab
        slab_gen = SlabGenerator(
            structure,
            self.miller_indices,
            self.min_slab_size,
            self.min_vacuum_size,
            center_slab=True
        )
        slab = slab_gen.get_slab()
        
        # Make supercell
        slab.make_supercell([supercell_size[0], supercell_size[1], 1])
        
        return slab

    def add_functional_groups(
        self,
        structure: Structure,
        num_carboxyl: int,
        num_hydroxyl: int,
        num_epoxy: int = 0
    ) -> Structure:
        """
        Add functional groups using adapted GOPY algorithms.
        """
        new_structure = structure.copy()
        
        # 1. Add hydroxyl groups (-OH)
        for _ in range(num_hydroxyl):
            new_structure = self._add_single_hydroxyl(new_structure)
            
        # 2. Add carboxyl groups (-COOH)
        for _ in range(num_carboxyl):
            new_structure = self._add_single_carboxyl(new_structure)
            
        return new_structure

    def _add_single_hydroxyl(self, structure: Structure) -> Structure:
        # Find a carbon atom that isn't already functionalized
        # (Simplified: just pick a random C)
        c_indices = [i for i, s in enumerate(structure.species) if s.symbol == "C"]
        if not c_indices:
            return structure
            
        target_idx = random.choice(c_indices)
        target_pos = structure.cart_coords[target_idx]
        
        # Direction: up (z-axis)
        normal = np.array([0, 0, 1])
        
        # O atom
        o_pos = target_pos + normal * self.bond_lengths['C-OH']
        # H atom (at an angle, simplified to up for now)
        h_pos = o_pos + np.array([0.5, 0, 0.8]) # Dummy offset
        h_pos = o_pos + (h_pos - o_pos) / np.linalg.norm(h_pos - o_pos) * self.bond_lengths['O-H']
        
        structure.append("O", o_pos, coords_are_cartesian=True)
        structure.append("H", h_pos, coords_are_cartesian=True)
        
        return structure

    def _add_single_carboxyl(self, structure: Structure) -> Structure:
        c_indices = [i for i, s in enumerate(structure.species) if s.symbol == "C"]
        if not c_indices:
            return structure
            
        target_idx = random.choice(c_indices)
        target_pos = structure.cart_coords[target_idx]
        
        normal = np.array([0, 0, 1])
        
        # C of COOH
        c_acid_pos = target_pos + normal * self.bond_lengths['C-COOH']
        # O (C=O)
        o1_pos = c_acid_pos + np.array([1.0, 0, 0.5])
        o1_pos = c_acid_pos + (o1_pos - c_acid_pos) / np.linalg.norm(o1_pos - c_acid_pos) * self.bond_lengths['C=O']
        # O (C-OH)
        o2_pos = c_acid_pos + np.array([-1.0, 0, 0.5])
        o2_pos = c_acid_pos + (o2_pos - c_acid_pos) / np.linalg.norm(o2_pos - c_acid_pos) * self.bond_lengths['C-O_acid']
        # H
        h_pos = o2_pos + np.array([0, 0, 0.97])
        
        structure.append("C", c_acid_pos, coords_are_cartesian=True)
        structure.append("O", o1_pos, coords_are_cartesian=True)
        structure.append("O", o2_pos, coords_are_cartesian=True)
        structure.append("H", h_pos, coords_are_cartesian=True)
        
        return structure

    def add_zn2_ion(
        self,
        structure: Structure,
        distance: float,
        position: Optional[Tuple[float, float]] = None
    ) -> Structure:
        """
        Add Zn2+ ion at specified distance from surface.
        """
        # Find center of mass or average position of surface
        c_coords = np.array([s.coords for s in structure if s.specie.symbol == "C"])
        avg_pos = np.mean(c_coords, axis=0)
        
        if position is not None:
            # Use provided fractional (x, y)
            frac_pos = [position[0], position[1], 0.5] # Start at middle z
            cart_pos = structure.lattice.get_cartesian_coords(frac_pos)
            zn_pos = np.array([cart_pos[0], cart_pos[1], avg_pos[2] + distance])
        else:
            # Default to center
            zn_pos = np.array([avg_pos[0], avg_pos[1], avg_pos[2] + distance])
            
        structure.append("Zn", zn_pos, coords_are_cartesian=True)
        return structure

    def find_nearest_surface_atom(
        self,
        structure: Structure,
        zn2_position: np.ndarray
    ) -> Tuple[int, float]:
        """
        Find nearest carbon/oxygen atom to Zn2+ position.
        """
        min_dist = float('inf')
        nearest_idx = -1
        
        for i, site in enumerate(structure):
            if site.specie.symbol in ["C", "O"]:
                dist = structure.lattice.get_distance_and_image(
                    structure.lattice.get_fractional_coords(zn2_position),
                    site.frac_coords
                )[0]
                if dist < min_dist:
                    min_dist = dist
                    nearest_idx = i
                    
        return nearest_idx, min_dist
