"""
Build functionalized carbon surfaces for adsorption studies.

This module creates carbon surface models with functional groups suitable
for Zn2+ adsorption energy calculations, using pymatgen for periodic
boundary conditions and adapting GOPY's functionalization logic.
"""

from __future__ import annotations

from typing import Dict, Iterable, List, Optional, Tuple, Union

# Type for FG-surface bond: (anchor_c_idx, fg_first_atom_idx, bond_length_angstrom)
FGBond = Tuple[int, int, float]
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
        
        # Bond lengths (Å): GOPY / literature values
        # Carboxylate (-COO⁻): C–C to surface, C–O ~1.26 (resonance)
        self.bond_lengths = {
            'C-C': 1.42,
            'C-OH': 1.49,
            'O-H': 0.97,
            'C-COOH': 1.52,  # alias; carboxylate uses C-COO
            'C-COO': 1.52,   # graphene C to carboxylate C
            'C-O_carboxylate': 1.26,  # C–O in -COO⁻
            'C=O': 1.21,
            'C-O_acid': 1.32,
        }
        # Cutoff for nearest-neighbor (avoid placing FGs on adjacent C)
        self._nn_cutoff = 1.5  # Å, slightly above C-C 1.42
        # Clustering around Zn: exclusion and annulus
        self._r_exclude = 1.2  # Å, exclude carbons directly below Zn
        self._r_cluster = 5.0  # Å, prefer FGs within this lateral radius

    def _get_carbon_indices_near_zn(
        self,
        structure: Structure,
        zn_center_xy: Tuple[float, float],
        n_graphene: int,
        forbidden: set[int],
        r_exclude: Optional[float] = None,
        r_cluster: Optional[float] = None,
    ) -> List[int]:
        """
        Return eligible carbon indices sorted by lateral distance from zn_center_xy.
        Excludes carbons with r < r_exclude (directly below Zn).
        Prefers carbons with r in [r_exclude, r_cluster]; falls back to farther sites.
        """
        r_ex = r_exclude if r_exclude is not None else self._r_exclude
        r_cl = r_cluster if r_cluster is not None else self._r_cluster
        cx, cy = zn_center_xy
        candidates: List[Tuple[float, int]] = []
        for i in range(n_graphene):
            if structure.species[i].symbol != "C" or i in forbidden:
                continue
            pos = structure.cart_coords[i]
            r = np.sqrt((pos[0] - cx) ** 2 + (pos[1] - cy) ** 2)
            if r < r_ex:
                continue
            candidates.append((r, i))
        # Sort by distance (closest first)
        candidates.sort(key=lambda x: x[0])
        # Prefer annulus [r_exclude, r_cluster]; if none, use all
        in_annulus = [i for r, i in candidates if r <= r_cl]
        eligible = in_annulus if in_annulus else [i for _, i in candidates]
        return eligible

    def _get_neighbor_indices(
        self, structure: Structure, site_indices: set, n_sites: int
    ) -> set:
        """Return indices of all sites that are nearest neighbors of given indices."""
        neighbors = set()
        for i in site_indices:
            for j in range(n_sites):
                if j == i:
                    continue
                if structure[i].specie.symbol != "C" or structure[j].specie.symbol != "C":
                    continue
                d = structure.get_distance(i, j)
                if d < self._nn_cutoff:
                    neighbors.add(j)
        return neighbors

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
        num_epoxy: int = 0,
        zn_center_xy: Optional[Tuple[float, float]] = None,
    ) -> Tuple[Structure, List[FGBond]]:
        """
        Add functional groups using adapted GOPY algorithms.
        Hydroxyl and carboxylate (-COO⁻) occupy distinct carbon sites (no overlap).
        FGs are placed to avoid nearest-neighbor carbons for uniform distribution.

        If zn_center_xy is provided, FGs are clustered around that (x,y) position:
        carbons directly below (r < r_exclude) are excluded; prefer sites in
        annulus [r_exclude, r_cluster]. When None, use random placement.

        Returns:
            (structure, fg_bonds) where fg_bonds is a list of (anchor_c_idx,
            fg_first_atom_idx, bond_length_angstrom) for each FG-surface bond.
        """
        new_structure = structure.copy()
        n_graphene = len(new_structure)
        used_carbon_indices: set[int] = set()
        fg_bonds: List[FGBond] = []

        # 1. Add hydroxyl groups (-OH)
        for _ in range(num_hydroxyl):
            new_structure, bond = self._add_single_hydroxyl(
                new_structure, used_carbon_indices, n_graphene,
                zn_center_xy=zn_center_xy,
            )
            fg_bonds.append(bond)

        # 2. Add carboxylate groups (-COO⁻)
        for _ in range(num_carboxyl):
            new_structure, bond = self._add_single_carboxylate(
                new_structure, used_carbon_indices, n_graphene,
                zn_center_xy=zn_center_xy,
            )
            fg_bonds.append(bond)

        return new_structure, fg_bonds

    def _add_single_hydroxyl(
        self,
        structure: Structure,
        used_carbon_indices: set[int],
        n_graphene: int,
        zn_center_xy: Optional[Tuple[float, float]] = None,
    ) -> Tuple[Structure, FGBond]:
        forbidden = used_carbon_indices | self._get_neighbor_indices(
            structure, used_carbon_indices, n_graphene
        )
        if zn_center_xy is not None:
            c_indices = self._get_carbon_indices_near_zn(
                structure, zn_center_xy, n_graphene, forbidden
            )
        else:
            c_indices = [
                i
                for i in range(n_graphene)
                if structure.species[i].symbol == "C" and i not in forbidden
            ]
        if not c_indices:
            raise ValueError(
                "Not enough distinct carbon sites for hydroxyl groups"
            )

        target_idx = random.choice(c_indices) if zn_center_xy is None else c_indices[0]
        used_carbon_indices.add(target_idx)
        target_pos = structure.cart_coords[target_idx]
        
        # Direction: up (z-axis)
        normal = np.array([0, 0, 1])
        bond_len = self.bond_lengths['C-OH']
        
        # O atom (first atom of this FG)
        o_pos = target_pos + normal * bond_len
        # H atom (at an angle, simplified to up for now)
        h_pos = o_pos + np.array([0.5, 0, 0.8]) # Dummy offset
        h_pos = o_pos + (h_pos - o_pos) / np.linalg.norm(h_pos - o_pos) * self.bond_lengths['O-H']
        
        structure.append("O", o_pos, coords_are_cartesian=True)
        structure.append("H", h_pos, coords_are_cartesian=True)
        # O is at index len-2 after append; with one hydroxyl, O at n_graphene
        fg_first_idx = len(structure) - 2
        return structure, (target_idx, fg_first_idx, bond_len)

    def _add_single_carboxylate(
        self,
        structure: Structure,
        used_carbon_indices: set[int],
        n_graphene: int,
        zn_center_xy: Optional[Tuple[float, float]] = None,
    ) -> Tuple[Structure, FGBond]:
        """Add a carboxylate group (-COO⁻): C + 2 O atoms, no H."""
        forbidden = used_carbon_indices | self._get_neighbor_indices(
            structure, used_carbon_indices, n_graphene
        )
        if zn_center_xy is not None:
            c_indices = self._get_carbon_indices_near_zn(
                structure, zn_center_xy, n_graphene, forbidden
            )
        else:
            c_indices = [
                i
                for i in range(n_graphene)
                if structure.species[i].symbol == "C" and i not in forbidden
            ]
        if not c_indices:
            raise ValueError(
                "Not enough distinct carbon sites for carboxylate groups"
            )

        target_idx = random.choice(c_indices) if zn_center_xy is None else c_indices[0]
        used_carbon_indices.add(target_idx)
        target_pos = structure.cart_coords[target_idx]
        
        normal = np.array([0, 0, 1])
        bond_len = self.bond_lengths['C-COO']
        c_o_len = self.bond_lengths['C-O_carboxylate']
        
        # Carboxylate C (first atom of this FG), above graphene
        c_acid_pos = target_pos + normal * bond_len
        # Two O atoms: O-C-O angle 120°, symmetric, tilted up for Zn²⁺ access
        o1_dir = np.array([0.866, 0, 0.5]) * c_o_len   # 30° from horizontal
        o2_dir = np.array([-0.866, 0, 0.5]) * c_o_len  # 120° from o1
        o1_pos = c_acid_pos + o1_dir
        o2_pos = c_acid_pos + o2_dir
        
        structure.append("C", c_acid_pos, coords_are_cartesian=True)
        structure.append("O", o1_pos, coords_are_cartesian=True)
        structure.append("O", o2_pos, coords_are_cartesian=True)
        # C_acid is at index len-3 after append (3 atoms: C, O, O)
        fg_first_idx = len(structure) - 3
        return structure, (target_idx, fg_first_idx, bond_len)

    def add_zn2_ion(
        self,
        structure: Structure,
        distance: float,
        position: Optional[Tuple[float, float]] = None,
        above_atom_idx: Optional[int] = None,
    ) -> Structure:
        """
        Add Zn2+ ion at specified distance above the surface.

        Args:
            structure: pymatgen Structure.
            distance: Height above the surface plane (Angstrom).
            position: (frac_x, frac_y) to place Zn. Ignored if above_atom_idx set.
            above_atom_idx: If given, place Zn directly above this atom's
                (x, y) at z = atom_z + distance. Used for per-site scans.
        """
        c_coords = np.array([s.coords for s in structure if s.specie.symbol == "C"])
        avg_z = float(np.mean(c_coords[:, 2]))

        if above_atom_idx is not None:
            target_pos = structure.cart_coords[above_atom_idx]
            zn_pos = np.array([target_pos[0], target_pos[1], target_pos[2] + distance])
        elif position is not None:
            frac_pos = [position[0], position[1], 0.5]
            cart_pos = structure.lattice.get_cartesian_coords(frac_pos)
            zn_pos = np.array([cart_pos[0], cart_pos[1], avg_z + distance])
        else:
            avg_pos = np.mean(c_coords, axis=0)
            zn_pos = np.array([avg_pos[0], avg_pos[1], avg_z + distance])

        structure.append("Zn", zn_pos, coords_are_cartesian=True)
        return structure

    def find_hexagon_ring(
        self,
        structure: Structure,
        center_xy: Tuple[float, float],
        n_graphene: int,
    ) -> List[int]:
        """
        Return indices of the 6 graphene carbons forming the hexagonal ring
        closest to center_xy. These are the nearest-neighbor shell at ~1.42 A
        lateral distance from the ring center.

        Args:
            structure: pymatgen Structure (with or without Zn/FGs).
            center_xy: (x, y) cartesian coordinates of the ring center.
            n_graphene: Number of graphene carbon atoms (indices 0..n_graphene-1).

        Returns:
            List of 6 carbon indices sorted by lateral distance.
        """
        cx, cy = center_xy
        dists: List[Tuple[float, int]] = []
        for i in range(n_graphene):
            if structure.species[i].symbol != "C":
                continue
            pos = structure.cart_coords[i]
            r = np.sqrt((pos[0] - cx) ** 2 + (pos[1] - cy) ** 2)
            dists.append((r, i))
        dists.sort()
        if len(dists) < 6:
            raise ValueError(
                f"Only {len(dists)} graphene carbons found; need at least 6 for hexagon ring."
            )
        return [idx for _, idx in dists[:6]]

    def find_central_carbon_away_from_fg(
        self,
        structure: Structure,
        center_xy: Tuple[float, float],
        n_graphene: int,
        fg_anchor_indices: Optional[Iterable[int]] = None,
    ) -> int:
        """
        Return the graphene carbon closest to center_xy that is not an FG
        attachment site and not adjacent to one. Targets the "most middle" site
        away from functional groups.

        Args:
            structure: pymatgen Structure (with or without Zn/FGs).
            center_xy: (x, y) cartesian coordinates of the surface center.
            n_graphene: Number of graphene carbon atoms.
            fg_anchor_indices: Indices of carbons with FGs attached. These and
                their nearest neighbors are excluded.

        Returns:
            Index of the preferred central carbon.
        """
        forbidden: set[int] = set()
        if fg_anchor_indices:
            anchors = set(fg_anchor_indices) & set(range(n_graphene))
            forbidden = anchors | self._get_neighbor_indices(
                structure, anchors, n_graphene
            )
        cx, cy = center_xy
        candidates: List[Tuple[float, int]] = []
        for i in range(n_graphene):
            if structure.species[i].symbol != "C" or i in forbidden:
                continue
            pos = structure.cart_coords[i]
            r = np.sqrt((pos[0] - cx) ** 2 + (pos[1] - cy) ** 2)
            candidates.append((r, i))
        if not candidates:
            raise ValueError(
                "No central carbon site away from functional groups. "
                "All carbons are FG anchors or neighbors. Try a larger surface."
            )
        candidates.sort()
        return candidates[0][1]

    def find_nearest_surface_atom(
        self,
        structure: Structure,
        zn2_position: np.ndarray,
        graphene_only: bool = False,
        n_graphene: Optional[int] = None,
    ) -> Tuple[int, float]:
        """
        Find nearest carbon/oxygen atom to Zn2+ position.

        Args:
            structure: pymatgen Structure.
            zn2_position: Cartesian coordinates of Zn.
            graphene_only: If True, consider only graphene carbons (indices < n_graphene).
            n_graphene: Required when graphene_only=True; number of graphene atoms.

        Returns:
            (nearest_index, distance).
        """
        min_dist = float('inf')
        nearest_idx = -1
        symbols = ["C", "O"] if not graphene_only else ["C"]

        for i, site in enumerate(structure):
            if site.specie.symbol not in symbols:
                continue
            if graphene_only and n_graphene is not None and i >= n_graphene:
                continue
            dist = structure.lattice.get_distance_and_image(
                structure.lattice.get_fractional_coords(zn2_position),
                site.frac_coords
            )[0]
            if dist < min_dist:
                min_dist = dist
                nearest_idx = i

        return nearest_idx, min_dist
