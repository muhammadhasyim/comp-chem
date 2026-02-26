"""
Nudged Elastic Band (NEB) calculator for Zn2+ adsorption pathways.

This module handles the complete workflow for NEB calculations including:
- Building initial and final structures with Zn²⁺ at different distances
- Constrained pre-optimization of endpoints
- NEB-TS input generation
- Path analysis
"""

from __future__ import annotations

from typing import Optional, Dict, Tuple, List, Union
from pathlib import Path
import numpy as np

try:
    from ase import Atoms
    from ase.io import read, write
    ASE_AVAILABLE = True
except ImportError:
    ASE_AVAILABLE = False
    Atoms = None  # type: ignore[misc, assignment]

try:
    from pymatgen.core import Structure
    PYMATGEN_AVAILABLE = True
except ImportError:
    PYMATGEN_AVAILABLE = False
    Structure = None  # type: ignore[misc, assignment]

from zn2_adsorption.orca_generator import OrcaInputGenerator
from zn2_adsorption.surface_builder import FunctionalizedGrapheneBuilder


class NebCalculator:
    """
    Orchestrate the Zn2+ NEB calculation workflow.
    
    Creates initial and final structures with Zn²⁺ at different distances,
    pre-optimizes both endpoints with constraints, and generates NEB-TS input.
    """
    
    def __init__(
        self,
        orca_generator: Optional[OrcaInputGenerator] = None,
        surface_builder: Optional[FunctionalizedGrapheneBuilder] = None,
        method: str = "B3LYP",
        basis_set: str = "def2-TZVP",
        solvent: str = "water",
        solvation_model: str = "CPCM"
    ):
        self.orca_gen = orca_generator or OrcaInputGenerator(
            method=method,
            basis_set=basis_set,
            solvent=solvent,
            solvation_model=solvation_model
        )
        self.surface_builder = surface_builder or FunctionalizedGrapheneBuilder()

    def prepare_neb_calculation(
        self,
        num_carboxyl: int = 0,
        num_hydroxyl: int = 0,
        surface_size: Tuple[int, int] = (4, 4),
        start_distance: float = 5.0,
        end_distance: float = 2.5,
        neb_images: Optional[int] = None,
        output_dir: Optional[Union[str, Path]] = None,
        constrain_endpoints: bool = True
    ) -> Dict[str, str]:
        """
        Prepare NEB calculation workflow.

        Args:
            num_carboxyl: Number of carboxyl groups
            num_hydroxyl: Number of hydroxyl groups
            surface_size: Surface supercell size (nx, ny)
            start_distance: Initial Zn²⁺ distance from surface (Å)
            end_distance: Final Zn²⁺ distance from surface (Å)
            neb_images: Number of NEB images (None = auto-calculate)
            output_dir: Directory for output files
            constrain_endpoints: If True, freeze carbon surface and fix Zn-surface
                distance during endpoint optimization (avoids steric hindrance
                by relaxing only functional groups). Default True.

        Returns:
            Dictionary with paths to generated input files
        """
        if output_dir:
            output_dir = Path(output_dir)
            output_dir.mkdir(parents=True, exist_ok=True)
        
        # 1. Build functionalized surface (common to both endpoints)
        pristine = self.surface_builder.build_pristine_surface(supercell_size=surface_size)
        c_coords = np.array([s.coords for s in pristine if s.specie.symbol == "C"])
        zn_center = (float(np.mean(c_coords[:, 0])), float(np.mean(c_coords[:, 1])))
        surface_structure, _ = self.surface_builder.add_functional_groups(
            pristine,
            num_carboxyl=num_carboxyl,
            num_hydroxyl=num_hydroxyl,
            zn_center_xy=zn_center,
        )
        
        # 2. Create initial structure (Zn²⁺ at start_distance)
        initial_structure = self.create_initial_structure(
            surface_structure.copy(),
            start_distance
        )
        
        # 3. Create final structure (Zn²⁺ at end_distance)
        final_structure = self.create_final_structure(
            surface_structure.copy(),
            end_distance
        )
        
        # 4. Ensure atom ordering is consistent (Zn²⁺ must be last)
        initial_structure = self._ensure_atom_ordering(initial_structure)
        final_structure = self._ensure_atom_ordering(final_structure)
        
        # 5. Convert to ASE
        initial_ase = self._pmg_to_ase(initial_structure)
        final_ase = self._pmg_to_ase(final_structure)
        
        # 6. Set up constraints for endpoint optimization
        # Carbon surface = indices 0..n_graphene-1; functional groups + Zn = flexible
        n_graphene = len(pristine)
        freeze_atoms = list(range(n_graphene)) if constrain_endpoints else None
        
        def _endpoint_constraints(atoms_ase, structure_pmg):
            zn_idx = len(atoms_ase) - 1
            zn_pos = structure_pmg.cart_coords[zn_idx]
            nearest_idx, actual_dist = self.surface_builder.find_nearest_surface_atom(
                structure_pmg, zn_pos
            )
            return zn_idx, nearest_idx, actual_dist
        
        initial_zn, initial_nearest, initial_dist = _endpoint_constraints(
            initial_ase, initial_structure
        )
        final_zn, final_nearest, final_dist = _endpoint_constraints(
            final_ase, final_structure
        )
        
        # 7. Generate optimization inputs for endpoints
        # Charge: Zn²⁺ (+2) + n carboxylate (-1 each) = 2 - num_carboxyl
        charge = 2 - num_carboxyl
        self.orca_gen.charge = charge
        self.orca_gen.multiplicity = self._calculate_multiplicity(initial_ase, charge=charge)
        
        # Initial endpoint: freeze surface + fix Zn-surface distance; relax functional groups
        initial_input = self.orca_gen.generate_from_ase(
            initial_ase,
            calc_type="opt",
            output_file=output_dir / "initial.inp" if output_dir else None,
            title="Initial endpoint (Zn²⁺ at start distance)",
            atom1_idx=initial_zn if constrain_endpoints else None,
            atom2_idx=initial_nearest if constrain_endpoints else None,
            distance=initial_dist if constrain_endpoints else None,
            freeze_atoms=freeze_atoms
        )
        
        # Final endpoint: same constraints
        final_input = self.orca_gen.generate_from_ase(
            final_ase,
            calc_type="opt",
            output_file=output_dir / "final.inp" if output_dir else None,
            title="Final endpoint (Zn²⁺ at end distance)",
            atom1_idx=final_zn if constrain_endpoints else None,
            atom2_idx=final_nearest if constrain_endpoints else None,
            distance=final_dist if constrain_endpoints else None,
            freeze_atoms=freeze_atoms
        )
        
        # 8. Export structures as XYZ files (unoptimized; will be overwritten after opt)
        if output_dir:
            self._export_xyz(initial_ase, output_dir / "initial.xyz")
            self._export_xyz(final_ase, output_dir / "final.xyz")
        
        # 9. Generate NEB input file (initially with unoptimized structures)
        # This will be regenerated after optimizations complete
        neb_input = self.orca_gen.generate_neb_input(
            initial_atoms=initial_ase,
            final_atoms=final_ase,
            output_file=output_dir / "neb.inp" if output_dir else None,
            neb_images=neb_images
        )
        
        return {
            "initial_input": initial_input,
            "final_input": final_input,
            "neb_input": neb_input,
            "initial_xyz": str(output_dir / "initial.xyz") if output_dir else None,
            "final_xyz": str(output_dir / "final.xyz") if output_dir else None
        }

    def prepare_endpoints_only(
        self,
        num_carboxyl: int = 0,
        num_hydroxyl: int = 0,
        surface_size: Tuple[int, int] = (4, 4),
        start_distance: float = 5.0,
        end_distance: float = 2.5,
        output_dir: Optional[Union[str, Path]] = None,
        constrain_endpoints: bool = True,
        generate_orca_inputs: bool = True,
    ) -> Dict[str, str]:
        """
        Prepare endpoint structures (and optionally ORCA inputs).

        Generates constrained endpoint optimizations: surface frozen, Zn-surface
        distance fixed, functional groups relaxed. No neb.inp file is created.

        Args:
            num_carboxyl: Number of carboxyl groups
            num_hydroxyl: Number of hydroxyl groups
            surface_size: Surface supercell size (nx, ny)
            start_distance: Initial Zn²⁺ distance from surface (Å)
            end_distance: Final Zn²⁺ distance from surface (Å)
            output_dir: Directory for output files
            constrain_endpoints: If True, freeze carbon surface and fix Zn-surface
                distance during endpoint optimization. Default True.
            generate_orca_inputs: If True, generate initial.inp and final.inp.
                If False (e.g. for UFF mode), only export initial.xyz and final.xyz.

        Returns:
            Dictionary with initial_xyz, final_xyz, and optionally initial_input,
            final_input when generate_orca_inputs is True.
        """
        if output_dir:
            output_dir = Path(output_dir)
            output_dir.mkdir(parents=True, exist_ok=True)

        # 1-8: Same as prepare_neb_calculation, but skip neb.inp
        pristine = self.surface_builder.build_pristine_surface(supercell_size=surface_size)
        c_coords = np.array([s.coords for s in pristine if s.specie.symbol == "C"])
        zn_center = (float(np.mean(c_coords[:, 0])), float(np.mean(c_coords[:, 1])))
        surface_structure, fg_bonds = self.surface_builder.add_functional_groups(
            pristine,
            num_carboxyl=num_carboxyl,
            num_hydroxyl=num_hydroxyl,
            zn_center_xy=zn_center,
        )
        initial_structure = self.create_initial_structure(
            surface_structure.copy(),
            start_distance
        )
        final_structure = self.create_final_structure(
            surface_structure.copy(),
            end_distance
        )
        initial_structure = self._ensure_atom_ordering(initial_structure)
        final_structure = self._ensure_atom_ordering(final_structure)
        initial_ase = self._pmg_to_ase(initial_structure)
        final_ase = self._pmg_to_ase(final_structure)

        n_graphene = len(pristine)
        zn_idx = len(initial_ase) - 1
        # Freeze carbon surface AND Zn ion so only functional groups optimize.
        # This keeps Zn height and surface fixed; functional groups relax to avoid steric hindrance.
        freeze_atoms = (
            list(range(n_graphene)) + [zn_idx] if constrain_endpoints else None
        )
        # Constrain FG-surface bonds so the attachment doesn't stretch/detach
        fg_bond_constraints = (
            [(i, j, d) for i, j, d in fg_bonds] if constrain_endpoints else []
        )

        if output_dir:
            self._export_xyz(initial_ase, output_dir / "initial.xyz")
            self._export_xyz(final_ase, output_dir / "final.xyz")

        result = {
            "initial_xyz": str(output_dir / "initial.xyz") if output_dir else None,
            "final_xyz": str(output_dir / "final.xyz") if output_dir else None,
            "n_graphene": n_graphene,
            "initial_ase": initial_ase,
            "final_ase": final_ase,
            "freeze_atoms": freeze_atoms,
            "fg_bond_constraints": fg_bond_constraints,
            "initial_bond_constraint": None,
            "final_bond_constraint": None,
        }
        if generate_orca_inputs:
            charge = 2 - num_carboxyl
            self.orca_gen.charge = charge
            self.orca_gen.multiplicity = self._calculate_multiplicity(
                initial_ase, charge=charge
            )
            initial_input = self.orca_gen.generate_from_ase(
                initial_ase,
                calc_type="opt",
                output_file=output_dir / "initial.inp" if output_dir else None,
                title="Initial endpoint (Zn²⁺ at start distance)",
                atom1_idx=None,
                atom2_idx=None,
                distance=None,
                freeze_atoms=freeze_atoms,
                bond_constraints=fg_bond_constraints,
            )
            final_input = self.orca_gen.generate_from_ase(
                final_ase,
                calc_type="opt",
                output_file=output_dir / "final.inp" if output_dir else None,
                title="Final endpoint (Zn²⁺ at end distance)",
                atom1_idx=None,
                atom2_idx=None,
                distance=None,
                freeze_atoms=freeze_atoms,
                bond_constraints=fg_bond_constraints,
            )
            result["initial_input"] = initial_input
            result["final_input"] = final_input
        return result

    def prepare_scan_base(
        self,
        num_carboxyl: int = 0,
        num_hydroxyl: int = 0,
        surface_size: Tuple[int, int] = (4, 4),
    ) -> Dict:
        """
        Build base surface with functional groups (no Zn) for distance scan.

        Returns surface structure and metadata needed to create structures
        at arbitrary Zn-surface distances.

        Returns:
            Dict with:
                - surface_structure: pymatgen Structure (FGs, no Zn)
                - n_graphene: number of carbon atoms
                - fg_bond_constraints: list of (i, j, d) for FG-surface bonds
        """
        pristine = self.surface_builder.build_pristine_surface(
            supercell_size=surface_size
        )
        c_coords = np.array([s.coords for s in pristine if s.specie.symbol == "C"])
        zn_center = (float(np.mean(c_coords[:, 0])), float(np.mean(c_coords[:, 1])))
        surface_structure, fg_bonds = self.surface_builder.add_functional_groups(
            pristine,
            num_carboxyl=num_carboxyl,
            num_hydroxyl=num_hydroxyl,
            zn_center_xy=zn_center,
        )
        n_graphene = len(pristine)
        fg_bond_constraints = [(i, j, d) for i, j, d in fg_bonds]
        return {
            "surface_structure": surface_structure,
            "n_graphene": n_graphene,
            "fg_bond_constraints": fg_bond_constraints,
        }

    def create_initial_structure(
        self,
        surface_structure,
        distance: float,
        position: Optional[Tuple[float, float]] = None
    ):
        """Create initial structure with Zn²⁺ at start distance."""
        return self.surface_builder.add_zn2_ion(
            surface_structure,
            distance=distance,
            position=position
        )

    def create_final_structure(
        self,
        surface_structure,
        distance: float,
        position: Optional[Tuple[float, float]] = None
    ):
        """Create final structure with Zn²⁺ at end distance."""
        return self.surface_builder.add_zn2_ion(
            surface_structure,
            distance=distance,
            position=position
        )

    def _ensure_atom_ordering(self, structure):
        """
        Ensure Zn²⁺ is the last atom in the structure (critical for NEB).
        """
        # Find Zn²⁺ index
        zn_indices = [i for i, site in enumerate(structure) if site.specie.symbol == "Zn"]
        
        if not zn_indices:
            return structure  # No Zn²⁺ found
        
        zn_idx = zn_indices[0]
        
        # If Zn²⁺ is already last, return as-is
        if zn_idx == len(structure) - 1:
            return structure
        
        # Reorder: move Zn²⁺ to end
        sites = list(structure)
        zn_site = sites.pop(zn_idx)
        sites.append(zn_site)
        
        # Create new structure with reordered sites
        if PYMATGEN_AVAILABLE:
            new_structure = Structure(
                structure.lattice,
                [site.specie for site in sites],
                [site.coords for site in sites],
                coords_are_cartesian=True
            )
        else:
            # Fallback: use structure's class
            new_structure = structure.__class__(
                structure.lattice,
                [site.specie for site in sites],
                [site.coords for site in sites],
                coords_are_cartesian=True
            )
        
        return new_structure

    def _calculate_multiplicity(self, atoms, charge: int) -> int:
        """
        Calculate correct multiplicity based on electron count.
        """
        total_electrons = sum(atoms.get_atomic_numbers())
        net_electrons = total_electrons - charge
        
        if net_electrons % 2 == 0:
            return 1  # Even → singlet
        else:
            return 2  # Odd → doublet
    
    def _pmg_to_ase(self, structure):
        """Convert pymatgen structure to ASE Atoms."""
        symbols = [site.specie.symbol for site in structure]
        positions = structure.cart_coords
        cell = structure.lattice.matrix
        return Atoms(symbols=symbols, positions=positions, cell=cell, pbc=True)
    
    def _export_xyz(self, atoms: Atoms, output_file: Path):
        """Export ASE Atoms to XYZ file format."""
        if not ASE_AVAILABLE:
            raise ImportError("ASE is required for XYZ export. Install with: pip install ase")
        
        write(str(output_file), atoms, format='xyz')
