"""
Nudged Elastic Band (NEB) calculator for Zn2+ adsorption pathways.

This module handles the complete workflow for NEB calculations including:
- Building initial and final structures with Zn²⁺ at different distances
- Constrained pre-optimization of endpoints
- NEB-TS input generation
- Path analysis
"""

from typing import Optional, Dict, Tuple, List, Union
from pathlib import Path
import numpy as np

try:
    from ase import Atoms
    from ase.io import read, write
    ASE_AVAILABLE = True
except ImportError:
    ASE_AVAILABLE = False

try:
    from pymatgen.core import Structure
    PYMATGEN_AVAILABLE = True
except ImportError:
    PYMATGEN_AVAILABLE = False

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
        output_dir: Optional[Union[str, Path]] = None
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
        
        Returns:
            Dictionary with paths to generated input files
        """
        if output_dir:
            output_dir = Path(output_dir)
            output_dir.mkdir(parents=True, exist_ok=True)
        
        # 1. Build functionalized surface (common to both endpoints)
        pristine = self.surface_builder.build_pristine_surface(supercell_size=surface_size)
        surface_structure = self.surface_builder.add_functional_groups(
            pristine,
            num_carboxyl=num_carboxyl,
            num_hydroxyl=num_hydroxyl
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
        
        # 6. Generate unconstrained optimization inputs for endpoints
        # According to ORCA documentation, endpoints should be optimized separately
        # without constraints. PREOPT_ENDS TRUE in NEB will handle re-optimization.
        # Set charge and multiplicity
        self.orca_gen.charge = 2
        self.orca_gen.multiplicity = self._calculate_multiplicity(initial_ase, charge=2)
        
        # Initial endpoint (unconstrained optimization)
        # This allows the structure to relax naturally while maintaining approximate distance
        initial_input = self.orca_gen.generate_from_ase(
            initial_ase,
            calc_type="opt",
            output_file=output_dir / "initial.inp" if output_dir else None,
            title="Initial endpoint (Zn²⁺ at start distance)",
            atom1_idx=None,  # No constraints
            atom2_idx=None,
            distance=None
        )
        
        # Final endpoint (unconstrained optimization)
        final_input = self.orca_gen.generate_from_ase(
            final_ase,
            calc_type="opt",
            output_file=output_dir / "final.inp" if output_dir else None,
            title="Final endpoint (Zn²⁺ at end distance)",
            atom1_idx=None,  # No constraints
            atom2_idx=None,
            distance=None
        )
        
        # 8. Export structures as XYZ files for NEB
        if output_dir:
            self._export_xyz(initial_ase, output_dir / "initial.xyz")
            self._export_xyz(final_ase, output_dir / "final.xyz")
        
        # 9. Generate NEB input file
        neb_input = self.orca_gen.generate_neb_input(
            initial_ase,
            final_ase,
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
