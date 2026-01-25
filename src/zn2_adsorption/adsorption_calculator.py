"""
Calculate adsorption energies for Zn2+ on functionalized carbon surfaces.

This module handles the complete workflow for adsorption energy calculations
including surface building, ORCA input generation with constraints,
and energy evaluation.
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

from zn2_adsorption.orca_generator import OrcaInputGenerator
from zn2_adsorption.surface_builder import FunctionalizedGrapheneBuilder


class AdsorptionEnergyCalculator:
    """
    Orchestrate the Zn2+ adsorption energy calculation workflow.
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

    def prepare_calculations(
        self,
        num_carboxyl: int = 0,
        num_hydroxyl: int = 0,
        surface_size: Tuple[int, int] = (4, 4),
        zn2_distance: float = 3.0,
        constrain_distance: bool = True,
        output_dir: Optional[Union[str, Path]] = None
    ) -> Dict[str, str]:
        """
        Prepare ORCA input files for the adsorption energy calculation.
        """
        if output_dir:
            output_dir = Path(output_dir)
            output_dir.mkdir(parents=True, exist_ok=True)
            
        # 1. Build functionalized surface
        pristine = self.surface_builder.build_pristine_surface(supercell_size=surface_size)
        surface_structure = self.surface_builder.add_functional_groups(
            pristine,
            num_carboxyl=num_carboxyl,
            num_hydroxyl=num_hydroxyl
        )
        
        # 2. Add Zn2+ ion
        zn2_surface_structure = self.surface_builder.add_zn2_ion(
            surface_structure,
            distance=zn2_distance
        )
        
        # 3. Find nearest atom for constraint if needed
        zn2_idx = len(zn2_surface_structure) - 1 # Last atom added
        zn2_pos = zn2_surface_structure[zn2_idx].coords
        
        nearest_idx, actual_dist = self.surface_builder.find_nearest_surface_atom(
            surface_structure, # Look in original surface
            zn2_pos
        )
        
        # 4. Generate ORCA inputs
        
        # Convert to ASE for OrcaInputGenerator
        # (Assuming Structure.to_ase_atoms() is available or manual conversion)
        zn2_surface_ase = self._pmg_to_ase(zn2_surface_structure)
        surface_ase = self._pmg_to_ase(surface_structure)
        zn2_only_ase = Atoms(symbols=['Zn'], positions=[[0, 0, 0]])
        
        # Zn2+@surface (charge +2)
        self.orca_gen.charge = 2
        self.orca_gen.multiplicity = self._calculate_multiplicity(zn2_surface_ase, charge=2)
        zn2_surface_input = self.orca_gen.generate_from_ase(
            zn2_surface_ase,
            calc_type="opt",
            output_file=output_dir / "zn2_surface.inp" if output_dir else None,
            title="Zn2+ on functionalized graphene",
            atom1_idx=zn2_idx if constrain_distance else None,
            atom2_idx=nearest_idx if constrain_distance else None,
            distance=actual_dist if constrain_distance else None
        )
        
        # Surface only (charge 0)
        self.orca_gen.charge = 0
        self.orca_gen.multiplicity = self._calculate_multiplicity(surface_ase, charge=0)
        surface_input = self.orca_gen.generate_from_ase(
            surface_ase,
            calc_type="opt",
            output_file=output_dir / "surface.inp" if output_dir else None,
            title="Functionalized graphene surface"
        )
        
        # Zn2+ only (charge +2)
        self.orca_gen.charge = 2
        self.orca_gen.multiplicity = self._calculate_multiplicity(zn2_only_ase, charge=2)
        zn2_input = self.orca_gen.generate_from_ase(
            zn2_only_ase,
            calc_type="sp", # Single point for isolated ion
            output_file=output_dir / "zn2_ion.inp" if output_dir else None,
            title="Isolated Zn2+ ion"
        )
        
        return {
            "zn2_surface_input": zn2_surface_input,
            "surface_input": surface_input,
            "zn2_input": zn2_input
        }

    def _calculate_multiplicity(self, atoms, charge: int) -> int:
        """
        Calculate correct multiplicity based on electron count.
        
        Rules:
        - If (total_electrons - charge) is even → multiplicity 1 (singlet, closed shell)
        - If (total_electrons - charge) is odd → multiplicity 2 (doublet, one unpaired electron)
        """
        # Get atomic numbers from ASE Atoms object
        total_electrons = sum(atoms.get_atomic_numbers())  # Sum of atomic numbers = total electrons
        net_electrons = total_electrons - charge
        
        if net_electrons % 2 == 0:
            return 1  # Even number of electrons → singlet (closed shell)
        else:
            return 2  # Odd number of electrons → doublet (one unpaired electron)
    
    def _pmg_to_ase(self, structure):
        """Convert pymatgen structure to ASE Atoms."""
        symbols = [site.specie.symbol for site in structure]
        positions = structure.cart_coords
        cell = structure.lattice.matrix
        return Atoms(symbols=symbols, positions=positions, cell=cell, pbc=True)
