"""
ORCA input file generator using ASE and pymatgen structures.

This module generates ORCA input files from ASE Atoms objects or pymatgen
Structures, with support for implicit solvent models (CPCM, SMD) and
distance constraints for geometry optimization.
"""

from __future__ import annotations

from typing import Optional, List, Dict, Union, Tuple
from pathlib import Path
import numpy as np

try:
    from ase import Atoms
    ASE_AVAILABLE = True
except ImportError:
    ASE_AVAILABLE = False
    Atoms = None  # type: ignore[misc, assignment]

try:
    from pymatgen.core import Structure, Molecule
    PYMATGEN_AVAILABLE = True
except ImportError:
    PYMATGEN_AVAILABLE = False


class OrcaInputGenerator:
    """
    Generate ORCA input files from molecular/crystal structures.
    
    Supports implicit solvent models (CPCM, SMD) and distance constraints.
    """
    
    def __init__(
        self,
        method: str = "B3LYP",
        basis_set: str = "def2-TZVP",
        charge: int = 0,
        multiplicity: int = 1,
        solvent: Optional[str] = None,
        solvation_model: str = "CPCM",
        nprocs: int = 4,
        memory: str = "16GB",  # Increased default for large systems
        dispersion: Optional[str] = "D3BJ",
        scf_convergence: str = "TightSCF",
        grid: str = "DefGrid4",  # ORCA 6 uses DefGrid4 instead of Grid4
        additional_keywords: Optional[List[str]] = None
    ):
        self.method = method
        self.basis_set = basis_set
        self.charge = charge
        self.multiplicity = multiplicity
        self.solvent = solvent
        self.solvation_model = solvation_model
        self.nprocs = nprocs
        self.memory = memory
        self.dispersion = dispersion
        self.scf_convergence = scf_convergence
        self.grid = grid
        self.additional_keywords = additional_keywords or []
        
        if solvation_model not in ["CPCM", "SMD"]:
            raise ValueError(f"solvation_model must be 'CPCM' or 'SMD', got {solvation_model}")
    
    def generate_from_ase(
        self,
        atoms: Union[Atoms, np.ndarray],  # Can also take symbols+positions later if needed
        calc_type: str = "opt",
        output_file: Optional[Union[str, Path]] = None,
        title: str = "ORCA Calculation",
        atom1_idx: Optional[int] = None,
        atom2_idx: Optional[int] = None,
        distance: Optional[float] = None,
        freeze_atoms: Optional[List[int]] = None,
        bond_constraints: Optional[List[Tuple[int, int, float]]] = None,
    ) -> str:
        """
        Generate ORCA input from ASE Atoms object.
        """
        if not ASE_AVAILABLE:
            raise ImportError("ASE is required. Install with: pip install ase")
        
        command = self._build_command(calc_type)
        input_content = self._build_input_file(
            atoms, command, title, atom1_idx, atom2_idx, distance, freeze_atoms,
            bond_constraints=bond_constraints,
        )
        
        if output_file:
            output_path = Path(output_file)
            output_path.parent.mkdir(parents=True, exist_ok=True)
            with open(output_path, 'w') as f:
                f.write(input_content)
        
        return input_content

    def _build_command(self, calc_type: str) -> str:
        parts = [f"!{self.method} {self.basis_set}"]
        
        if self.dispersion:
            parts.append(self.dispersion)
        
        if self.solvent:
            if self.solvation_model == "CPCM":
                parts.append(f"CPCM({self.solvent})")
            elif self.solvation_model == "SMD":
                parts.append(f"SMD(solvent={self.solvent})")
        
        if calc_type == "opt":
            parts.append("OPT")
        elif calc_type == "freq":
            parts.append("FREQ")
        elif calc_type == "ts":
            parts.append("OPT TS")
        
        if self.scf_convergence:
            parts.append(self.scf_convergence)
        
        # Grid settings are handled in %method block for ORCA 6
        # Don't add grid to simple input line
        
        parts.extend(self.additional_keywords)
        
        return " ".join(parts)

    def _build_input_file(
        self,
        atoms: Atoms,
        command: str,
        title: str,
        atom1_idx: Optional[int] = None,
        atom2_idx: Optional[int] = None,
        distance: Optional[float] = None,
        freeze_atoms: Optional[List[int]] = None,
        bond_constraints: Optional[List[Tuple[int, int, float]]] = None,
    ) -> str:
        lines = [
            command,
            "",
            f"%pal nprocs {self.nprocs} end",
            f"%maxcore {self._parse_memory()}",
        ]
        
        # Add grid settings in %method block (ORCA 6 format)
        if self.grid:
            # Extract grid number from grid string (e.g., "Grid4" -> 4, "DefGrid4" -> 4)
            grid_num = "4"  # default
            if "Grid" in self.grid or "DefGrid" in self.grid:
                import re
                match = re.search(r'(\d+)', self.grid)
                if match:
                    grid_num = match.group(1)
            lines.extend([
                "",
                "%method",
                f"  Grid {grid_num}",
                "end",
            ])
        
        # Note: CPCM solvent is already specified in the simple input line as CPCM(solvent)
        # Only add %cpcm block if we need additional CPCM settings (not needed for basic usage)
        # For SMD, we might need a %smd block, but SMD(solvent=...) in simple line should work
            
        geom_block = self._build_geom_constraints(
            atoms, atom1_idx, atom2_idx, distance, freeze_atoms,
            bond_constraints=bond_constraints,
        )
        if geom_block:
            lines.extend(["", geom_block])
        
        lines.extend([
            "",
            f"* xyz {self.charge} {self.multiplicity}",
        ])
        
        symbols = atoms.get_chemical_symbols()
        positions = atoms.get_positions()
        
        for symbol, pos in zip(symbols, positions):
            lines.append(
                f"  {symbol:4s} {pos[0]:15.10f} {pos[1]:15.10f} {pos[2]:15.10f}"
            )
        
        lines.append("*")
        
        return "\n".join(lines)

    def add_distance_constraint(
        self,
        atoms: Atoms,
        atom1_idx: int,
        atom2_idx: int,
        distance: float
    ) -> str:
        """
        Generate ORCA distance constraint block.
        Uses unit matrix Hessian (InHess 1) for stability with constraints.
        """
        # ORCA uses 0-based indexing for constraints
        n1 = atom1_idx
        n2 = atom2_idx
        
        constraint_block = f"""%geom
 InHess Unit
 Constraints
  {{ B {n1} {n2} {distance:.6f} C }}
 end
end"""
        return constraint_block

    def _build_geom_constraints(
        self,
        atoms: Atoms,
        atom1_idx: Optional[int],
        atom2_idx: Optional[int],
        distance: Optional[float],
        freeze_atoms: Optional[List[int]],
        bond_constraints: Optional[List[Tuple[int, int, float]]] = None,
    ) -> str:
        """
        Build %geom block with optional distance constraint and/or atom freezing.
        Freezes atoms (carbon surface), optionally constrains Zn-surface distance,
        and constrains FG-surface bonds so attachment points don't detach.
        Uses InHess Unit for stability with constraints.
        """
        # ORCA uses 0-based indexing for constraints
        constraints: List[str] = []
        if freeze_atoms:
            for idx in sorted(set(freeze_atoms)):
                if 0 <= idx < len(atoms):
                    constraints.append(f"  {{ C {idx} C }}")
        if atom1_idx is not None and atom2_idx is not None and distance is not None:
            n1, n2 = atom1_idx, atom2_idx
            constraints.append(f"  {{ B {n1} {n2} {distance:.6f} C }}")
        if bond_constraints:
            for i, j, d in bond_constraints:
                if 0 <= i < len(atoms) and 0 <= j < len(atoms):
                    constraints.append(f"  {{ B {i} {j} {d:.6f} C }}")
        if not constraints:
            return ""
        constraint_lines = "\n".join(constraints)
        return f"""%geom
 InHess Unit
 Constraints
{constraint_lines}
 end
end"""

    def generate_neb_input(
        self,
        initial_atoms: Optional[Atoms] = None,
        final_atoms: Optional[Atoms] = None,
        output_file: Optional[Union[str, Path]] = None,
        neb_images: Optional[int] = None,
        initial_xyz_file: Optional[Union[str, Path]] = None,
        final_xyz_file: Optional[Union[str, Path]] = None
    ) -> str:
        """
        Generate ORCA NEB-TS input file.
        
        Args:
            initial_atoms: Initial structure (reactant) as ASE Atoms (optional if initial_xyz_file provided)
            final_atoms: Final structure (product) as ASE Atoms (optional if final_xyz_file provided)
            output_file: Path to save input file
            neb_images: Number of NEB images (optional, ORCA will auto-determine if not specified)
            initial_xyz_file: Path to initial XYZ file (preferred - uses optimized structure)
            final_xyz_file: Path to final XYZ file (preferred - uses optimized structure)
        
        Returns:
            Input file content as string
        """
        if not ASE_AVAILABLE:
            raise ImportError("ASE is required. Install with: pip install ase")
        
        # Determine XYZ file paths
        if output_file:
            output_path = Path(output_file)
            if initial_xyz_file:
                initial_xyz_path = Path(initial_xyz_file)
            else:
                initial_xyz_path = output_path.parent / "initial.xyz"
            
            if final_xyz_file:
                final_xyz_path = Path(final_xyz_file)
            else:
                final_xyz_path = output_path.parent / "final.xyz"
            
            # Use relative paths if in same directory
            if initial_xyz_path.parent == output_path.parent:
                initial_xyz_ref = initial_xyz_path.name
            else:
                initial_xyz_ref = str(initial_xyz_path)
            
            if final_xyz_path.parent == output_path.parent:
                final_xyz_ref = final_xyz_path.name
            else:
                final_xyz_ref = str(final_xyz_path)
        else:
            initial_xyz_ref = initial_xyz_file if initial_xyz_file else "initial.xyz"
            final_xyz_ref = final_xyz_file if final_xyz_file else "final.xyz"
        
        # If XYZ files are provided, use them; otherwise verify atom ordering consistency
        if initial_atoms is not None and final_atoms is not None:
            initial_symbols = initial_atoms.get_chemical_symbols()
            final_symbols = final_atoms.get_chemical_symbols()
            
            if len(initial_symbols) != len(final_symbols):
                raise ValueError(f"Initial and final structures must have same number of atoms. "
                               f"Got {len(initial_symbols)} vs {len(final_symbols)}")
            
            if initial_symbols != final_symbols:
                raise ValueError("Initial and final structures must have identical atom ordering")
        
        # Build command line with NEB-TS keyword
        command_parts = [f"!{self.method} {self.basis_set}"]
        
        if self.dispersion:
            command_parts.append(self.dispersion)
        
        if self.solvent:
            if self.solvation_model == "CPCM":
                command_parts.append(f"CPCM({self.solvent})")
            elif self.solvation_model == "SMD":
                command_parts.append(f"SMD(solvent={self.solvent})")
        
        command_parts.append("NEB-TS")
        
        if self.scf_convergence:
            command_parts.append(self.scf_convergence)
        
        command = " ".join(command_parts)
        
        # Build input file
        lines = [
            command,
            "",
            f"%pal nprocs {self.nprocs} end",
            f"%maxcore {self._parse_memory()}",
        ]
        
        # Add grid settings
        if self.grid:
            grid_num = "4"
            if "Grid" in self.grid or "DefGrid" in self.grid:
                import re
                match = re.search(r'(\d+)', self.grid)
                if match:
                    grid_num = match.group(1)
            lines.extend([
                "",
                "%method",
                f"  Grid {grid_num}",
                "end",
            ])
        
        # NEB block
        neb_block = [
            "",
            "%NEB",
            f'  NEB_END_XYZFILE "{final_xyz_ref}"',
            "  PREOPT_ENDS TRUE",
            "END"
        ]
        lines.extend(neb_block)
        
        # Use XYZFILE to read initial structure from file (preferred for optimized structures)
        lines.extend([
            "",
            f"* XYZFILE {self.charge} {self.multiplicity} {initial_xyz_ref}",
        ])
        
        input_content = "\n".join(lines)
        if not input_content.endswith("\n"):
            input_content += "\n"
        
        # Save input file
        if output_file:
            output_path = Path(output_file)
            output_path.parent.mkdir(parents=True, exist_ok=True)
            with open(output_path, 'w', encoding='utf-8', newline='\n') as f:
                f.write(input_content)
            
            # Export structures to XYZ files if provided and files don't exist
            if initial_atoms is not None:
                initial_xyz_path = output_path.parent / "initial.xyz"
                if not initial_xyz_path.exists():
                    self._export_xyz(initial_atoms, initial_xyz_path)
            
            if final_atoms is not None:
                final_xyz_path = output_path.parent / "final.xyz"
                if not final_xyz_path.exists():
                    self._export_xyz(final_atoms, final_xyz_path)
        
        return input_content
    
    def _export_xyz(self, atoms: Atoms, output_file: Path):
        """Export ASE Atoms to XYZ file format."""
        if not ASE_AVAILABLE:
            raise ImportError("ASE is required for XYZ export. Install with: pip install ase")
        
        try:
            from ase.io import write
            write(str(output_file), atoms, format='xyz')
        except ImportError:
            # Fallback: manual XYZ format
            symbols = atoms.get_chemical_symbols()
            positions = atoms.get_positions()
            
            with open(output_file, 'w') as f:
                f.write(f"{len(symbols)}\n")
                f.write("Generated for ORCA NEB\n")
                for symbol, pos in zip(symbols, positions):
                    f.write(f"{symbol:4s} {pos[0]:15.10f} {pos[1]:15.10f} {pos[2]:15.10f}\n")

    def _parse_memory(self) -> int:
        memory_str = self.memory.upper().replace(" ", "")
        if "GB" in memory_str:
            return int(float(memory_str.replace("GB", "")) * 1000)
        elif "MB" in memory_str:
            return int(memory_str.replace("MB", ""))
        else:
            return 8000
