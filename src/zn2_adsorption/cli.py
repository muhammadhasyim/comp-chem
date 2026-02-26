import argparse
import copy
import json
import os
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import List, Optional, Tuple, Dict, Any
from pathlib import Path

from zn2_adsorption.neb_calculator import NebCalculator
from zn2_adsorption.orca_generator import OrcaInputGenerator

def calculate_neb_images(start_dist: float, end_dist: float) -> int:
    """
    Calculate number of NEB images based on distance difference.
    
    Args:
        start_dist: Initial Zn²⁺ distance from surface (Å)
        end_dist: Final Zn²⁺ distance from surface (Å)
    
    Returns:
        Number of NEB images (minimum 5, maximum 20)
    """
    distance_diff = abs(end_dist - start_dist)
    # ~1 image per 0.5 Å of distance change, minimum 5, maximum 20
    images = max(5, min(20, int(distance_diff / 0.5) + 1))
    return images

def parse_arguments(args: List[str]) -> argparse.Namespace:
    """
    Parse command-line arguments for the adsorption calculator.
    """
    parser = argparse.ArgumentParser(
        description="Zn2+ MEP pipeline: FSM (Freezing String) or NEB for adsorption path"
    )
    
    # Surface parameters
    parser.add_argument(
        "--carboxyl", "-c",
        type=int,
        default=0,
        help="Number of carboxylate (-COO⁻) groups (default: 0)"
    )
    parser.add_argument(
        "--hydroxyl", "-y",
        type=int,
        default=0,
        help="Number of hydroxyl (-OH) groups (default: 0)"
    )
    parser.add_argument(
        "--size", "-s",
        type=str,
        default="4x4",
        help="Surface size as 'NxM' or single integer for square (default: '4x4')"
    )
    
    # NEB parameters
    parser.add_argument(
        "--start-distance",
        type=float,
        required=True,
        help="Initial Zn2+ distance from surface in Angstroms (required)"
    )
    parser.add_argument(
        "--end-distance",
        type=float,
        required=True,
        help="Final Zn2+ distance from surface in Angstroms (required)"
    )
    parser.add_argument(
        "--neb-images",
        type=int,
        default=None,
        help="Number of NEB images (default: auto-calculated based on distance difference)"
    )
    
    # ORCA parameters
    parser.add_argument(
        "--method",
        type=str,
        default="B3LYP",
        help="DFT method (default: 'B3LYP')"
    )
    parser.add_argument(
        "--basis",
        type=str,
        default="def2-TZVP",
        help="Basis set (default: 'def2-TZVP')"
    )
    parser.add_argument(
        "--solvent",
        type=str,
        default="water",
        help="Solvent name (default: 'water')"
    )
    parser.add_argument(
        "--solvation",
        type=str,
        default="CPCM",
        choices=["CPCM", "SMD"],
        help="Solvation model 'CPCM' or 'SMD' (default: 'CPCM')"
    )
    parser.add_argument(
        "--memory",
        type=str,
        default="16GB",
        help="Memory allocation for ORCA (e.g., '16GB', '32GB') (default: '16GB')"
    )
    parser.add_argument(
        "--nprocs",
        type=int,
        default=4,
        help="Number of processors for ORCA (default: 4)"
    )
    parser.add_argument(
        "--no-mpi",
        action="store_true",
        default=False,
        help="Disable MPI/parallel execution - run ORCA in serial mode (nprocs=1) (default: False)"
    )
    parser.add_argument(
        "--scf-convergence",
        type=str,
        default="TightSCF",
        choices=["TightSCF", "NormalSCF", "LooseSCF"],
        help="SCF convergence criteria (default: 'TightSCF')"
    )
    parser.add_argument(
        "--grid",
        type=int,
        default=4,
        choices=[1, 2, 3, 4, 5, 6, 7],
        help="Integration grid quality (1-7, higher is better but slower) (default: 4)"
    )
    parser.add_argument(
        "--chg",
        type=int,
        default=None,
        help="Total charge. If omitted: 2 - n_carboxyl (Zn²⁺ + carboxylate -1 each)"
    )
    parser.add_argument(
        "--mult",
        type=int,
        default=1,
        help="Spin multiplicity (default: 1)"
    )
    parser.add_argument(
        "--orca-command",
        type=str,
        default=None,
        help="Path to ORCA executable (default: auto-detect from repo or PATH)"
    )
    parser.add_argument(
        "--orca-simple",
        type=str,
        default=None,
        help="ORCA '!' line override (e.g. 'B3LYP def2-TZVP D3BJ CPCM(water)'). Default: built from --method, --basis, --solvent"
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        default=False,
        help="Verbose output (passed to FSM when --path-method fsm)"
    )
    
    # Endpoints-only mode (no pathway)
    parser.add_argument(
        "--endpoints-only",
        action="store_true",
        default=False,
        help="Generate and optionally run only initial and final endpoint optimizations (no NEB/FSM pathway). Uses UFF by default for fast runs."
    )
    parser.add_argument(
        "--calculator",
        type=str,
        default="uff",
        choices=["uff", "orca", "rowan"],
        help="Calculator: 'uff' (fast), 'orca' (DFT), or 'rowan' (cloud NNP). Default: uff."
    )
    parser.add_argument(
        "--model",
        type=str,
        default="uma_s_omol",
        help="[Rowan] Neural network potential model. Examples: "
             "uma_s_omol (UMA Small OMol25, default), uma_m_omol (UMA Medium), "
             "uma_s_omat, uma_m_omat, uma_s_omc, uma_m_omc, "
             "mace_mp_0, mace_mp_0b2_l, aimnet2_wb97md3, "
             "orb_v3_conservative_omol, orb_v3_conservative_inf_omat, "
             "omol25_conserving_s, egret_1, egret_1e, egret_1t, "
             "gfn2_xtb, gfn1_xtb, gfn0_xtb, gfn_ff, g_xtb.",
    )
    # Path method: NEB, FSM, or distance scan
    parser.add_argument(
        "--path-method",
        type=str,
        default="fsm",
        choices=["neb", "fsm", "scan"],
        help="Method: 'fsm' (Freezing String), 'neb' (Nudged Elastic Band), or 'scan' (distance scan) (default: 'fsm')"
    )
    # Scan-specific options (used when --path-method scan)
    parser.add_argument(
        "--scan-points",
        type=int,
        default=None,
        help="[SCAN] Number of distances to sample (default: auto from step 0.5 Å)"
    )
    parser.add_argument(
        "--scan-step",
        type=float,
        default=None,
        help="[SCAN] Step size in Angstrom; overrides --scan-points if set"
    )
    parser.add_argument(
        "--fg-sweep",
        action="store_true",
        default=False,
        help="[SCAN] Sweep over (carboxyl, hydroxyl) combinations from (0,0) to max; overrides --carboxyl and --hydroxyl",
    )
    parser.add_argument(
        "--fg-sweep-max-carboxyl",
        type=int,
        default=2,
        help="[SCAN] Max carboxylate count when --fg-sweep (default: 2)",
    )
    parser.add_argument(
        "--fg-sweep-max-hydroxyl",
        type=int,
        default=2,
        help="[SCAN] Max hydroxyl count when --fg-sweep (default: 2)",
    )
    # FSM-specific options (used when --path-method fsm)
    parser.add_argument(
        "--nnodes-min",
        type=int,
        default=10,
        help="[FSM] Minimum number of nodes on the string (default: 10)"
    )
    parser.add_argument(
        "--fsm-maxiter",
        type=int,
        default=2,
        help="[FSM] Max optimization iterations per node (default: 2)"
    )
    parser.add_argument(
        "--fsm-maxls",
        type=int,
        default=2,
        help="[FSM] Max line-search iterations (default: 2)"
    )
    parser.add_argument(
        "--fsm-suffix",
        type=str,
        default=None,
        help="[FSM] Suffix for FSM output directory (default: None)"
    )
    
    # Execution and output parameters
    parser.add_argument(
        "--output-dir", "-O",
        type=str,
        default="./orca_inputs",
        metavar="DIR",
        help="Output directory for ORCA/FSM files; fsm_reaction is created under DIR (default: './orca_inputs')"
    )
    parser.add_argument(
        "--run",
        action="store_true",
        default=False,
        help="Run optimization (UFF when --calculator uff, ORCA when --calculator orca). Default: False."
    )
    parser.add_argument(
        "--run-orca",
        action="store_true",
        default=False,
        help="(Alias for --run) Run the calculation. Kept for backward compatibility."
    )
    parser.add_argument(
        "--skip-endpoint-opt",
        action="store_true",
        default=False,
        help="[NEB] Skip pre-optimization of endpoints (let PREOPT_ENDS handle it in NEB). Faster but may converge slower. (default: False)"
    )
    parser.add_argument(
        "--no-constrain-endpoints",
        action="store_true",
        default=False,
        help="[NEB] Disable constrained endpoint optimization (freeze surface + fix Zn distance). By default, carbon surface is frozen and Zn-surface distance fixed while functional groups relax to avoid steric hindrance."
    )
    parser.add_argument(
        "--relax-surface",
        action="store_true",
        default=False,
        help="[Rowan scan] Allow graphene surface to relax during scan. "
             "Default: surface frozen. Angle constraints still enforce "
             "perpendicular Zn approach.",
    )
    parser.add_argument(
        "--freeze-functional-groups",
        action="store_true",
        default=False,
        help="[Scan] Fully freeze functional groups (no motion). "
             "Default: C-FG bonds constrained; FGs can flex/rotate.",
    )
    parser.add_argument(
        "--json-output",
        type=str,
        default="./mep_results.json",
        help="JSON output file path (default: './mep_results.json')"
    )
    parser.add_argument(
        "--plot",
        action="store_true",
        default=False,
        help="Plot NEB energy profile (saves path.png in output dir when NEB completes)"
    )
    
    parsed = parser.parse_args(args)
    parsed.run_orca = parsed.run_orca or parsed.run
    return parsed

def validate_inputs(args: argparse.Namespace) -> Tuple[Tuple[int, int], float, float]:
    """
    Validate user inputs and return normalized parameters.
    """
    # Parse size
    size_str = args.size.lower()
    if 'x' in size_str:
        try:
            nx, ny = map(int, size_str.split('x'))
            size = (nx, ny)
        except ValueError:
            raise ValueError(f"Invalid surface size format: {args.size}. Use 'NxM' or 'N'.")
    else:
        try:
            n = int(size_str)
            size = (n, n)
        except ValueError:
            raise ValueError(f"Invalid surface size format: {args.size}. Use 'NxM' or 'N'.")
            
    # Validate distances
    if args.start_distance <= 0:
        raise ValueError(f"Start distance must be positive, got {args.start_distance}")
    if args.end_distance <= 0:
        raise ValueError(f"End distance must be positive, got {args.end_distance}")
    if args.start_distance == args.end_distance:
        raise ValueError(f"Start and end distances must be different, both are {args.start_distance}")
        
    return size, args.start_distance, args.end_distance

def check_orca_available() -> Tuple[bool, Optional[str]]:
    """
    Check if ORCA is available in common locations:
    - orca/ folder
    - Folders containing "orca" in the name (case-insensitive)
    - Project root directories
    - PATH (avoiding system screen reader)
    
    Returns:
        Tuple of (is_available, path_to_orca_executable)
    """
    def find_orca_in_dir(directory: Path) -> Optional[str]:
        """Helper to find ORCA executable in a directory."""
        if not directory.exists() or not directory.is_dir():
            return None
        
        # Look for common ORCA executable names
        possible_names = ["orca", "orca_5", "orca_5_0_3", "orca_5_0_4", "orca_5_0_5", "orca_6"]
        for name in possible_names:
            orca_exe = directory / name
            if orca_exe.exists() and os.access(orca_exe, os.X_OK):
                return str(orca_exe)
        
        # If no standard name found, look for any executable file named "orca"
        for item in directory.iterdir():
            if item.is_file() and item.name == "orca" and os.access(item, os.X_OK):
                return str(item)
        
        return None
    
    # FIRST: Check for local orca/ folder in current working directory
    cwd = Path.cwd()
    result = find_orca_in_dir(cwd / "orca")
    if result:
        return (True, result)
    
    # SECOND: Check for ANY folder containing "orca" in name (case-insensitive) in current directory
    for item in cwd.iterdir():
        if item.is_dir() and "orca" in item.name.lower():
            result = find_orca_in_dir(item)
            if result:
                return (True, result)
    
    # THIRD: Check in project root (where setup.py/pyproject.toml is)
    # Try to find project root by looking for setup.py or pyproject.toml
    current = Path.cwd()
    for parent in [current] + list(current.parents):
        if (parent / "setup.py").exists() or (parent / "pyproject.toml").exists():
            # Check for "orca" folder
            result = find_orca_in_dir(parent / "orca")
            if result:
                return (True, result)
            
            # Check for ANY folder containing "orca" in name (case-insensitive) in project root
            try:
                for item in parent.iterdir():
                    if item.is_dir() and "orca" in item.name.lower():
                        result = find_orca_in_dir(item)
                        if result:
                            return (True, result)
            except PermissionError:
                # Skip directories we can't read
                pass
            break  # Found project root, no need to go further
    
    # FOURTH: Check parent directories up to 2 levels for orca folders
    # (in case project is in a subdirectory)
    for parent in list(current.parents)[:2]:
        try:
            for item in parent.iterdir():
                if item.is_dir() and "orca" in item.name.lower():
                    result = find_orca_in_dir(item)
                    if result:
                        return (True, result)
        except (PermissionError, OSError):
            # Skip directories we can't read
            continue
    
    # LAST: Check PATH for specific ORCA quantum chemistry executables
    # (Skip "orca" to avoid finding GNOME screen reader)
    orca_names = ["orca_5", "orca_5_0_3", "orca_5_0_4", "orca_5_0_5", "orca_6"]
    for name in orca_names:
        orca_path = shutil.which(name)
        if orca_path:
            return (True, orca_path)
    
    return (False, None)

def run_orca_command(input_file: Path, output_file: Path, orca_path: str, no_mpi: bool = False) -> bool:
    """
    Execute ORCA command.
    
    Args:
        input_file: Path to ORCA input file
        output_file: Path to ORCA output file
        orca_path: Path to ORCA executable
        no_mpi: If True, disable MPI execution by setting environment variables
    
    Returns:
        True if ORCA started successfully, False if there was an input error.
    Note: Long-running calculations may be interrupted; check output file for status.
    """
    print(f"  Running ORCA for {input_file.name}...")
    try:
        # Use absolute path for input file
        input_path = input_file.resolve()
        
        # Find ORCA directory and set up library path
        orca_exe = Path(orca_path).resolve()
        orca_dir = orca_exe.parent
        orca_lib_dir = orca_dir / "lib"
        
        # Set up environment with library path
        env = os.environ.copy()
        if orca_lib_dir.exists():
            # Add ORCA lib directory to LD_LIBRARY_PATH
            current_ld_path = env.get("LD_LIBRARY_PATH", "")
            if current_ld_path:
                env["LD_LIBRARY_PATH"] = f"{orca_lib_dir}:{current_ld_path}"
            else:
                env["LD_LIBRARY_PATH"] = str(orca_lib_dir)
        
        # If --no-mpi is set, configure environment to prevent MPI usage
        if no_mpi:
            # Force OpenMPI to use only 1 process
            env["OMPI_COMM_WORLD_SIZE"] = "1"
            env["OMPI_COMM_WORLD_RANK"] = "0"
            # Disable OpenMP parallelization as well
            env["OMP_NUM_THREADS"] = "1"
            # Prevent ORCA from detecting multiple processors
            env["OMPI_UNIVERSE_SIZE"] = "1"
            # Some ORCA installations check this
            env["ORCA_NPROCS"] = "1"
            # Try to find serial-only helper executables
            # Check for serial versions of helper executables (e.g., orca_prop instead of orca_prop_mpi)
            orca_bin_dir = orca_dir
            # If helper executables exist, ORCA will prefer serial versions when MPI is disabled
            # This is handled by ORCA's internal logic, but we ensure environment is set correctly
        
        # Change to input file's directory so relative paths in input file work correctly
        # ORCA interprets relative paths relative to its working directory, not the input file location
        input_dir = input_path.parent
        output_file_resolved = output_file.resolve()
        # If output file is in the same directory as input file, use it directly
        # Otherwise, place it in the input file's directory
        if output_file_resolved.parent == input_dir:
            output_path_in_dir = output_file_resolved
        else:
            output_path_in_dir = input_dir / output_file.name
        
        with open(output_path_in_dir, 'w') as out_f:
            result = subprocess.run(
                [str(orca_exe), str(input_path)],  # Use absolute path for input file
                stdout=out_f,
                stderr=subprocess.STDOUT,
                check=False,  # Don't raise on non-zero exit - check output file instead
                env=env,  # Use environment with library path
                cwd=str(input_dir)  # Run ORCA from the input file's directory so relative paths work
            )
        
        # Check output file for actual errors (not just exit code)
        # ORCA may return non-zero if interrupted, but still produce valid output
        # Use the output path in the input directory (where we actually wrote it)
        if output_path_in_dir.exists() and output_path_in_dir.stat().st_size > 0:
            try:
                # Try reading as text with error handling for binary data
                with open(output_path_in_dir, 'r', encoding='utf-8', errors='replace') as f:
                    content = f.read()
            except Exception as e:
                # Fallback: try reading as binary and decode with error handling
                try:
                    with open(output_path_in_dir, 'rb') as f:
                        content = f.read().decode('utf-8', errors='replace')
                except Exception:
                    # If we can't read it at all, assume it's still running or corrupted
                    print(f"  Warning: Could not read {output_path_in_dir.name} (may contain binary data). Assuming calculation is running.")
                    return True
            
            # Check for actual input errors
            if "INPUT ERROR" in content or "UNRECOGNIZED" in content:
                print(f"  Error: ORCA input error for {input_file.name}. Check {output_path_in_dir.name} for details.")
                return False
            # Check for segmentation fault
            if "segmentation fault" in content.lower() or "SIGSEGV" in content:
                print(f"  Error: ORCA crashed (segmentation fault) for {input_file.name}. Check {output_path_in_dir.name} for details.")
                return False
            # If ORCA started processing, consider it successful
            if "Starting time:" in content and ("Geometry Optimization" in content or "SCF" in content or "ORCA TERMINATED" in content or "NEB" in content):
                if result.returncode == 0:
                    print(f"  ✓ ORCA completed successfully for {input_file.name}")
                elif result.returncode == -11:
                    # SIGSEGV - segmentation fault
                    print(f"  ✗ ORCA crashed (segmentation fault) for {input_file.name}. Check {output_path_in_dir.name} for details.")
                    return False
                else:
                    print(f"  ⚠ ORCA was interrupted for {input_file.name} (exit code {result.returncode}). Check {output_path_in_dir.name} for progress.")
                return True
        
        # If we get here, something went wrong
        if result.returncode != 0:
            if result.returncode == -11:
                print(f"  ✗ ORCA crashed (segmentation fault) for {input_file.name}. Check {output_path_in_dir.name} for details.")
            else:
                print(f"  Error: ORCA failed for {input_file.name} (exit code {result.returncode}). Check {output_path_in_dir.name} for details.")
            return False
        
        return True
    except Exception as e:
        print(f"  Error running ORCA: {e}")
        return False

def extract_energy_from_orca(output_file: Path) -> Optional[float]:
    """
    Extract final SCF energy from ORCA output file.
    """
    if not output_file.exists():
        return None
        
    energy = None
    try:
        with open(output_file, 'r', encoding='utf-8', errors='replace') as f:
            for line in f:
                if "FINAL SINGLE POINT ENERGY" in line:
                    parts = line.split()
                    if len(parts) >= 5:
                        energy = float(parts[4])
    except Exception as e:
        print(f"  Error parsing {output_file.name}: {e}")
        
    return energy

def extract_optimized_geometry_from_orca(output_file: Path, xyz_file: Path) -> bool:
    """
    Extract the final optimized geometry from ORCA output and write to XYZ file.
    
    Args:
        output_file: Path to ORCA output file (.out)
        xyz_file: Path to write XYZ file
        
    Returns:
        True if successful, False otherwise
    """
    if not output_file.exists():
        return False
    
    try:
        with open(output_file, 'r', encoding='utf-8', errors='replace') as f:
            lines = f.readlines()
        
        # Find the last occurrence of "CARTESIAN COORDINATES (ANGSTROEM)"
        last_coord_start = -1
        for i in range(len(lines) - 1, -1, -1):
            if "CARTESIAN COORDINATES (ANGSTROEM)" in lines[i]:
                last_coord_start = i
                break
        
        if last_coord_start == -1:
            # Try to read from existing XYZ file if ORCA created one
            xyz_auto = output_file.with_suffix('.xyz')
            if xyz_auto.exists():
                import shutil
                shutil.copy2(xyz_auto, xyz_file)
                return True
            return False
        
        # Parse coordinates (skip separator line, read until blank line or next section)
        coords = []
        i = last_coord_start + 2  # Skip header and separator line
        while i < len(lines):
            line = lines[i].strip()
            if not line or line.startswith('-') or line.startswith('='):
                # Check if we've read all coordinates (next section starting)
                if coords and (line.startswith('CARTESIAN') or 'NO LB' in line or 'FINAL' in line.upper()):
                    break
                if not line:
                    i += 1
                    continue
                if line.startswith('-') and coords:
                    break
            
            # Parse coordinate line: "  C      0.019688    1.388033   15.004487"
            parts = line.split()
            if len(parts) >= 4:
                try:
                    symbol = parts[0]
                    x = float(parts[1])
                    y = float(parts[2])
                    z = float(parts[3])
                    coords.append((symbol, x, y, z))
                except (ValueError, IndexError):
                    pass
            
            i += 1
        
        if not coords:
            # Fallback: try to read from auto-generated XYZ file
            xyz_auto = output_file.with_suffix('.xyz')
            if xyz_auto.exists():
                import shutil
                shutil.copy2(xyz_auto, xyz_file)
                return True
            return False
        
        # Write XYZ file
        with open(xyz_file, 'w') as f:
            f.write(f"{len(coords)}\n")
            f.write(f"Optimized geometry from {output_file.name}\n")
            for symbol, x, y, z in coords:
                f.write(f"{symbol:4s} {x:15.10f} {y:15.10f} {z:15.10f}\n")
        
        return True
        
    except Exception as e:
        print(f"  Warning: Could not extract optimized geometry from {output_file.name}: {e}")
        # Fallback: try to read from auto-generated XYZ file
        xyz_auto = output_file.with_suffix('.xyz')
        if xyz_auto.exists():
            try:
                import shutil
                shutil.copy2(xyz_auto, xyz_file)
                return True
            except Exception:
                pass
        return False

def parse_neb_output(output_file: Path) -> Dict[str, Optional[float]]:
    """
    Parse NEB output file to extract transition state and path information.
    
    Returns:
        Dictionary with keys:
        - ts_energy: Transition state energy (Hartree)
        - initial_energy: Initial endpoint energy (Hartree)
        - final_energy: Final endpoint energy (Hartree)
        - barrier_height: Barrier height (E_TS - E_initial) in kcal/mol
        - reaction_energy: Reaction energy (E_final - E_initial) in kcal/mol
    """
    if not output_file.exists():
        return {
            "ts_energy": None,
            "initial_energy": None,
            "final_energy": None,
            "barrier_height": None,
            "reaction_energy": None
        }
    
    result = {
        "ts_energy": None,
        "initial_energy": None,
        "final_energy": None,
        "barrier_height": None,
        "reaction_energy": None
    }
    
    try:
        with open(output_file, 'r', encoding='utf-8', errors='replace') as f:
            content = f.read()
            lines = content.split('\n')
            
            # Look for transition state energy
            for i, line in enumerate(lines):
                if "TRANSITION STATE ENERGY" in line or "TS ENERGY" in line:
                    parts = line.split()
                    for j, part in enumerate(parts):
                        try:
                            energy = float(part)
                            result["ts_energy"] = energy
                            break
                        except ValueError:
                            continue
                
                # Look for initial endpoint energy (from PREOPT_ENDS)
                if "INITIAL ENDPOINT ENERGY" in line or "REACTANT ENERGY" in line:
                    parts = line.split()
                    for part in parts:
                        try:
                            energy = float(part)
                            result["initial_energy"] = energy
                            break
                        except ValueError:
                            continue
                
                # Look for final endpoint energy (from PREOPT_ENDS)
                if "FINAL ENDPOINT ENERGY" in line or "PRODUCT ENERGY" in line:
                    parts = line.split()
                    for part in parts:
                        try:
                            energy = float(part)
                            result["final_energy"] = energy
                            break
                        except ValueError:
                            continue
                
                # Also check for standard FINAL SINGLE POINT ENERGY patterns
                if "FINAL SINGLE POINT ENERGY" in line:
                    parts = line.split()
                    if len(parts) >= 5:
                        try:
                            energy = float(parts[4])
                            # If we haven't found TS energy yet, this might be it
                            # Or it could be an endpoint energy
                            if result["ts_energy"] is None:
                                result["ts_energy"] = energy
                        except ValueError:
                            pass
        
        # Calculate barrier height and reaction energy if we have the energies
        if result["ts_energy"] is not None and result["initial_energy"] is not None:
            barrier_hartree = result["ts_energy"] - result["initial_energy"]
            result["barrier_height"] = barrier_hartree * 627.509  # Convert to kcal/mol
        
        if result["initial_energy"] is not None and result["final_energy"] is not None:
            reaction_hartree = result["final_energy"] - result["initial_energy"]
            result["reaction_energy"] = reaction_hartree * 627.509  # Convert to kcal/mol
            
    except Exception as e:
        print(f"  Error parsing NEB output {output_file.name}: {e}")
    
    return result

def _prepare_fsm_reaction_dir(
    output_dir: Path,
    initial_xyz: Path,
    final_xyz: Path,
    charge: int = 2,
    multiplicity: int = 1,
) -> Path:
    """
    Create fsm_reaction directory with combined initial.xyz (two frames) and chg, mult.
    Returns path to fsm_reaction directory.
    """
    fsm_dir = output_dir / "fsm_reaction"
    fsm_dir.mkdir(parents=True, exist_ok=True)
    # Read initial and final XYZ (standard: line1 natoms, line2 comment, then coords)
    ini_lines = initial_xyz.read_text().strip().splitlines()
    fin_lines = final_xyz.read_text().strip().splitlines()
    if not ini_lines or not fin_lines:
        raise ValueError(f"Empty XYZ: {initial_xyz} or {final_xyz}")
    n_ini = int(ini_lines[0])
    n_fin = int(fin_lines[0])
    if n_ini != n_fin:
        raise ValueError(f"Atom count mismatch: initial {n_ini} vs final {n_fin}")
    # Combined initial.xyz for ML-FSM: frame1 = initial, frame2 = final
    combined = "\n".join(ini_lines) + "\n" + "\n".join(fin_lines) + "\n"
    (fsm_dir / "initial.xyz").write_text(combined)
    (fsm_dir / "chg").write_text(str(charge) + "\n")
    (fsm_dir / "mult").write_text(str(multiplicity) + "\n")
    return fsm_dir


def _parse_fsm_vfile(vfile_path: Path) -> Tuple[List[Tuple[float, float]], int, float]:
    """
    Parse FSM vfile_*.xyz; return (list of (s, E_rel), ts_node_1based, ts_energy_eV).
    """
    lines = vfile_path.read_text().strip().splitlines()
    data: List[Tuple[float, float]] = []
    i = 0
    n = len(lines)
    while i < n:
        try:
            natoms = int(lines[i])
        except ValueError:
            i += 1
            continue
        i += 1
        if i >= n:
            break
        parts = lines[i].split()
        if len(parts) >= 2:
            s = float(parts[0])
            e = float(parts[1])
            data.append((s, e))
        i += 1 + natoms
    if not data:
        return [], 1, 0.0
    ts_idx = max(range(len(data)), key=lambda j: data[j][1])
    return data, ts_idx + 1, data[ts_idx][1]


def run_fsm_calculation(args: argparse.Namespace) -> None:
    """End-to-end FSM (Freezing String Method) MEP workflow for Zn2+ adsorption."""
    try:
        size, start_distance, end_distance = validate_inputs(args)
    except ValueError as e:
        print(f"Error: {e}")
        sys.exit(1)

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    nprocs = 1 if args.no_mpi else args.nprocs
    if args.no_mpi:
        print("Running in serial mode (no MPI/parallel execution)")

    orca_gen = OrcaInputGenerator(
        method=args.method,
        basis_set=args.basis,
        solvent=args.solvent,
        solvation_model=args.solvation,
        memory=args.memory,
        nprocs=nprocs,
        scf_convergence=args.scf_convergence,
        grid=f"DefGrid{args.grid}",
    )
    calc = NebCalculator(orca_generator=orca_gen)

    # Prepare endpoint structures (writes initial.xyz, final.xyz to output_dir)
    neb_images = args.neb_images if args.neb_images is not None else calculate_neb_images(start_distance, end_distance)
    results = calc.prepare_neb_calculation(
        num_carboxyl=args.carboxyl,
        num_hydroxyl=args.hydroxyl,
        surface_size=size,
        start_distance=start_distance,
        end_distance=end_distance,
        neb_images=neb_images,
        output_dir=output_dir,
        constrain_endpoints=not getattr(args, "no_constrain_endpoints", False),
    )

    chg = getattr(args, "chg", None)
    if chg is None:
        chg = 2 - getattr(args, "carboxyl", 0)
    mult = getattr(args, "mult", 1)
    fsm_reaction_dir = _prepare_fsm_reaction_dir(
        output_dir,
        output_dir / "initial.xyz",
        output_dir / "final.xyz",
        charge=chg,
        multiplicity=mult,
    )
    print(f"\n--- FSM reaction directory prepared: {fsm_reaction_dir} ---")
    print(f"  initial.xyz (two frames: reactant, product)")
    print(f"  chg={chg}, mult={mult}")

    orca_path = getattr(args, "orca_command", None)
    if not orca_path:
        _found, orca_path = check_orca_available()
        if not _found:
            orca_path = None
    orca_available = orca_path is not None
    if not orca_available:
        print("\nORCA not found. FSM reaction files are ready; run FSM manually, e.g.:")
        print(f"  python ML-FSM/examples/fsm_example.py {fsm_reaction_dir} --calculator orca --orca_command /path/to/orca ...")
        _write_fsm_json(args, output_dir, fsm_reaction_dir, None, None, None)
        return

    orca_simple = getattr(args, "orca_simple", None) or f"{args.method} {args.basis} D3BJ CPCM({args.solvent})"
    repo_root = Path(__file__).resolve().parent.parent.parent
    fsm_script = repo_root / "ML-FSM" / "examples" / "fsm_example.py"
    if not fsm_script.exists():
        print(f"\nFSM script not found: {fsm_script}")
        _write_fsm_json(args, output_dir, fsm_reaction_dir, None, None, None)
        return

    cmd = [
        sys.executable,
        str(fsm_script),
        str(fsm_reaction_dir),
        "--calculator", "orca",
        "--orca_command", orca_path,
        "--orca_simple", orca_simple,
        "--nt", str(nprocs),
        "--chg", str(chg),
        "--mult", str(mult),
        "--nnodes_min", str(args.nnodes_min),
        "--maxiter", str(args.fsm_maxiter),
        "--maxls", str(args.fsm_maxls),
    ]
    if args.fsm_suffix:
        cmd.extend(["--suffix", args.fsm_suffix])
    if getattr(args, "verbose", False):
        cmd.append("--verbose")

    if not args.run_orca:
        print("\nFSM reaction directory is ready. Run with --run to execute FSM, e.g.:")
        print("  " + " ".join(cmd))
        _write_fsm_json(args, output_dir, fsm_reaction_dir, None, None, None)
        return

    print("\n--- Running FSM with ORCA ---")
    print(f"  ORCA: {orca_path}")
    print(f"  nnodes_min={args.nnodes_min}, maxiter={args.fsm_maxiter}, maxls={args.fsm_maxls}")
    try:
        result = subprocess.run(cmd, cwd=str(repo_root), check=False)
        if result.returncode != 0:
            print(f"\nFSM exited with code {result.returncode}. Check output above.")
    except Exception as e:
        print(f"\nError running FSM: {e}")
        _write_fsm_json(args, output_dir, fsm_reaction_dir, None, None, None)
        return

    # Parse FSM output: find latest vfile in fsm_reaction/fsm_interp_*
    fsm_out_dirs = sorted(fsm_reaction_dir.glob("fsm_interp_*"))
    path_data = None
    ts_node = None
    ts_energy_eV = None
    ngrad = None
    latest_out = None
    if fsm_out_dirs:
        latest_out = fsm_out_dirs[-1]
        vfiles = sorted(latest_out.glob("vfile_*.xyz"))
        if vfiles:
            path_data, ts_node, ts_energy_eV = _parse_fsm_vfile(vfiles[-1])
        ngrad_file = latest_out / "ngrad.txt"
        if ngrad_file.exists():
            ngrad = int(ngrad_file.read_text().strip())

    _write_fsm_json(args, output_dir, fsm_reaction_dir, path_data, ts_node, ts_energy_eV, ngrad)
    _print_fsm_summary(path_data, ts_node, ts_energy_eV, ngrad)
    # Auto-plot path if FSM output dir exists and plot script is available
    if latest_out is not None and path_data:
        _run_fsm_plot(repo_root, latest_out)
    return


def _write_fsm_json(
    args: argparse.Namespace,
    output_dir: Path,
    fsm_reaction_dir: Path,
    path_data: Optional[List[Tuple[float, float]]],
    ts_node: Optional[int],
    ts_energy_eV: Optional[float],
    ngrad: Optional[int] = None,
) -> None:
    """Write JSON summary for FSM run."""
    out = {
        "path_method": "fsm",
        "parameters": vars(args),
        "fsm_reaction_dir": str(fsm_reaction_dir),
        "path_energies": [{"s_Ang": s, "E_rel_eV": e} for (s, e) in (path_data or [])],
        "ts_guess_node": ts_node,
        "ts_guess_energy_eV": ts_energy_eV,
        "ngrad": ngrad,
    }
    try:
        with open(args.json_output, "w") as f:
            json.dump(out, f, indent=2)
        print(f"\nResults saved to: {args.json_output}")
    except Exception as e:
        print(f"\nWarning: Could not save JSON: {e}")


def _run_fsm_plot(repo_root: Path, fsm_out_dir: Path) -> None:
    """Run plot_fsm_path.py on FSM output dir and save path.png (no-op if script missing)."""
    plot_script = repo_root / "scripts" / "plot_fsm_path.py"
    out_png = fsm_out_dir / "path.png"
    if not plot_script.exists():
        return
    try:
        subprocess.run(
            [
                sys.executable,
                str(plot_script),
                str(fsm_out_dir),
                "-o",
                str(out_png),
                "--no-show",
            ],
            cwd=str(repo_root),
            check=False,
            capture_output=True,
        )
        if out_png.exists():
            print(f"Path plot saved: {out_png}")
    except Exception:
        pass


def _print_fsm_summary(
    path_data: Optional[List[Tuple[float, float]]],
    ts_node: Optional[int],
    ts_energy_eV: Optional[float],
    ngrad: Optional[int],
) -> None:
    """Print FSM summary to console."""
    print("\n" + "=" * 50)
    print("   FSM MEP SUMMARY")
    print("=" * 50)
    if path_data:
        print(f"Path nodes: {len(path_data)}")
        if ts_node is not None and ts_energy_eV is not None:
            print(f"TS guess: node {ts_node}  E_rel = {ts_energy_eV:.4f} eV")
        if ngrad is not None:
            print(f"Gradient calls: {ngrad}")
    else:
        print("(No path data parsed; check FSM output directory.)")
    print("=" * 50)


def _run_neb_uff(args: argparse.Namespace) -> None:
    """NEB workflow using UFF (in-process, fast)."""
    from ase.io import write
    from zn2_adsorption.neb_runner import run_neb_uff
    from zn2_adsorption.uff_optimizer import optimize_with_uff

    try:
        size, start_distance, end_distance = validate_inputs(args)
    except ValueError as e:
        print(f"Error: {e}")
        sys.exit(1)

    neb_images = (
        args.neb_images
        if args.neb_images is not None
        else calculate_neb_images(start_distance, end_distance)
    )
    if args.neb_images is None:
        print(f"Auto-calculated NEB images: {neb_images}")

    print("--- Preparing NEB Calculation (UFF) ---")
    print(f"Surface size: {size[0]}x{size[1]}")
    print(f"Functional groups: {args.carboxyl} carboxylate, {args.hydroxyl} hydroxyl")
    print(f"Zn2+ start distance: {start_distance} Angstroms")
    print(f"Zn2+ end distance: {end_distance} Angstroms")
    print(f"NEB images: {neb_images}")
    print("Calculator: UFF (Universal Force Field, fast)")

    calc = NebCalculator()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    constrain = not getattr(args, "no_constrain_endpoints", False)

    results = calc.prepare_endpoints_only(
        num_carboxyl=args.carboxyl,
        num_hydroxyl=args.hydroxyl,
        surface_size=size,
        start_distance=start_distance,
        end_distance=end_distance,
        output_dir=output_dir,
        constrain_endpoints=constrain,
        generate_orca_inputs=False,
    )

    initial_ase = results["initial_ase"]
    final_ase = results["final_ase"]
    freeze = results["freeze_atoms"] or []
    fg_bonds = results.get("fg_bond_constraints") or []

    if args.run_orca:
        print("\n--- Step 1: UFF optimization (initial endpoint) ---")
        if constrain:
            print("  (constrained: surface + Zn frozen, FG-surface bonds fixed)")
        try:
            if optimize_with_uff(
                initial_ase,
                freeze_indices=freeze,
                bond_constraints=fg_bonds,
            ):
                write(output_dir / "initial.xyz", initial_ase, format="xyz")
                print("  Updated initial.xyz")
            else:
                print("  Warning: UFF endpoint optimization failed")
        except ImportError as e:
            print(f"  Error: {e}")
            sys.exit(1)

        print("\n--- Step 2: UFF optimization (final endpoint) ---")
        try:
            if optimize_with_uff(
                final_ase,
                freeze_indices=freeze,
                bond_constraints=fg_bonds,
            ):
                write(output_dir / "final.xyz", final_ase, format="xyz")
                print("  Updated final.xyz")
            else:
                print("  Warning: UFF endpoint optimization failed")
        except ImportError as e:
            print(f"  Error: {e}")
            sys.exit(1)

        print("\n--- Step 3: Running NEB (UFF) ---")
        try:
            neb_result = run_neb_uff(
                initial_ase,
                final_ase,
                n_images=neb_images,
                freeze_indices=freeze,
                fg_bond_constraints=fg_bonds,
                interpolation="idpp",
                fmax=0.05,
                max_steps=200,
                output_dir=output_dir,
            )
        except ImportError as e:
            print(f"  Error: {e}")
            sys.exit(1)

        energies = neb_result["energies"]
        ts_idx = neb_result["ts_index"]
        converged = neb_result["converged"]

        EV_TO_KCAL = 23.0609
        e_init = energies[0]
        e_final = energies[-1]
        e_ts = energies[ts_idx]
        barrier = (e_ts - e_init) * EV_TO_KCAL if e_init is not None else None
        reaction = (e_final - e_init) * EV_TO_KCAL if e_init is not None else None

        print(f"  Converged: {converged}")
        print(f"  E(initial): {e_init:.4f} eV")
        print(f"  E(final):   {e_final:.4f} eV")
        print(f"  E(TS):      {e_ts:.4f} eV (image {ts_idx})")
        if barrier is not None:
            print(f"  Barrier:    {barrier:.2f} kcal/mol")
        if reaction is not None:
            print(f"  Reaction:   {reaction:.2f} kcal/mol")

        print("\n" + "=" * 50)
        print("   NEB (UFF) SUMMARY")
        print("=" * 50)
        print(f"  Path images written to {output_dir}/neb_*.xyz")
        print("  (UFF energies are qualitative; use ORCA for DFT accuracy)")
        print("=" * 50)

        if getattr(args, "plot", False):
            try:
                from zn2_adsorption.plotting import plot_neb_energy_profile
                plot_path = output_dir / "path.png"
                plot_neb_energy_profile(energies, ts_idx, plot_path, title="NEB (UFF) Energy Profile")
                print(f"\n  Energy profile plot saved: {plot_path}")
            except ImportError as e:
                print(f"\n  Warning: Could not plot ({e})")
            except Exception as e:
                print(f"\n  Warning: Plot failed: {e}")

        output_data = {
            "parameters": {k: v for k, v in vars(args).items() if isinstance(v, (str, int, float, bool, type(None)))},
            "calculator": "uff",
            "neb_run": True,
            "neb_converged": bool(converged),
            "neb_images": neb_images,
            "ts_index": int(ts_idx),
            "energies_eV": [float(e) for e in energies],
            "barrier_height_kcal_mol": float(barrier) if barrier is not None else None,
            "reaction_energy_kcal_mol": float(reaction) if reaction is not None else None,
        }
    else:
        print("\nFiles generated in:", args.output_dir)
        print("  - initial.xyz, final.xyz")
        print("\nTo run UFF NEB:")
        print(
            f"  run-adsorption --path-method neb --calculator uff "
            f"--start-distance {args.start_distance} --end-distance {args.end_distance} "
            f"--run -O {args.output_dir}"
        )
        output_data = {
            "parameters": dict(vars(args)),
            "calculator": "uff",
            "neb_run": False,
        }

    try:
        with open(args.json_output, "w") as f:
            json.dump(output_data, f, indent=2)
        print(f"\nResults summary saved to: {args.json_output}")
    except Exception as e:
        print(f"\nWarning: Could not save JSON output: {e}")


def _run_scan_uff(args: argparse.Namespace) -> None:
    """Distance scan workflow: optimize FGs at each Zn-surface distance (UFF)."""
    import numpy as np

    from zn2_adsorption.neb_calculator import NebCalculator
    from zn2_adsorption.scan_runner import run_distance_scan_uff

    try:
        size, start_distance, end_distance = validate_inputs(args)
    except ValueError as e:
        print(f"Error: {e}")
        sys.exit(1)

    scan_points = getattr(args, "scan_points", None)
    scan_step = getattr(args, "scan_step", None)
    if scan_step is not None:
        n = max(2, int(np.round((start_distance - end_distance) / scan_step)) + 1)
        distances = list(np.linspace(start_distance, end_distance, n))
    elif scan_points is not None:
        distances = list(np.linspace(start_distance, end_distance, scan_points))
    else:
        n = max(2, int(np.round((start_distance - end_distance) / 0.5)) + 1)
        distances = list(np.linspace(start_distance, end_distance, n))

    freeze_fg = getattr(args, "freeze_functional_groups", False)
    fg_mode = "FGs frozen" if freeze_fg else "C-FG bonds constrained, FGs flex"
    print("--- Preparing Distance Scan (UFF) ---")
    print(f"Surface size: {size[0]}x{size[1]}")
    print(f"Functional groups: {args.carboxyl} carboxylate, {args.hydroxyl} hydroxyl")
    print(f"Zn2+ distance range: {start_distance} to {end_distance} Angstroms")
    print(f"Scan points: {len(distances)}")
    print(f"Calculator: UFF (constrained: surface + Zn fixed, {fg_mode})")

    calc = NebCalculator()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    base = calc.prepare_scan_base(
        num_carboxyl=args.carboxyl,
        num_hydroxyl=args.hydroxyl,
        surface_size=size,
    )

    if not args.run_orca:
        print("\nFiles will be generated when --run is used.")
        print(
            f"  run-adsorption --path-method scan --calculator uff "
            f"--start-distance {args.start_distance} --end-distance {args.end_distance} "
            f"--scan-points {len(distances)} --run -O {args.output_dir}"
        )
        output_data = {
            "path_method": "scan",
            "parameters": dict(vars(args)),
            "calculator": "uff",
            "scan_run": False,
            "distances": distances,
            "energies": None,
        }
        try:
            with open(args.json_output, "w") as f:
                json.dump(output_data, f, indent=2)
            print(f"\nResults summary saved to: {args.json_output}")
        except Exception as e:
            print(f"\nWarning: Could not save JSON output: {e}")
        return

    print("\n--- Running distance scan ---")
    result = run_distance_scan_uff(
        surface_structure=base["surface_structure"],
        n_graphene=base["n_graphene"],
        fg_bond_constraints=base["fg_bond_constraints"],
        distances=distances,
        surface_builder=calc.surface_builder,
        neb_calculator=calc,
        freeze_functional_groups=freeze_fg,
        output_dir=output_dir,
        verbose=True,
        save_geometries=True,
    )

    distances_out = result["distances"]
    energies_out = result["energies"]

    print("\n==================================================")
    print("   DISTANCE SCAN (UFF) SUMMARY")
    print("==================================================")
    valid = [e for e in energies_out if e == e]
    if valid:
        print(f"  E(min): {min(valid):.4f} eV  E(max): {max(valid):.4f} eV")
        print(f"  Optimized geometries written to {output_dir}/scan_*.xyz")
    print("==================================================")

    if args.plot:
        plot_path = output_dir / "path.png"
        try:
            from zn2_adsorption.plotting import plot_scan_profile
            plot_scan_profile(distances_out, energies_out, plot_path, title="Distance Scan (UFF)")
            print(f"\n  Energy profile plot saved: {plot_path}")
        except ImportError as e:
            print(f"\n  Could not generate plot: {e}")

    output_data = {
        "path_method": "scan",
        "parameters": dict(vars(args)),
        "calculator": "uff",
        "scan_run": True,
        "distances": [float(d) for d in distances_out],
        "energies": [float(e) if e == e else None for e in energies_out],
    }
    try:
        with open(args.json_output, "w") as f:
            json.dump(output_data, f, indent=2)
        print(f"\nResults summary saved to: {args.json_output}")
    except Exception as e:
        print(f"\nWarning: Could not save JSON output: {e}")


def _run_scan_fg_sweep(args: argparse.Namespace, calculator: str) -> None:
    """Run scans over all (carboxyl, hydroxyl) combinations from (0,0) to max.

    Submits all combinations in parallel (Rowan) or runs them in parallel (UFF).
    """
    base_dir = Path(args.output_dir)
    max_c = getattr(args, "fg_sweep_max_carboxyl", 2)
    max_y = getattr(args, "fg_sweep_max_hydroxyl", 2)
    combos = [(c, y) for c in range(max_c + 1) for y in range(max_y + 1)]
    n_combos = len(combos)
    print(f"\n--- FG SWEEP: {n_combos} combinations (carboxyl 0–{max_c}, hydroxyl 0–{max_y}), parallel ---")

    def _run_one(c: int, y: int) -> None:
        subdir = base_dir / f"carboxyl_{c}_hydroxyl_{y}"
        subdir.mkdir(parents=True, exist_ok=True)
        sweep_args = copy.deepcopy(args)
        sweep_args.carboxyl = c
        sweep_args.hydroxyl = y
        sweep_args.output_dir = str(subdir)
        sweep_args.json_output = str(subdir / "scan_results.json")
        if calculator == "rowan":
            _run_scan_rowan(sweep_args)
        else:
            _run_scan_uff(sweep_args)

    max_workers = min(n_combos, 9)
    with ThreadPoolExecutor(max_workers=max_workers) as pool:
        futures = {pool.submit(_run_one, c, y): (c, y) for c, y in combos}
        for future in as_completed(futures):
            c, y = futures[future]
            try:
                future.result()
            except Exception as exc:
                print(f"\n  ERROR (carboxyl={c}, hydroxyl={y}): {exc}")


def _run_scan_rowan(args: argparse.Namespace) -> None:
    """Distance scan workflow via Rowan cloud (UMA Small).

    Targets the central carbon away from functional groups (single site).
    """
    import numpy as np

    from zn2_adsorption.neb_calculator import NebCalculator
    from zn2_adsorption.rowan_scan_runner import run_hexagon_scans_rowan, _resolve_method

    try:
        size, start_distance, end_distance = validate_inputs(args)
    except ValueError as e:
        print(f"Error: {e}")
        sys.exit(1)

    scan_points = getattr(args, "scan_points", None)
    scan_step = getattr(args, "scan_step", None)
    if scan_step is not None:
        n = max(2, int(np.round((start_distance - end_distance) / scan_step)) + 1)
        distances = list(np.linspace(start_distance, end_distance, n))
    elif scan_points is not None:
        distances = list(np.linspace(start_distance, end_distance, scan_points))
    else:
        n = max(2, int(np.round((start_distance - end_distance) / 0.5)) + 1)
        distances = list(np.linspace(start_distance, end_distance, n))

    charge = getattr(args, "chg", None)
    if charge is None:
        charge = 2 - getattr(args, "carboxyl", 0)
    multiplicity = getattr(args, "mult", 1)
    freeze_surface = not getattr(args, "relax_surface", False)
    freeze_functional_groups = getattr(args, "freeze_functional_groups", False)

    model_name = getattr(args, "model", "uma_s_omol")
    try:
        resolved_method = _resolve_method(model_name)
    except ValueError as e:
        print(f"Error: {e}")
        sys.exit(1)

    surface_mode = "surface frozen" if freeze_surface else "surface relaxed"
    fg_mode = "FGs frozen" if freeze_functional_groups else "C-FG bonds constrained, FGs flex"
    print(f"--- Preparing Distance Scan (Rowan / {resolved_method.name}) ---")
    print(f"Surface size: {size[0]}x{size[1]}")
    print(f"Functional groups: {args.carboxyl} carboxylate, {args.hydroxyl} hydroxyl")
    print(f"Zn2+ distance range: {start_distance} to {end_distance} Angstroms")
    print(f"Scan points per site: {len(distances)}")
    print(f"Calculator: Rowan ({resolved_method.name}, central site away from FGs, {surface_mode}, {fg_mode}, perpendicular approach)")

    try:
        calc = NebCalculator()
        output_dir = Path(args.output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        base = calc.prepare_scan_base(
            num_carboxyl=args.carboxyl,
            num_hydroxyl=args.hydroxyl,
            surface_size=size,
        )

        if not args.run_orca:
            print("\nRowan scan requires --run to execute. Use --run to submit to Rowan cloud.")
            print(
                f"  run-adsorption --path-method scan --calculator rowan "
                f"--start-distance {args.start_distance} --end-distance {args.end_distance} "
                f"--scan-points {len(distances)} --run -O {args.output_dir}"
            )
            return

        print(f"\n--- Running central-site scan on Rowan ({resolved_method.name}) ---")
        result = run_hexagon_scans_rowan(
            surface_structure=base["surface_structure"],
            n_graphene=base["n_graphene"],
            fg_bond_constraints=base["fg_bond_constraints"],
            distances=distances,
            neb_calculator=calc,
            charge=charge,
            multiplicity=multiplicity,
            output_dir=output_dir,
            verbose=True,
            freeze_surface=freeze_surface,
            freeze_functional_groups=freeze_functional_groups,
            method=resolved_method,
        )

        print("\n==================================================")
        print(f"   CENTRAL-SITE SCAN (Rowan {resolved_method.name}) SUMMARY")
        print("==================================================")
        for site in result.get("sites", []):
            c_label = f"C{site['target_c_idx'] + 1}"
            energies = site.get("energies", [])
            valid = [e for e in energies if e == e]
            if valid:
                print(f"  {c_label}: E(min)={min(valid):.4f} eV  E(max)={max(valid):.4f} eV")
            else:
                err = site.get("error", "no energies")
                print(f"  {c_label}: {err}")
            url = site.get("workflow_url", "N/A")
            print(f"         {url}")
        print("==================================================")

        output_data = {
            "path_method": "scan",
            "parameters": dict(vars(args)),
            "calculator": "rowan",
            "scan_run": True,
            "hexagon_indices": result.get("hexagon_indices"),
            "sites": result.get("sites"),
        }
        try:
            with open(args.json_output, "w") as f:
                json.dump(output_data, f, indent=2, default=str)
            print(f"\nResults summary saved to: {args.json_output}")
        except Exception as e:
            print(f"\nWarning: Could not save JSON output: {e}")

    except ImportError as e:
        print(f"Error: {e}")
        print("Install Rowan support: pip install zn2-adsorption[rowan]")
        sys.exit(1)
    except RuntimeError as e:
        if "API key" in str(e):
            print(f"Error: {e}")
            print("Set ROWAN_API_KEY or rowan.api_key before running.")
            sys.exit(1)
        raise


def _run_endpoints_only_uff(args: argparse.Namespace) -> None:
    """Endpoints-only workflow using UFF (fast, in-process)."""
    from ase.io import write
    from zn2_adsorption.uff_optimizer import optimize_with_uff

    try:
        size, start_distance, end_distance = validate_inputs(args)
    except ValueError as e:
        print(f"Error: {e}")
        sys.exit(1)

    print("--- Preparing Endpoints-Only Calculation (UFF) ---")
    print(f"Surface size: {size[0]}x{size[1]}")
    print(f"Functional groups: {args.carboxyl} carboxylate, {args.hydroxyl} hydroxyl")
    print(f"Zn2+ start distance: {start_distance} Angstroms")
    print(f"Zn2+ end distance: {end_distance} Angstroms")
    print("Calculator: UFF (Universal Force Field, fast)")

    calc = NebCalculator()
    output_dir = Path(args.output_dir)
    constrain = not getattr(args, "no_constrain_endpoints", False)

    results = calc.prepare_endpoints_only(
        num_carboxyl=args.carboxyl,
        num_hydroxyl=args.hydroxyl,
        surface_size=size,
        start_distance=start_distance,
        end_distance=end_distance,
        output_dir=output_dir,
        constrain_endpoints=constrain,
        generate_orca_inputs=False,
    )

    print(f"\nFiles generated in: {args.output_dir}")
    print(f"  - initial.xyz (initial structure)")
    print(f"  - final.xyz (final structure)")

    output_data = {
        "parameters": dict(vars(args)),
        "calculator": "uff",
        "files": {
            "initial_xyz": str(output_dir / "initial.xyz"),
            "final_xyz": str(output_dir / "final.xyz"),
        },
        "optimization_run": False,
        "energies": {"initial": None, "final": None},
    }

    if args.run_orca:
        initial_ase = results["initial_ase"]
        final_ase = results["final_ase"]
        freeze = results["freeze_atoms"]
        fg_bonds = results.get("fg_bond_constraints") or []

        print("\n--- Step 1: UFF optimization (initial endpoint) ---")
        if constrain:
            print("  (constrained: surface + Zn frozen, FG-surface bonds fixed, only FG relax)")
        try:
            if optimize_with_uff(
                initial_ase,
                freeze_indices=freeze,
                bond_constraints=fg_bonds,
            ):
                write(output_dir / "initial.xyz", initial_ase, format="xyz")
                print("  Updated initial.xyz with optimized geometry")
            else:
                print("  Warning: UFF optimization did not complete successfully")
        except ImportError as e:
            print(f"  Error: {e}")
            print("  Install Open Babel: pip install openbabel-wheel")

        print("\n--- Step 2: UFF optimization (final endpoint) ---")
        try:
            if optimize_with_uff(
                final_ase,
                freeze_indices=freeze,
                bond_constraints=fg_bonds,
            ):
                write(output_dir / "final.xyz", final_ase, format="xyz")
                print("  Updated final.xyz with optimized geometry")
            else:
                print("  Warning: UFF optimization did not complete successfully")
        except ImportError as e:
            print(f"  Error: {e}")

        output_data["optimization_run"] = True

    if output_data["optimization_run"]:
        print("\n" + "=" * 50)
        print("   ENDPOINTS-ONLY SUMMARY (UFF)")
        print("=" * 50)
        print("  Geometries optimized. UFF energies not reported (use ORCA for DFT energies).")
        print("=" * 50)

    try:
        with open(args.json_output, "w") as f:
            json.dump(output_data, f, indent=2)
        print(f"\nResults summary saved to: {args.json_output}")
    except Exception as e:
        print(f"\nWarning: Could not save JSON output: {e}")

    if not output_data["optimization_run"]:
        print("\nTo run UFF optimization:")
        print(
            f"  run-adsorption --endpoints-only --calculator uff "
            f"--start-distance {args.start_distance} --end-distance {args.end_distance} "
            f"--run -O {args.output_dir}"
        )


def _run_endpoints_only_orca(args: argparse.Namespace) -> None:
    """Endpoints-only workflow using ORCA (DFT, slower)."""
    try:
        size, start_distance, end_distance = validate_inputs(args)
    except ValueError as e:
        print(f"Error: {e}")
        sys.exit(1)

    basis = "def2-SVP" if (args.endpoints_only and args.basis == "def2-TZVP") else args.basis

    print("--- Preparing Endpoints-Only Calculation (ORCA) ---")
    print(f"Surface size: {size[0]}x{size[1]}")
    print(f"Functional groups: {args.carboxyl} carboxylate, {args.hydroxyl} hydroxyl")
    print(f"Zn2+ start distance: {start_distance} Angstroms")
    print(f"Zn2+ end distance: {end_distance} Angstroms")
    print(f"Basis set: {basis}")

    nprocs = 1 if args.no_mpi else args.nprocs
    if args.no_mpi:
        print("Running in serial mode (no MPI/parallel execution)")

    orca_gen = OrcaInputGenerator(
        method=args.method,
        basis_set=basis,
        solvent=args.solvent,
        solvation_model=args.solvation,
        memory=args.memory,
        nprocs=nprocs,
        scf_convergence=args.scf_convergence,
        grid=f"DefGrid{args.grid}",
    )
    calc = NebCalculator(orca_generator=orca_gen)
    output_dir = Path(args.output_dir)

    results = calc.prepare_endpoints_only(
        num_carboxyl=args.carboxyl,
        num_hydroxyl=args.hydroxyl,
        surface_size=size,
        start_distance=start_distance,
        end_distance=end_distance,
        output_dir=output_dir,
        constrain_endpoints=not getattr(args, "no_constrain_endpoints", False),
        generate_orca_inputs=True,
    )

    print(f"\nORCA input files generated in: {args.output_dir}")
    print(f"  - initial.inp (initial endpoint optimization)")
    print(f"  - final.inp (final endpoint optimization)")
    print(f"  - initial.xyz (initial structure)")
    print(f"  - final.xyz (final structure)")

    orca_available, orca_path = check_orca_available()
    output_data = {
        "parameters": dict(vars(args), basis=basis),
        "calculator": "orca",
        "files": {
            "initial_input": str(output_dir / "initial.inp"),
            "final_input": str(output_dir / "final.inp"),
            "initial_xyz": str(output_dir / "initial.xyz"),
            "final_xyz": str(output_dir / "final.xyz"),
        },
        "orca_available": orca_available,
        "orca_run": False,
        "energies": {"initial": None, "final": None},
    }

    if args.run_orca:
        if not orca_available:
            print("\nWarning: ORCA not found. Skipping execution.")
        else:
            print("\nORCA found at:", orca_path)
            output_data["orca_run"] = True

            initial_inp = output_dir / "initial.inp"
            initial_out = output_dir / "initial.out"
            initial_xyz = output_dir / "initial.xyz"
            final_inp = output_dir / "final.inp"
            final_out = output_dir / "final.out"
            final_xyz = output_dir / "final.xyz"

            print("\n--- Step 1: Optimizing initial endpoint ---")
            if not getattr(args, "no_constrain_endpoints", False):
                print("  (constrained: surface frozen, Zn-surface distance fixed, functional groups relaxed)")
            if run_orca_command(initial_inp, initial_out, orca_path, no_mpi=args.no_mpi):
                initial_energy = extract_energy_from_orca(initial_out)
                output_data["energies"]["initial"] = initial_energy
                if initial_energy is not None:
                    print(f"  Initial endpoint energy: {initial_energy:.8f} Hartree")
                if extract_optimized_geometry_from_orca(initial_out, initial_xyz):
                    print(f"  Updated {initial_xyz.name} with optimized geometry")

            print("\n--- Step 2: Optimizing final endpoint ---")
            if run_orca_command(final_inp, final_out, orca_path, no_mpi=args.no_mpi):
                final_energy = extract_energy_from_orca(final_out)
                output_data["energies"]["final"] = final_energy
                if final_energy is not None:
                    print(f"  Final endpoint energy: {final_energy:.8f} Hartree")
                if extract_optimized_geometry_from_orca(final_out, final_xyz):
                    print(f"  Updated {final_xyz.name} with optimized geometry")

    if output_data["orca_run"]:
        print("\n" + "=" * 50)
        print("   ENDPOINTS-ONLY SUMMARY (ORCA)")
        print("=" * 50)
        e_initial = output_data["energies"]["initial"]
        e_final = output_data["energies"]["final"]
        if e_initial is not None:
            print(f"E(initial):  {e_initial:15.8f} Hartree")
        else:
            print(f"E(initial):  {'(calculation incomplete)':>15}")
        if e_final is not None:
            print(f"E(final):    {e_final:15.8f} Hartree")
        else:
            print(f"E(final):    {'(calculation incomplete)':>15}")
        print("=" * 50)

    try:
        with open(args.json_output, "w") as f:
            json.dump(output_data, f, indent=2)
        print(f"\nResults summary saved to: {args.json_output}")
    except Exception as e:
        print(f"\nWarning: Could not save JSON output: {e}")

    if not output_data["orca_run"]:
        print("\nTo run the calculations manually:")
        print(f"  cd {args.output_dir}")
        print("  orca initial.inp > initial.out")
        print("  orca final.inp > final.out")


def run_endpoints_only_calculation(args: argparse.Namespace) -> None:
    """
    Endpoints-only workflow: generate and optionally run initial and final endpoint optimizations.
    Uses UFF by default (fast); use --calculator orca for DFT.
    """
    calculator = getattr(args, "calculator", "uff")
    if calculator == "uff":
        _run_endpoints_only_uff(args)
    else:
        _run_endpoints_only_orca(args)


def run_calculation(args: argparse.Namespace):
    """Main MEP workflow: FSM (default), NEB, or endpoints-only."""
    if args.endpoints_only:
        run_endpoints_only_calculation(args)
        return
    if args.path_method == "fsm":
        run_fsm_calculation(args)
        return

    calculator = getattr(args, "calculator", "orca")

    # Distance scan
    if args.path_method == "scan":
        if getattr(args, "fg_sweep", False):
            _run_scan_fg_sweep(args, calculator)
            return
        if calculator == "rowan":
            _run_scan_rowan(args)
            return
        if calculator != "uff":
            print("Warning: Scan with non-Rowan calculator uses UFF only. Ignoring --calculator orca.")
        _run_scan_uff(args)
        return

    # NEB with UFF (in-process, fast)
    if args.path_method == "neb" and calculator == "uff":
        _run_neb_uff(args)
        return

    # NEB with ORCA (DFT, external)
    try:
        size, start_distance, end_distance = validate_inputs(args)
    except ValueError as e:
        print(f"Error: {e}")
        sys.exit(1)

    # Auto-calculate NEB images if not specified
    if args.neb_images is None:
        neb_images = calculate_neb_images(start_distance, end_distance)
        print(f"Auto-calculated NEB images: {neb_images}")
    else:
        neb_images = args.neb_images
        
    print(f"--- Preparing Zn2+ NEB Calculation ---")
    print(f"Surface size: {size[0]}x{size[1]}")
    print(f"Functional groups: {args.carboxyl} carboxylate, {args.hydroxyl} hydroxyl")
    print(f"Zn2+ start distance: {start_distance} Angstroms")
    print(f"Zn2+ end distance: {end_distance} Angstroms")
    print(f"NEB images: {neb_images}")
    
    # Initialize calculator
    # If --no-mpi is set, force nprocs to 1 for serial execution
    nprocs = 1 if args.no_mpi else args.nprocs
    if args.no_mpi:
        print("Running in serial mode (no MPI/parallel execution)")
    
    orca_gen = OrcaInputGenerator(
        method=args.method,
        basis_set=args.basis,
        solvent=args.solvent,
        solvation_model=args.solvation,
        memory=args.memory,
        nprocs=nprocs,
        scf_convergence=args.scf_convergence,
        grid=f"DefGrid{args.grid}"
    )
    calc = NebCalculator(orca_generator=orca_gen)
    
    # Prepare files
    output_dir = Path(args.output_dir)
    results = calc.prepare_neb_calculation(
        num_carboxyl=args.carboxyl,
        num_hydroxyl=args.hydroxyl,
        surface_size=size,
        start_distance=start_distance,
        end_distance=end_distance,
        neb_images=neb_images,
        output_dir=output_dir,
        constrain_endpoints=not getattr(args, "no_constrain_endpoints", False),
    )
    
    print(f"\nORCA input files generated in: {args.output_dir}")
    print(f"  - initial.inp (initial endpoint optimization)")
    print(f"  - final.inp (final endpoint optimization)")
    print(f"  - neb.inp (NEB-TS calculation with PREOPT_ENDS)")
    print(f"  - initial.xyz (initial structure for NEB)")
    print(f"  - final.xyz (final structure for NEB)")
    
    # Check ORCA availability
    orca_available, orca_path = check_orca_available()
    
    output_data = {
        "parameters": vars(args),
        "files": {
            "initial_input": str(output_dir / "initial.inp"),
            "final_input": str(output_dir / "final.inp"),
            "neb_input": str(output_dir / "neb.inp"),
            "initial_xyz": str(output_dir / "initial.xyz"),
            "final_xyz": str(output_dir / "final.xyz")
        },
        "neb_images": neb_images,
        "orca_available": orca_available,
        "orca_run": False,
        "energies": {
            "initial": None,
            "final": None,
            "ts_energy": None,
            "barrier_height_kcal_mol": None,
            "reaction_energy_kcal_mol": None
        }
    }
    
    if args.run_orca:
        if not orca_available:
            print("\nWarning: ORCA not found. Skipping execution.")
        else:
            print("\nORCA found at:", orca_path)
            output_data["orca_run"] = True
            
            initial_inp = output_dir / "initial.inp"
            initial_out = output_dir / "initial.out"
            initial_xyz = output_dir / "initial.xyz"
            final_inp = output_dir / "final.inp"
            final_out = output_dir / "final.out"
            final_xyz = output_dir / "final.xyz"
            neb_inp = output_dir / "neb.inp"
            
            # Pre-optimize endpoints unless skipped
            if not args.skip_endpoint_opt:
                # Step 1: Pre-optimize initial endpoint
                print("\n--- Step 1: Pre-optimizing initial endpoint ---")
                if not getattr(args, "no_constrain_endpoints", False):
                    print("  (constrained: surface frozen, Zn-surface distance fixed, functional groups relaxed)")
                print("  (This can be skipped with --skip-endpoint-opt; PREOPT_ENDS will handle it)")
                if run_orca_command(initial_inp, initial_out, orca_path, no_mpi=args.no_mpi):
                    initial_energy = extract_energy_from_orca(initial_out)
                    output_data["energies"]["initial"] = initial_energy
                    if initial_energy is not None:
                        print(f"  Initial endpoint energy: {initial_energy:.8f} Hartree")
                    
                    # Extract optimized geometry
                    if extract_optimized_geometry_from_orca(initial_out, initial_xyz):
                        print(f"  Updated {initial_xyz.name} with optimized geometry")
                
                # Step 2: Pre-optimize final endpoint
                print("\n--- Step 2: Pre-optimizing final endpoint ---")
                if run_orca_command(final_inp, final_out, orca_path, no_mpi=args.no_mpi):
                    final_energy = extract_energy_from_orca(final_out)
                    output_data["energies"]["final"] = final_energy
                    if final_energy is not None:
                        print(f"  Final endpoint energy: {final_energy:.8f} Hartree")
                    
                    # Extract optimized geometry
                    if extract_optimized_geometry_from_orca(final_out, final_xyz):
                        print(f"  Updated {final_xyz.name} with optimized geometry")
                
                # Regenerate neb.inp using optimized XYZ files
                print("\n--- Regenerating NEB input with optimized endpoints ---")
                if initial_xyz.exists() and final_xyz.exists():
                    neb_input = orca_gen.generate_neb_input(
                        output_file=neb_inp,
                        neb_images=neb_images,
                        initial_xyz_file=initial_xyz,
                        final_xyz_file=final_xyz
                    )
                    print(f"  Regenerated {neb_inp.name} using optimized structures")
                else:
                    print(f"  Warning: Could not regenerate {neb_inp.name} - XYZ files missing")
            else:
                print("\n--- Skipping endpoint pre-optimization ---")
                print("  Using unoptimized structures; PREOPT_ENDS will optimize during NEB")
                # Use the original unoptimized XYZ files (already created by prepare_neb_calculation)
                if not initial_xyz.exists() or not final_xyz.exists():
                    print(f"  Warning: XYZ files not found. They should have been created during file preparation.")
            
            # Step 3: Run NEB calculation
            print("\n--- Step 3: Running NEB-TS calculation ---")
            neb_out = output_dir / "neb.out"
            if run_orca_command(neb_inp, neb_out, orca_path, no_mpi=args.no_mpi):
                neb_results = parse_neb_output(neb_out)
                output_data["energies"]["ts_energy"] = neb_results["ts_energy"]
                output_data["energies"]["barrier_height_kcal_mol"] = neb_results["barrier_height"]
                output_data["energies"]["reaction_energy_kcal_mol"] = neb_results["reaction_energy"]
                
                # Use endpoint energies from NEB if available, otherwise use pre-opt values
                if neb_results["initial_energy"] is not None:
                    output_data["energies"]["initial"] = neb_results["initial_energy"]
                if neb_results["final_energy"] is not None:
                    output_data["energies"]["final"] = neb_results["final_energy"]
    
    # Print summary if ORCA was run
    if output_data["orca_run"]:
        print("\n" + "="*50)
        print("   NEB CALCULATION SUMMARY")
        print("="*50)
        
        e_initial = output_data["energies"]["initial"]
        e_final = output_data["energies"]["final"]
        e_ts = output_data["energies"]["ts_energy"]
        barrier = output_data["energies"]["barrier_height_kcal_mol"]
        reaction = output_data["energies"]["reaction_energy_kcal_mol"]
        
        if e_initial is not None:
            print(f"E(initial):      {e_initial:15.8f} Hartree")
        else:
            print(f"E(initial):      {'(calculation incomplete)':>15}")
        
        if e_final is not None:
            print(f"E(final):        {e_final:15.8f} Hartree")
        else:
            print(f"E(final):         {'(calculation incomplete)':>15}")
        
        if e_ts is not None:
            print(f"E(TS):           {e_ts:15.8f} Hartree")
        else:
            print(f"E(TS):            {'(calculation incomplete)':>15}")
        
        print("-" * 50)
        
        if barrier is not None:
            print(f"Barrier Height:   {barrier:15.4f} kcal/mol")
        else:
            print(f"Barrier Height:   {'(cannot calculate)':>15}")
        
        if reaction is not None:
            print(f"Reaction Energy:  {reaction:15.4f} kcal/mol")
        else:
            print(f"Reaction Energy:  {'(cannot calculate)':>15}")
        
        print("="*50)
            
    # Save JSON output
    try:
        with open(args.json_output, 'w') as f:
            json.dump(output_data, f, indent=2)
        print(f"\nResults summary saved to: {args.json_output}")
    except Exception as e:
        print(f"\nWarning: Could not save JSON output: {e}")
    
    if not output_data["orca_run"]:
        print("\nTo run the calculations manually:")
        print(f"  cd {args.output_dir}")
        print("  orca initial.inp > initial.out")
        print("  orca final.inp > final.out")
        print("  orca neb.inp > neb.out")

def main():
    """Main entry point for the CLI."""
    args = parse_arguments(sys.argv[1:])
    run_calculation(args)

if __name__ == "__main__":
    main()
