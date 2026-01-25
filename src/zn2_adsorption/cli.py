import argparse
import sys
import os
import json
import shutil
import subprocess
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
        description="Zn2+ Nudged Elastic Band (NEB) Path Calculator CLI"
    )
    
    # Surface parameters
    parser.add_argument(
        "--carboxyl", "-c",
        type=int,
        default=0,
        help="Number of carboxyl (-COOH) groups (default: 0)"
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
    
    # Execution and output parameters
    parser.add_argument(
        "--output-dir", "-O",
        type=str,
        default="./orca_inputs",
        help="Output directory for ORCA files (default: './orca_inputs')"
    )
    parser.add_argument(
        "--run-orca",
        action="store_true",
        default=False,
        help="Attempt to run ORCA if available (default: False)"
    )
    parser.add_argument(
        "--json-output",
        type=str,
        default="./neb_results.json",
        help="JSON output file path (default: './neb_results.json')"
    )
    
    return parser.parse_args(args)

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

def run_orca_command(input_file: Path, output_file: Path, orca_path: str) -> bool:
    """
    Execute ORCA command.
    Returns True if ORCA started successfully, False if there was an input error.
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
        
        with open(output_file, 'w') as out_f:
            result = subprocess.run(
                [str(orca_exe), str(input_path)],
                stdout=out_f,
                stderr=subprocess.STDOUT,
                check=False,  # Don't raise on non-zero exit - check output file instead
                env=env  # Use environment with library path
            )
        
        # Check output file for actual errors (not just exit code)
        # ORCA may return non-zero if interrupted, but still produce valid output
        if output_file.exists() and output_file.stat().st_size > 0:
            try:
                # Try reading as text with error handling for binary data
                with open(output_file, 'r', encoding='utf-8', errors='replace') as f:
                    content = f.read()
            except Exception as e:
                # Fallback: try reading as binary and decode with error handling
                try:
                    with open(output_file, 'rb') as f:
                        content = f.read().decode('utf-8', errors='replace')
                except Exception:
                    # If we can't read it at all, assume it's still running or corrupted
                    print(f"  Warning: Could not read {output_file.name} (may contain binary data). Assuming calculation is running.")
                    return True
            
            # Check for actual input errors
            if "INPUT ERROR" in content or "UNRECOGNIZED" in content:
                print(f"  Error: ORCA input error for {input_file.name}. Check {output_file.name} for details.")
                return False
            # Check for segmentation fault
            if "segmentation fault" in content.lower() or "SIGSEGV" in content:
                print(f"  Error: ORCA crashed (segmentation fault) for {input_file.name}. Check {output_file.name} for details.")
                return False
            # If ORCA started processing, consider it successful
            if "Starting time:" in content and ("Geometry Optimization" in content or "SCF" in content or "ORCA TERMINATED" in content or "NEB" in content):
                if result.returncode == 0:
                    print(f"  ✓ ORCA completed successfully for {input_file.name}")
                elif result.returncode == -11:
                    # SIGSEGV - segmentation fault
                    print(f"  ✗ ORCA crashed (segmentation fault) for {input_file.name}. Check {output_file.name} for details.")
                    return False
                else:
                    print(f"  ⚠ ORCA was interrupted for {input_file.name} (exit code {result.returncode}). Check {output_file.name} for progress.")
                return True
        
        # If we get here, something went wrong
        if result.returncode != 0:
            if result.returncode == -11:
                print(f"  ✗ ORCA crashed (segmentation fault) for {input_file.name}. Check {output_file.name} for details.")
            else:
                print(f"  Error: ORCA failed for {input_file.name} (exit code {result.returncode}). Check {output_file.name} for details.")
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

def run_calculation(args: argparse.Namespace):
    """Main NEB workflow orchestration."""
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
    print(f"Functional groups: {args.carboxyl} carboxyl, {args.hydroxyl} hydroxyl")
    print(f"Zn2+ start distance: {start_distance} Angstroms")
    print(f"Zn2+ end distance: {end_distance} Angstroms")
    print(f"NEB images: {neb_images}")
    
    # Initialize calculator
    orca_gen = OrcaInputGenerator(
        method=args.method,
        basis_set=args.basis,
        solvent=args.solvent,
        solvation_model=args.solvation,
        memory=args.memory,
        nprocs=args.nprocs,
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
        output_dir=output_dir
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
            
            # Step 1: Pre-optimize initial endpoint (unconstrained)
            print("\n--- Step 1: Pre-optimizing initial endpoint ---")
            initial_inp = output_dir / "initial.inp"
            initial_out = output_dir / "initial.out"
            initial_xyz = output_dir / "initial.xyz"
            if run_orca_command(initial_inp, initial_out, orca_path):
                initial_energy = extract_energy_from_orca(initial_out)
                output_data["energies"]["initial"] = initial_energy
                if initial_energy is not None:
                    print(f"  Initial endpoint energy: {initial_energy:.8f} Hartree")
                
                # Extract optimized geometry
                if extract_optimized_geometry_from_orca(initial_out, initial_xyz):
                    print(f"  Updated {initial_xyz.name} with optimized geometry")
            
            # Step 2: Pre-optimize final endpoint (unconstrained)
            print("\n--- Step 2: Pre-optimizing final endpoint ---")
            final_inp = output_dir / "final.inp"
            final_out = output_dir / "final.out"
            final_xyz = output_dir / "final.xyz"
            if run_orca_command(final_inp, final_out, orca_path):
                final_energy = extract_energy_from_orca(final_out)
                output_data["energies"]["final"] = final_energy
                if final_energy is not None:
                    print(f"  Final endpoint energy: {final_energy:.8f} Hartree")
                
                # Extract optimized geometry
                if extract_optimized_geometry_from_orca(final_out, final_xyz):
                    print(f"  Updated {final_xyz.name} with optimized geometry")
            
            # Regenerate neb.inp using optimized XYZ files
            print("\n--- Regenerating NEB input with optimized endpoints ---")
            neb_inp = output_dir / "neb.inp"
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
            
            # Step 3: Run NEB calculation
            print("\n--- Step 3: Running NEB-TS calculation ---")
            neb_out = output_dir / "neb.out"
            if run_orca_command(neb_inp, neb_out, orca_path):
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
