# Zn2+ Nudged Elastic Band (NEB) Path Calculator CLI

A plug-and-play command-line script for students to calculate minimum energy paths for Zn2+ adsorption on functionalized carbon surfaces using ORCA's NEB-TS method.

## Overview

This tool automates the workflow for finding minimum energy paths (MEP) and transition states for Zn²⁺ adsorption on functionalized graphene surfaces using Nudged Elastic Band (NEB) calculations. It uses a hybrid approach:
- **pymatgen**: Provides periodic boundary conditions to eliminate edge effects.
- **GOPY Logic**: Implements reliable functionalization patterns for carboxyl (-COOH) and hydroxyl (-OH) groups.
- **ORCA**: Generates NEB-TS input files with constrained endpoint optimizations.

## Installation

### Prerequisites

- Python 3.9 or higher
- ORCA 6.1.1 or later (optional, for running calculations)
- pip package manager

### Step 1: Clone the Repository

```bash
git clone https://github.com/muhammadhasyim/comp-chem.git
cd comp-chem
```

### Step 2: Install Python Dependencies

Install the package and its dependencies in editable mode:

```bash
pip install -e .
```

This will install:
- `pymatgen` - Materials structure generation
- `ase` - Atomic Simulation Environment
- `numpy`, `scipy` - Numerical computing
- `matplotlib` - Plotting (optional)
- `rdkit` - Molecular chemistry toolkit

### Step 3: Verify Installation

Check that the CLI tool is available:

```bash
run-adsorption --help
```

You should see the help message with all available command-line options.

### Step 4: Set Up ORCA (Optional)

If you want to run calculations automatically (not just generate input files):

1. **Download ORCA**: Get ORCA 6.1.1 or later from [ORCA Forum](https://orcaforum.kofo.mpg.de)

2. **Extract ORCA**: Place the ORCA folder in your project directory or system PATH

3. **Verify ORCA Detection**: The tool will automatically detect ORCA in:
   - Project directory (folders containing "orca" in the name)
   - System PATH (excluding system screen reader)

   To test:
   ```bash
   run-adsorption --start-distance 5.0 --end-distance 2.5 --run-orca
   ```

### Troubleshooting Installation

**Issue: `ERROR: Could not find a version that satisfies the requirement install`**
- **Solution**: Use exactly `pip install -e .` or `pip install .`. Do not combine commands (e.g., avoid `pip install . install -e .`).

**Issue: `ModuleNotFoundError: No module named 'pymatgen'`**
- **Solution**: Install dependencies explicitly:
  ```bash
  pip install -r requirements.txt
  pip install -e .
  ```

**Issue: `ERROR: Could not find a version that satisfies the requirement rdkit-pypi`**
- **Solution**: The package name is `rdkit`, not `rdkit-pypi`. Update `requirements.txt` if needed.

**Issue: ORCA not detected**
- **Solution**: 
  - Place ORCA folder in the project directory (e.g., `comp-chem/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg/`)
  - Or add ORCA to your system PATH
  - The tool will automatically set `LD_LIBRARY_PATH` for ORCA's shared libraries

## Usage

The tool provides a command-line interface `run-adsorption` for NEB calculations.

### Basic Usage: Generate Input Files Only

Generate ORCA input files without running calculations:

```bash
run-adsorption \
    --carboxyl 2 \
    --hydroxyl 1 \
    --size 4x4 \
    --start-distance 5.0 \
    --end-distance 2.5
```

This creates:
- `initial.inp` - Initial endpoint optimization
- `final.inp` - Final endpoint optimization  
- `neb.inp` - NEB-TS calculation
- `initial.xyz` / `final.xyz` - Structure files

### Run with ORCA (Automatic Execution)

To automatically run the calculations (requires ORCA to be installed):

```bash
run-adsorption \
    --carboxyl 2 \
    --hydroxyl 1 \
    --size 4x4 \
    --start-distance 5.0 \
    --end-distance 2.5 \
    --run-orca
```

The workflow will:
1. Generate input files
2. Run initial endpoint optimization
3. Run final endpoint optimization
4. Extract optimized geometries
5. Regenerate NEB input with optimized structures
6. Run NEB-TS calculation
7. Generate summary JSON file

### Laptop-Optimized Settings

For faster calculations on limited resources:

```bash
run-adsorption \
    --carboxyl 0 \
    --hydroxyl 0 \
    --size 2x2 \
    --start-distance 5.0 \
    --end-distance 2.5 \
    --basis def2-SVP \
    --grid 3 \
    --scf-convergence NormalSCF \
    --memory 8GB \
    --nprocs 2 \
    --run-orca
```

**Settings explanation:**
- `--basis def2-SVP`: Smaller basis set (faster than def2-TZVP)
- `--grid 3`: Lower integration grid (faster, slightly less accurate)
- `--scf-convergence NormalSCF`: Less strict convergence (faster)
- `--memory 8GB`: Lower memory allocation
- `--nprocs 2`: Fewer processors

### Advanced Options

Full control over all parameters:

```bash
run-adsorption \
    --carboxyl 3 \
    --hydroxyl 2 \
    --size 6x6 \
    --start-distance 6.0 \
    --end-distance 2.0 \
    --neb-images 15 \
    --method PBE0 \
    --basis def2-TZVP \
    --solvent water \
    --solvation CPCM \
    --scf-convergence TightSCF \
    --grid 4 \
    --memory 32GB \
    --nprocs 8 \
    --output-dir ./neb_calc \
    --run-orca \
    --json-output ./neb_results.json
```

### Running Calculations Manually

If you prefer to run ORCA manually:

1. **Generate input files** (without `--run-orca`):
   ```bash
   run-adsorption --start-distance 5.0 --end-distance 2.5
   ```

2. **Run initial endpoint optimization**:
   ```bash
   cd orca_inputs
   export LD_LIBRARY_PATH=/path/to/orca/lib:$LD_LIBRARY_PATH
   /path/to/orca initial.inp > initial.out 2>&1
   ```

3. **Run final endpoint optimization**:
   ```bash
   /path/to/orca final.inp > final.out 2>&1
   ```

4. **Regenerate NEB input** (the tool will do this automatically if you use `--run-orca`):
   ```bash
   cd ..
   run-adsorption --start-distance 5.0 --end-distance 2.5
   # This will regenerate neb.inp with optimized structures
   ```

5. **Run NEB calculation**:
   ```bash
   cd orca_inputs
   /path/to/orca neb.inp > neb.out 2>&1
   ```

### Command-Line Flags

| Flag | Short | Description | Default |
|------|-------|-------------|---------|
| `--carboxyl` | `-c` | Number of carboxyl (-COOH) groups | 0 |
| `--hydroxyl` | `-y` | Number of hydroxyl (-OH) groups | 0 |
| `--size` | `-s` | Surface size (e.g., '4x4' or '4') | '4x4' |
| `--start-distance` | | Initial Zn²⁺ distance from surface (Å) | **Required** |
| `--end-distance` | | Final Zn²⁺ distance from surface (Å) | **Required** |
| `--neb-images` | | Number of NEB images (auto-calculated if not specified) | Auto |
| `--method` | | DFT functional (e.g., B3LYP, PBE0) | 'B3LYP' |
| `--basis` | | Basis set (e.g., def2-TZVP) | 'def2-TZVP' |
| `--solvent` | | Implicit solvent name | 'water' |
| `--solvation` | | Solvation model ('CPCM' or 'SMD') | 'CPCM' |
| `--memory` | | Memory allocation (e.g., '16GB', '32GB') | '16GB' |
| `--nprocs` | | Number of processors | 4 |
| `--no-mpi` | | Disable MPI/parallel execution - run ORCA in serial mode (nprocs=1) | False |
| `--scf-convergence` | | SCF convergence criteria ('TightSCF', 'NormalSCF', 'LooseSCF') | 'TightSCF' |
| `--grid` | | Integration grid quality (1-7, higher is better but slower) | 4 |
| `--output-dir` | `-O` | Directory for ORCA files | './orca_inputs' |
| `--run-orca` | | Attempt to run ORCA if available | False |
| `--skip-endpoint-opt` | | Skip pre-optimization of endpoints (PREOPT_ENDS will handle it). Faster but may converge slower. | False |
| `--json-output` | | Path for results JSON file | './neb_results.json' |

## Workflow

The NEB calculation proceeds in three main steps (unless `--skip-endpoint-opt` is used):

1. **Pre-optimize initial endpoint** (optional, skipped with `--skip-endpoint-opt`): Unconstrained geometry optimization starting from Zn²⁺ at `--start-distance`
   - Generates `initial.inp` and runs optimization
   - Extracts optimized geometry to `initial.xyz`
   - Energy saved to results

2. **Pre-optimize final endpoint** (optional, skipped with `--skip-endpoint-opt`): Unconstrained geometry optimization starting from Zn²⁺ at `--end-distance`
   - Generates `final.inp` and runs optimization
   - Extracts optimized geometry to `final.xyz`
   - Energy saved to results

3. **NEB-TS calculation**: Find minimum energy path and transition state between endpoints
   - Uses optimized `initial.xyz` and `final.xyz` files (if pre-optimized) or unoptimized structures (if `--skip-endpoint-opt`)
   - Uses `* XYZFILE` to read structures
   - Uses `PREOPT_ENDS TRUE` to optimize endpoints during NEB (if not pre-optimized)
   - Finds transition state and calculates barrier height

**Note**: 
- By default, endpoints are pre-optimized separately for better starting geometries and separate energy values.
- With `--skip-endpoint-opt`, ORCA's `PREOPT_ENDS TRUE` keyword will optimize endpoints during the NEB calculation itself, saving computation time but potentially requiring more NEB iterations.
- The endpoint optimizations are unconstrained (no distance constraints).

## Output Files

The tool generates the following files in the output directory:

- `initial.inp` / `initial.out` - Initial endpoint optimization
- `final.inp` / `final.out` - Final endpoint optimization
- `neb.inp` / `neb.out` - NEB-TS calculation
- `initial.xyz` - Initial structure (for NEB)
- `final.xyz` - Final structure (for NEB)
- `neb_results.json` - Summary of results (energies, barrier height, reaction energy)

## Auto-calculation of NEB Images

If `--neb-images` is not specified, the number of images is automatically calculated based on the distance difference:

- Formula: `images = max(5, min(20, int(distance_diff / 0.5) + 1))`
- Minimum: 5 images
- Maximum: 20 images
- Approximately 1 image per 0.5 Å of distance change

## Project Structure

- `src/zn2_adsorption/`: Core logic and modules.
  - `surface_builder.py`: Hybrid pymatgen + GOPY logic for surface construction.
  - `orca_generator.py`: ORCA input generation with constraints and NEB-TS support.
  - `neb_calculator.py`: NEB workflow orchestration.
  - `cli.py`: Command-line interface and validation.
- `scripts/`: Executable scripts for users.
- `tests/`: Unit and integration tests.

## References

- GOPY: [muhammadhasyim/GOPY](https://github.com/muhammadhasyim/GOPY)
- Pymatgen: [pymatgen.org](https://pymatgen.org)
- ORCA: [orcaforum.kofo.mpg.de](https://orcaforum.kofo.mpg.de)
- NEB Method: See ORCA manual for NEB-TS implementation details
