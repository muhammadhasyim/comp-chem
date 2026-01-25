# Zn2+ Nudged Elastic Band (NEB) Path Calculator CLI

A plug-and-play command-line script for students to calculate minimum energy paths for Zn2+ adsorption on functionalized carbon surfaces using ORCA's NEB-TS method.

## Overview

This tool automates the workflow for finding minimum energy paths (MEP) and transition states for Zn²⁺ adsorption on functionalized graphene surfaces using Nudged Elastic Band (NEB) calculations. It uses a hybrid approach:
- **pymatgen**: Provides periodic boundary conditions to eliminate edge effects.
- **GOPY Logic**: Implements reliable functionalization patterns for carboxyl (-COOH) and hydroxyl (-OH) groups.
- **ORCA**: Generates NEB-TS input files with constrained endpoint optimizations.

## Installation

1. Clone the repository:
   ```bash
   git clone https://github.com/muhammadhasyim/comp-chem.git
   cd comp-chem
   ```

2. Install the package in editable mode:
   ```bash
   pip install -e .
   ```

3. Verify installation:
   ```bash
   run-adsorption --help
   ```

### Troubleshooting

If you encounter an error like `ERROR: Could not find a version that satisfies the requirement install`, make sure you are using exactly `pip install -e .` or `pip install .`. Do not combine them (e.g., avoid `pip install . install -e .`).

## Usage

The tool provides a command-line interface `run-adsorption` for NEB calculations.

### Basic Example

Calculate NEB path from Zn²⁺ at 5.0 Å to 2.5 Å from the surface:

```bash
run-adsorption --carboxyl 2 --hydroxyl 1 --size 4x4 \
    --start-distance 5.0 --end-distance 2.5
```

### Run with ORCA

To automatically run the calculations (if ORCA is available):

```bash
run-adsorption --carboxyl 2 --hydroxyl 1 --size 4x4 \
    --start-distance 5.0 --end-distance 2.5 --run-orca
```

### Advanced Options

```bash
run-adsorption \
    --carboxyl 3 \
    --hydroxyl 2 \
    --size 6x6 \
    --start-distance 6.0 \
    --end-distance 2.0 \
    --neb-images 15 \
    --method PBE0 \
    --basis def2-SVP \
    --solvent water \
    --solvation CPCM \
    --memory 32GB \
    --nprocs 8 \
    --output-dir ./neb_calc \
    --run-orca \
    --json-output ./neb_results.json
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
| `--output-dir` | `-O` | Directory for ORCA files | './orca_inputs' |
| `--run-orca` | | Attempt to run ORCA if available | False |
| `--json-output` | | Path for results JSON file | './neb_results.json' |

## Workflow

The NEB calculation proceeds in three steps:

1. **Pre-optimize initial endpoint**: Constrained geometry optimization with Zn²⁺ at `--start-distance`
2. **Pre-optimize final endpoint**: Constrained geometry optimization with Zn²⁺ at `--end-distance`
3. **NEB-TS calculation**: Find minimum energy path and transition state between endpoints

Both endpoints are optimized with distance constraints to ensure they maintain their specified distances before the NEB calculation.

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
