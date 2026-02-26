# Scripts

Executable scripts for the Zn²⁺ MEP pipeline. Run from repo root with `python scripts/<script>.py` or `pip install -e .` for the main CLI.

## MEP / Adsorption Pipeline

| Script | Description |
|--------|-------------|
| **run_adsorption.py** | Main entry point. End-to-end MEP pipeline (FSM by default, NEB with `--path-method neb`). Prepares structures, optionally runs ORCA. |
| **run_fsm_orca.py** | Run ML-FSM with ORCA calculator. Use when you have a prepared `fsm_reaction/` directory. |

## Structure Preparation (No Optimization)

| Script | Description |
|--------|-------------|
| **prepare_endpoints.py** | Build `initial.xyz` and `final.xyz` at start/end Zn²⁺ distances. No optimization. |
| **prepare_initial_final_xyz.py** | Same as prepare_endpoints; alternate entry point. |
| **prepare_path_xyz.py** | Build a series of XYZ files (default 15) along the path with Zn²⁺ at evenly spaced distances. No optimization. |

## Analysis & Utilities

| Script | Description |
|--------|-------------|
| **neb_interp_to_csv.py** | Extract the latest iteration from ORCA NEB `.interp` file to CSV. |
| **summarize_fsm_output.py** | Summarize FSM run output. |
| **plot_fsm_path.py** | Plot FSM path results. |
