# MEP Pipeline: Zn²⁺ Adsorption on Functionalized Carbon

End-to-end workflow for computing the minimum energy path (MEP) of Zn²⁺ adsorbing to a functionalized carbon surface.

## Default: FSM (Freezing String Method)

The pipeline uses **FSM with ORCA** by default:

1. **Build endpoints** – Surface (size, carboxyl/hydroxyl), initial and final Zn²⁺ distances.
2. **Prepare fsm_reaction** – `<output-dir>/fsm_reaction/` with combined `initial.xyz` (two frames), `chg`, `mult`.
3. **Run FSM** – ML-FSM with ORCA finds the path and TS guess (when `--run-orca` is set).
4. **Results** – `mep_results.json`, path plot `path.png` in the FSM output dir, and console summary.

### Quick run

**Note:** There is no positional path. Use **-O** / **--output-dir** for the output directory (fsm_reaction is created under it).

```bash
# Prepare only (no ORCA run)
python scripts/run_adsorption.py --start-distance 5 --end-distance 2.5 -O my_run

# Full run (prepare + FSM with ORCA)
python scripts/run_adsorption.py --start-distance 5 --end-distance 2.5 --run-orca -O my_run

# Example: write to orca_inputs, 16 threads, custom ORCA path and "!" line
python scripts/run_adsorption.py --start-distance 5 --end-distance 0.5 --run-orca -O orca_inputs \
  --orca-command /path/to/orca --orca-simple "B3LYP def2-TZVP D3BJ CPCM(water)" \
  --nprocs 16 --chg 2 --fsm-suffix run1 --verbose
```

### Main options (FSM)

| Option | Default | Description |
|--------|--------|-------------|
| `--path-method` | fsm | Use `fsm` or `neb`. |
| `--start-distance` | (required) | Zn²⁺ initial distance (Å). |
| `--end-distance` | (required) | Zn²⁺ final distance (Å). |
| `--size` | 4x4 | Surface size (e.g. 2x2, 4x4). |
| `--carboxyl` / `--hydroxyl` | 0 | Functional group counts. |
| `--chg` | 2 | Total charge (Zn²⁺@surface). |
| `--mult` | 1 | Spin multiplicity. |
| `--nprocs` | 4 | Number of processors/threads for ORCA (passed as `--nt` to FSM). |
| `--orca-command` | (auto) | Path to ORCA executable (override auto-detect). |
| `--orca-simple` | (built) | ORCA "!" line override (e.g. `B3LYP def2-TZVP D3BJ CPCM(water)`). Default: built from `--method`, `--basis`, `--solvent`. |
| `--verbose` | False | Verbose output (passed to FSM). |
| `--nnodes-min` | 10 | FSM: minimum nodes on the string. |
| `--fsm-maxiter` | 2 | FSM: max opt iterations per node. |
| `--run-orca` | False | Run ORCA (and FSM when path-method=fsm). |
| `-O` / `--output-dir` | ./orca_inputs | Output directory. |
| `--json-output` | ./mep_results.json | JSON summary path. |

ORCA is auto-detected (e.g. `orca*` folder in repo or PATH) unless `--orca-command` is set. All ORCA options (`--method`, `--basis`, `--solvent`, `--nprocs`, etc.) apply to both FSM and NEB.

## NEB (Nudged Elastic Band)

To use the original NEB workflow:

```bash
python scripts/run_adsorption.py --path-method neb --start-distance 5 --end-distance 2.5 --run-orca -O my_run
```

This prepares `initial.inp`, `final.inp`, `neb.inp`, and optionally runs endpoint optimizations and NEB-TS.

## Output layout (FSM)

- **`<output-dir>/`** – `initial.xyz`, `final.xyz`, ORCA inputs (from shared endpoint prep).
- **`<output-dir>/fsm_reaction/`** – `initial.xyz` (two-frame), `chg`, `mult`.
- **`<output-dir>/fsm_reaction/fsm_interp_.../`** – FSM results: `vfile_*.xyz`, `ngrad.txt`, **`path.png`** (auto-generated when pipeline runs with `--run-orca`).
- **`mep_results.json`** – Path energies, TS guess node, ngrad, parameters.

See [FSM_output_guide.md](FSM_output_guide.md) for interpreting vfiles and the path plot.

## Post-processing

- **Path plot:** `python scripts/plot_fsm_path.py <fsm_output_dir> -o path.png`
- **Text summary:** `python scripts/summarize_fsm_output.py <fsm_output_dir>`
