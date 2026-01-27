# Git Commit Guide for NEB Implementation

## Files Changed/Created

### New Files:
- `src/zn2_adsorption/neb_calculator.py` - NEB workflow orchestration
- `tests/test_neb_calculator.py` - NEB calculator tests
- `.gitignore` - Exclude ORCA inputs and software

### Modified Files:
- `src/zn2_adsorption/cli.py` - Updated for NEB workflow (replaced adsorption energy)
- `src/zn2_adsorption/orca_generator.py` - Added NEB-TS input generation
- `tests/test_cli.py` - Updated for new CLI arguments
- `README.md` - Updated documentation for NEB workflow

## Git Commands

```bash
# Add all the relevant files
git add .gitignore
git add src/zn2_adsorption/neb_calculator.py
git add src/zn2_adsorption/cli.py
git add src/zn2_adsorption/orca_generator.py
git add tests/test_neb_calculator.py
git add tests/test_cli.py
git add README.md

# Verify what will be committed
git status

# Commit with descriptive message
git commit -m "Implement NEB workflow for Zn2+ adsorption path calculations

- Replace adsorption energy workflow with Nudged Elastic Band (NEB) calculations
- Add NebCalculator class for NEB workflow orchestration
- Update CLI: replace --distance with --start-distance and --end-distance
- Add NEB-TS input generation to OrcaInputGenerator
- Implement auto-calculation of NEB images based on distance difference
- Add PREOPT_ENDS support for automatic endpoint optimization
- Remove constraints from endpoint optimizations (ORCA handles automatically)
- Add comprehensive tests for NEB workflow
- Update README with NEB usage examples and documentation
- Add .gitignore to exclude ORCA inputs and software folder"
```

## Commit Message (Alternative - Single Line)

```bash
git commit -m "Implement NEB workflow: replace adsorption energy with NEB-TS path finding for Zn2+ adsorption"
```
