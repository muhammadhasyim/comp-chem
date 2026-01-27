#!/bin/bash
# Complete git setup and push script

cd /home/mh7373/GitRepos/comp-chem

echo "=== Step 1: Initialize Git ==="
git init

echo ""
echo "=== Step 2: Add Remote ==="
git remote add origin https://github.com/muhammadhasyim/comp-chem.git

echo ""
echo "=== Step 3: Stage Files ==="
git add .gitignore
git add src/zn2_adsorption/neb_calculator.py
git add src/zn2_adsorption/cli.py
git add src/zn2_adsorption/orca_generator.py
git add tests/test_neb_calculator.py
git add tests/test_cli.py
git add README.md
git add pyproject.toml setup.py requirements.txt
git add src/zn2_adsorption/*.py
git add tests/*.py
git add scripts/*.py

echo ""
echo "=== Step 4: Check Status ==="
git status

echo ""
echo "=== Step 5: Commit ==="
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

echo ""
echo "=== Step 6: Set Branch and Push ==="
git branch -M main
git push -u origin main

echo ""
echo "=== Done! ==="
echo "Repository pushed to: https://github.com/muhammadhasyim/comp-chem"
