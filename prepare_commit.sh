#!/bin/bash
# Git commit script for NEB implementation

echo "=== Checking git status ==="
git status --short

echo ""
echo "=== Files to be committed ==="
git status --short | grep -E "^[AM]" | awk '{print $2}'

echo ""
echo "=== Adding files ==="
git add .gitignore
git add src/zn2_adsorption/neb_calculator.py
git add src/zn2_adsorption/cli.py
git add src/zn2_adsorption/orca_generator.py
git add tests/test_neb_calculator.py
git add tests/test_cli.py
git add README.md

echo ""
echo "=== Staged files ==="
git status --short

echo ""
echo "=== Ready to commit ==="
echo "Run: git commit -m 'Implement NEB workflow for Zn2+ adsorption path calculations'"
