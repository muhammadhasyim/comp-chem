# Git Push Guide for NEB Implementation

## Complete Workflow: Commit and Push to GitHub

### Step 1: Check Git Status

```bash
cd /home/mh7373/GitRepos/comp-chem
git status
```

### Step 2: Initialize Git (if not already initialized)

```bash
# Only if git is not initialized
git init
```

### Step 3: Add Remote Repository

```bash
# Check if remote already exists
git remote -v

# If no remote exists, add it:
git remote add origin https://github.com/muhammadhasyim/comp-chem.git

# If remote exists but points to wrong URL, update it:
# git remote set-url origin https://github.com/muhammadhasyim/comp-chem.git
```

### Step 4: Stage Files for Commit

```bash
# Add all the NEB implementation files
git add .gitignore
git add src/zn2_adsorption/neb_calculator.py
git add src/zn2_adsorption/cli.py
git add src/zn2_adsorption/orca_generator.py
git add tests/test_neb_calculator.py
git add tests/test_cli.py
git add README.md

# Also add other project files if they exist and aren't ignored
git add pyproject.toml
git add setup.py
git add requirements.txt
git add src/zn2_adsorption/__init__.py
git add src/zn2_adsorption/surface_builder.py
git add src/zn2_adsorption/adsorption_calculator.py
git add tests/*.py
git add scripts/*.py

# Verify what will be committed
git status
```

### Step 5: Commit Changes

```bash
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

### Step 6: Push to GitHub

```bash
# Push to main branch (or master if that's your default)
git branch -M main
git push -u origin main

# If you get an error about the branch, try:
# git push -u origin main --force  # Only if repository is empty on GitHub
```

### Alternative: If Repository Already Has Content

If the GitHub repository already has commits, you may need to pull first:

```bash
# Pull and merge (if needed)
git pull origin main --allow-unrelated-histories

# Then push
git push -u origin main
```

## Troubleshooting

### If you get authentication errors:
- Use a Personal Access Token instead of password
- Or set up SSH keys: https://docs.github.com/en/authentication/connecting-to-github-with-ssh

### If you get "refusing to merge unrelated histories":
```bash
git pull origin main --allow-unrelated-histories
# Resolve any conflicts, then:
git push -u origin main
```

### Verify Push Success:
Visit: https://github.com/muhammadhasyim/comp-chem
