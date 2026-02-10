#!/usr/bin/env python3
"""
End-to-end MEP pipeline for Zn²⁺ adsorbing to a functionalized carbon surface.

By default uses FSM (Freezing String Method) with ORCA. Use --path-method neb for NEB.

Usage (from repo root, after: pip install -e .  or  PYTHONPATH=src):
  python scripts/run_adsorption.py --start-distance 5 --end-distance 2.5 --run-orca
  python scripts/run_adsorption.py --path-method neb --start-distance 5 --end-distance 2.5 --run-orca

Or use the installed entry point:
  run-adsorption --start-distance 5 --end-distance 2.5 --run-orca
"""
import sys
from pathlib import Path

# Allow running from repo root when package is not installed (src layout)
_script_dir = Path(__file__).resolve().parent
_repo_root = _script_dir.parent
_src = _repo_root / "src"
if _src.exists() and str(_src) not in sys.path:
    sys.path.insert(0, str(_src))

from zn2_adsorption.cli import main

if __name__ == "__main__":
    main()
