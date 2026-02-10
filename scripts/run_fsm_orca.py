#!/usr/bin/env python3
"""
Run ML-FSM with the ORCA calculator.

Finds the ORCA chemistry executable (not the GNOME orca) via ORCA_COMMAND
or PATH (orca_6, orca_5, etc.), then runs fsm_example.py with --calculator orca.

Usage:
  export ORCA_COMMAND=/path/to/orca_5   # optional if ORCA is in PATH
  python scripts/run_fsm_orca.py orca_inputs/fsm_reaction [fsm_example.py options]

Example:
  python scripts/run_fsm_orca.py orca_inputs/fsm_reaction --nnodes_min 6 --nt 1
"""

import os
import shutil
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
ML_FSM_EXAMPLES = REPO_ROOT / "ML-FSM" / "examples"
FSM_SCRIPT = ML_FSM_EXAMPLES / "fsm_example.py"

ORCA_NAMES = ["orca_6", "orca_5", "orca_5_0_5", "orca_5_0_4", "orca_5_0_3"]


def find_orca() -> str | None:
    if os.environ.get("ORCA_COMMAND"):
        path = Path(os.environ["ORCA_COMMAND"]).resolve()
        if path.exists():
            return str(path)
    for name in ORCA_NAMES:
        path = shutil.which(name)
        if path:
            return path
    return None


def main() -> int:
    if not FSM_SCRIPT.exists():
        print(f"FSM script not found: {FSM_SCRIPT}", file=sys.stderr)
        return 1
    args = sys.argv[1:]
    orca_path = None
    i = 0
    while i < len(args):
        if args[i] == "--orca_command" and i + 1 < len(args):
            orca_path = args[i + 1]
            args = args[:i] + args[i + 2:]
            break
        i += 1
    if not orca_path:
        orca_path = find_orca()
    if not orca_path:
        print(
            "ORCA chemistry executable not found. Set ORCA_COMMAND, pass --orca_command /path/to/orca, or add orca_6/orca_5 to PATH.",
            file=sys.stderr,
        )
        return 1
    orca_path = str(Path(orca_path).resolve())
    argv = [
        sys.executable,
        str(FSM_SCRIPT),
        *args,
        "--calculator",
        "orca",
        "--orca_command",
        orca_path,
    ]
    return subprocess.run(argv, cwd=str(REPO_ROOT)).returncode


if __name__ == "__main__":
    sys.exit(main())
