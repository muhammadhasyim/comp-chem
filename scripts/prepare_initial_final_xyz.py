#!/usr/bin/env python3
"""
Prepare only initial.xyz and final.xyz (no optimization, no NEB, no ORCA).

Single-purpose script: builds Zn²⁺@surface at start/end distances and writes
exactly two XYZ files with charge and multiplicity in the comment line.

Usage:
  python scripts/prepare_initial_final_xyz.py -o orca_inputs
  python scripts/prepare_initial_final_xyz.py --size 2x2 --start 5.0 --end 2.5 -o out/

Output:
  <output_dir>/initial.xyz   - reactant (Zn²⁺ at start_distance)
  <output_dir>/final.xyz     - product (Zn²⁺ at end_distance)
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

_script_dir = Path(__file__).resolve().parent
_repo_root = _script_dir.parent
_src = _repo_root / "src"
if _src.exists() and str(_src) not in sys.path:
    sys.path.insert(0, str(_src))
if str(_script_dir) not in sys.path:
    sys.path.insert(0, str(_script_dir))

from prepare_endpoints import run as prepare_run


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Prepare only initial.xyz and final.xyz (no optimization)."
    )
    parser.add_argument(
        "-o", "--output-dir",
        type=Path,
        default=Path("orca_inputs"),
        help="Output directory for initial.xyz and final.xyz (default: orca_inputs)",
    )
    parser.add_argument(
        "--size",
        type=str,
        default="2x2",
        metavar="NxM",
        help="Surface supercell size (default: 2x2)",
    )
    parser.add_argument(
        "--start",
        "--start-distance",
        type=float,
        default=5.0,
        dest="start_distance",
        metavar="Å",
        help="Zn²⁺ distance from surface for initial endpoint (Å, default: 5.0)",
    )
    parser.add_argument(
        "--end",
        "--end-distance",
        type=float,
        default=2.5,
        dest="end_distance",
        metavar="Å",
        help="Zn²⁺ distance from surface for final endpoint (Å, default: 2.5)",
    )
    parser.add_argument(
        "--carboxyl",
        type=int,
        default=0,
        metavar="N",
        help="Number of carboxyl groups (default: 0)",
    )
    parser.add_argument(
        "--hydroxyl",
        type=int,
        default=0,
        metavar="N",
        help="Number of hydroxyl groups (default: 0)",
    )
    parser.add_argument(
        "--charge",
        type=int,
        default=2,
        help="Total charge (default: 2)",
    )
    parser.add_argument(
        "--multiplicity",
        type=int,
        default=1,
        choices=(1, 2, 3, 4),
        help="Spin multiplicity (default: 1). Ignored if --auto-multiplicity is set.",
    )
    parser.add_argument(
        "--auto-multiplicity",
        action="store_true",
        help="Set multiplicity for closed-shell when possible (even electrons → 1, odd → 2).",
    )
    args = parser.parse_args()

    try:
        nx, ny = map(int, args.size.lower().split("x"))
        size = (nx, ny)
    except Exception:
        parser.error("--size must be of the form NxM (e.g. 2x2, 4x4)")

    multiplicity: int | None = None if args.auto_multiplicity else args.multiplicity
    return prepare_run(
        args.output_dir,
        size=size,
        start_distance=args.start_distance,
        end_distance=args.end_distance,
        carboxyl=args.carboxyl,
        hydroxyl=args.hydroxyl,
        charge=args.charge,
        multiplicity=multiplicity,
    )


if __name__ == "__main__":
    sys.exit(main())
