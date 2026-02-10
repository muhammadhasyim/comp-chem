#!/usr/bin/env python3
"""
Plot FSM path: arc length (s) vs relative energy from vfile_*.xyz.

Reads the latest vfile_*.xyz in an FSM output directory (or a specific file),
plots s [Å] vs E_rel [eV], and marks the TS guess (max energy).

Usage:
  python scripts/plot_fsm_path.py orca_inputs/fsm_reaction/fsm_interp_..._run1
  python scripts/plot_fsm_path.py orca_inputs/fsm_reaction/fsm_interp_.../vfile_03.xyz -o path.png
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def parse_vfile(path: Path) -> list[tuple[float, float]]:
    """Return list of (arc_length_Ang, rel_energy_eV) from vfile multi-frame XYZ."""
    lines = path.read_text().strip().splitlines()
    data = []
    i = 0
    n = len(lines)
    while i < n:
        try:
            natoms = int(lines[i])
        except ValueError:
            i += 1
            continue
        i += 1
        if i >= n:
            break
        parts = lines[i].split()
        if len(parts) >= 2:
            s = float(parts[0])
            e = float(parts[1])
            data.append((s, e))
        i += 1 + natoms
    return data


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot FSM path: s vs E_rel from vfile")
    parser.add_argument(
        "path",
        type=Path,
        help="FSM output directory or path to a vfile_*.xyz",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=None,
        help="Save figure to this path (e.g. path.png or path.pdf)",
    )
    parser.add_argument(
        "--no-show",
        action="store_true",
        help="Do not open interactive plot window",
    )
    parser.add_argument(
        "--title",
        type=str,
        default=None,
        help="Plot title (default: filename or dir name)",
    )
    args = parser.parse_args()
    path = args.path.resolve()

    if path.is_dir():
        vfiles = sorted(path.glob("vfile_*.xyz"))
        if not vfiles:
            raise SystemExit(f"No vfile_*.xyz in {path}")
        vfile = vfiles[-1]
        title = args.title or path.name
    elif path.is_file():
        vfile = path
        title = args.title or vfile.name
    else:
        raise SystemExit(f"Not found: {path}")

    data = parse_vfile(vfile)
    if not data:
        raise SystemExit(f"No frames in {vfile}")

    s = np.array([d[0] for d in data])
    e = np.array([d[1] for d in data])
    ts_idx = int(np.argmax(e))

    fig, ax = plt.subplots()
    ax.plot(s, e, "o-", color="C0", markersize=8, label="FSM path")
    ax.plot(s[ts_idx], e[ts_idx], "s", color="C3", markersize=12, label="TS guess")
    ax.set_xlabel("Arc length $s$ (Å)")
    ax.set_ylabel("Relative energy (eV)")
    ax.set_title(title)
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()

    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(args.output, dpi=150)
        print(f"Saved {args.output}")

    if not args.no_show:
        plt.show()


if __name__ == "__main__":
    main()
