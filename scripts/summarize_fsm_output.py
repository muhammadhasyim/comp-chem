#!/usr/bin/env python3
"""
Summarize FSM output: read the latest vfile_*.xyz in an FSM output directory
and print path summary + TS guess (frame with max relative energy).

Usage:
  python scripts/summarize_fsm_output.py orca_inputs/fsm_reaction/fsm_interp_..._run1
"""

import argparse
from pathlib import Path


def parse_vfile(path: Path) -> list[tuple[float, float]]:
    """Return list of (arc_length, rel_energy_eV) from vfile multi-frame XYZ."""
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
    parser = argparse.ArgumentParser(description="Summarize FSM vfile output")
    parser.add_argument("outdir", type=Path, help="FSM output directory (e.g. fsm_interp_..._run1)")
    args = parser.parse_args()
    outdir = args.outdir.resolve()
    if not outdir.is_dir():
        raise SystemExit(f"Not a directory: {outdir}")

    vfiles = sorted(outdir.glob("vfile_*.xyz"))
    if not vfiles:
        print(f"No vfile_*.xyz in {outdir}")
        return
    latest = vfiles[-1]
    data = parse_vfile(latest)
    if not data:
        print(f"No frames in {latest}")
        return

    print(f"FSM output: {outdir.name}")
    print(f"Latest vfile: {latest.name} ({len(data)} nodes)")
    print()
    print("Path (arc length s [Å], relative energy [eV]):")
    for i, (s, e) in enumerate(data):
        print(f"  node {i + 1}: s = {s:.4f}   E_rel = {e:.4f}")
    ts_idx = max(range(len(data)), key=lambda i: data[i][1])
    s_ts, e_ts = data[ts_idx]
    print()
    print(f"TS guess: node {ts_idx + 1}  (s = {s_ts:.4f} Å, E_rel = {e_ts:.4f} eV)")
    ngrad = outdir / "ngrad.txt"
    if ngrad.exists():
        n = ngrad.read_text().strip()
        print(f"Gradient calls (ngrad.txt): {n}")


if __name__ == "__main__":
    main()
