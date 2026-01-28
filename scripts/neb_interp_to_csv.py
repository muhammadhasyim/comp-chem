#!/usr/bin/env python3
"""
Extract the latest iteration from an ORCA NEB .interp file and write it to CSV.

Usage:
  python scripts/neb_interp_to_csv.py [interp_file] [-o output.csv]
  python scripts/neb_interp_to_csv.py orca_inputs/neb.interp -o neb_latest.csv

Output columns: iteration, section, fraction, distance_bohr, energy_eh
  - section is "images" or "interp"
  - fraction: reaction coordinate 0..1
  - distance_bohr: path length (Bohr)
  - energy_eh: energy (Hartree)
"""

import argparse
import csv
import re
import sys
from pathlib import Path


def parse_interp_file(path: Path) -> list[dict]:
    """Parse neb.interp and return list of blocks, each with iteration number and data rows."""
    text = path.read_text()
    blocks = []
    # Split by "Iteration: N" (keep the header in each block)
    parts = re.split(r"\n(?=Iteration:\s*\d+)", text)
    for part in parts:
        part = part.strip()
        if not part:
            continue
        m = re.match(r"Iteration:\s*(\d+)", part)
        if not m:
            continue
        iteration = int(m.group(1))
        rest = part[m.end() :].strip()
        # Parse "Images: ..." section
        images = []
        interp = []
        current = None
        for line in rest.splitlines():
            line = line.strip()
            if not line:
                continue
            if line.startswith("Images:"):
                current = "images"
                continue
            if line.startswith("Interp."):
                current = "interp"
                continue
            # Data line: three floats
            tokens = line.split()
            if len(tokens) >= 3 and current:
                try:
                    frac = float(tokens[0])
                    dist = float(tokens[1])
                    ene = float(tokens[2])
                    row = {
                        "iteration": iteration,
                        "section": current,
                        "fraction": frac,
                        "distance_bohr": dist,
                        "energy_eh": ene,
                    }
                    if current == "images":
                        images.append(row)
                    else:
                        interp.append(row)
                except ValueError:
                    pass
        blocks.append({"iteration": iteration, "images": images, "interp": interp})
    return blocks


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Extract latest NEB iteration from .interp file to CSV"
    )
    parser.add_argument(
        "interp_file",
        nargs="?",
        default="orca_inputs/neb.interp",
        type=Path,
        help="Path to neb.interp file (default: orca_inputs/neb.interp)",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=None,
        help="Output CSV path (default: <interp_stem>_latest.csv)",
    )
    parser.add_argument(
        "--section",
        choices=["interp", "images", "both"],
        default="interp",
        help="Which section to export: interp (dense), images, or both (default: interp)",
    )
    args = parser.parse_args()

    if not args.interp_file.exists():
        print(f"Error: file not found: {args.interp_file}", file=sys.stderr)
        sys.exit(1)

    blocks = parse_interp_file(args.interp_file)
    if not blocks:
        print("Error: no iterations found in file", file=sys.stderr)
        sys.exit(1)

    latest = blocks[-1]
    rows = []
    if args.section in ("interp", "both"):
        rows.extend(latest["interp"])
    if args.section in ("images", "both"):
        rows.extend(latest["images"])

    if not rows:
        print("Error: no data rows in latest iteration", file=sys.stderr)
        sys.exit(1)

    out = args.output or (args.interp_file.parent / f"{args.interp_file.stem}_latest.csv")
    fieldnames = ["iteration", "section", "fraction", "distance_bohr", "energy_eh"]
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)

    print(f"Iteration {latest['iteration']}: {len(rows)} rows -> {out}")


if __name__ == "__main__":
    main()
