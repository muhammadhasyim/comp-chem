#!/usr/bin/env python3
"""
Prepare a series of XYZ files along the Zn²⁺ adsorption path (no optimization).

Creates N structures with Zn²⁺ at evenly spaced distances from --start to --end.
No geometry optimization—only the initial geometries at each distance step.

Usage:
  python scripts/prepare_path_xyz.py -o path_xyz
  python scripts/prepare_path_xyz.py --size 2x2 --start 5.0 --end 2.5 -o path_xyz
  python scripts/prepare_path_xyz.py --n-images 15 -o path_xyz

Output:
  <output_dir>/image_000.xyz  through image_<N-1>.xyz

XYZ format: standard XYZ with line 2 = "charge=2 multiplicity=1" (or as set).
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

# Allow running from repo root when package is not installed (src layout)
_script_dir = Path(__file__).resolve().parent
_repo_root = _script_dir.parent
_src = _repo_root / "src"
if _src.exists() and str(_src) not in sys.path:
    sys.path.insert(0, str(_src))

try:
    from zn2_adsorption.neb_calculator import NebCalculator
except ImportError:
    NebCalculator = None

try:
    from pymatgen.core import Element
except ImportError:
    Element = None


def closed_shell_multiplicity(symbols: list[str], charge: int) -> int:
    """Return multiplicity for closed-shell when possible."""
    if Element is None:
        raise RuntimeError(
            "pymatgen is required for --auto-multiplicity. Install with: pip install pymatgen"
        )
    n_electrons = sum(Element(s).Z for s in symbols) - charge
    if n_electrons < 0:
        raise ValueError(
            f"Total electrons would be negative (charge={charge} too large for composition)"
        )
    return 1 if n_electrons % 2 == 0 else 2


def write_xyz_with_charge(
    output_path: Path,
    symbols: list[str],
    positions: list[tuple[float, float, float]],
    charge: int = 2,
    multiplicity: int = 1,
    comment: str | None = None,
) -> None:
    """Write XYZ file with charge and multiplicity in the comment line (line 2)."""
    n = len(symbols)
    if comment is None:
        comment = f"charge={charge} multiplicity={multiplicity}"
    else:
        comment = f"charge={charge} multiplicity={multiplicity}  {comment}".strip()

    lines = [str(n), comment]
    for sym, (x, y, z) in zip(symbols, positions):
        lines.append(f"  {sym:4s} {x:18.12f} {y:18.12f} {z:18.12f}")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(lines) + "\n")


def run(
    output_dir: Path,
    *,
    size: tuple[int, int] = (2, 2),
    start_distance: float = 5.0,
    end_distance: float = 2.5,
    n_images: int = 15,
    carboxyl: int = 0,
    hydroxyl: int = 0,
    charge: int = 2,
    multiplicity: int | None = 1,
) -> int:
    """
    Create n_images XYZ files along the path from start_distance to end_distance.

    Returns 0 on success, 1 on error.
    """
    if NebCalculator is None:
        print(
            "Error: zn2_adsorption is not installed. Install the package or run from repo root.",
            file=sys.stderr,
        )
        return 1

    calc = NebCalculator()
    pristine = calc.surface_builder.build_pristine_surface(supercell_size=size)
    surface, _ = calc.surface_builder.add_functional_groups(
        pristine,
        num_carboxyl=carboxyl,
        num_hydroxyl=hydroxyl,
    )

    # Distances for each image (linear interpolation)
    distances = np.linspace(start_distance, end_distance, n_images)

    # Build first structure to get symbols and set multiplicity
    first_pmg = calc.create_initial_structure(surface.copy(), float(distances[0]))
    first_pmg = calc._ensure_atom_ordering(first_pmg)
    symbols = [site.specie.symbol for site in first_pmg]

    if multiplicity is None:
        try:
            multiplicity = closed_shell_multiplicity(symbols, charge)
        except (ValueError, RuntimeError) as e:
            print(f"Error: {e}", file=sys.stderr)
            return 1
        print(f"Auto multiplicity (closed-shell when possible): {multiplicity}")

    output_dir.mkdir(parents=True, exist_ok=True)
    for i, d in enumerate(distances):
        pmg = calc.create_initial_structure(surface.copy(), float(d))
        pmg = calc._ensure_atom_ordering(pmg)
        positions = [tuple(map(float, site.coords)) for site in pmg]
        out_path = output_dir / f"image_{i:03d}.xyz"
        write_xyz_with_charge(
            out_path,
            symbols,
            positions,
            charge=charge,
            multiplicity=multiplicity,
            comment=f"Zn distance {d:.3f} Å (image {i+1}/{n_images})",
        )
        print(f"  Wrote {out_path} (d_Zn = {d:.3f} Å)")

    print(
        f"\nWrote {n_images} XYZ files to {output_dir}/ "
        f"(charge={charge}, multiplicity={multiplicity})"
    )
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Prepare a series of XYZ files along the Zn²⁺ adsorption path (no optimization)."
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        type=Path,
        default=Path("path_xyz"),
        help="Output directory for image_000.xyz, image_001.xyz, ... (default: path_xyz)",
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
        help="Zn²⁺ distance from surface for first image (Å, default: 5.0)",
    )
    parser.add_argument(
        "--end",
        "--end-distance",
        type=float,
        default=2.5,
        dest="end_distance",
        metavar="Å",
        help="Zn²⁺ distance from surface for last image (Å, default: 2.5)",
    )
    parser.add_argument(
        "--n-images",
        type=int,
        default=15,
        metavar="N",
        help="Number of XYZ files along the path (default: 15)",
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
        help="Total charge of the system (default: 2 for Zn²⁺@surface)",
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
        help="Set multiplicity for closed-shell when possible.",
    )
    args = parser.parse_args()

    try:
        nx, ny = map(int, args.size.lower().split("x"))
        size = (nx, ny)
    except Exception:
        parser.error("--size must be of the form NxM (e.g. 2x2, 4x4)")

    if args.n_images < 1:
        parser.error("--n-images must be at least 1")

    multiplicity: int | None = None if args.auto_multiplicity else args.multiplicity
    return run(
        args.output_dir,
        size=size,
        start_distance=args.start_distance,
        end_distance=args.end_distance,
        n_images=args.n_images,
        carboxyl=args.carboxyl,
        hydroxyl=args.hydroxyl,
        charge=args.charge,
        multiplicity=multiplicity,
    )


if __name__ == "__main__":
    sys.exit(main())
