#!/usr/bin/env python3
"""
Prepare initial and final NEB endpoint structures (no optimization).

Builds the Zn²⁺@surface initial and final geometries at specified distances,
ensures consistent atom ordering (Zn last), and writes XYZ files with charge
and multiplicity in the comment line for use with ORCA or other codes.

Usage:
  python scripts/prepare_endpoints.py -o orca_inputs
  python scripts/prepare_endpoints.py --size 2x2 --start 5.0 --end 2.5 -o out/
  python scripts/prepare_endpoints.py --charge 2 --multiplicity 1 -o xyz/
  python scripts/prepare_endpoints.py -o orca_inputs --auto-multiplicity  # closed-shell when possible

Output:
  <output_dir>/initial.xyz   - reactant (Zn²⁺ at start_distance)
  <output_dir>/final.xyz     - product (Zn²⁺ at end_distance)

XYZ format: standard XYZ with line 2 = "charge=2 multiplicity=1" (or as set).
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

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
    """
    Return multiplicity for closed-shell when possible (even electrons → 1, odd → 2).

    Total electrons = sum(atomic numbers) - charge. Even → singlet, odd → doublet.
    """
    if Element is None:
        raise RuntimeError("pymatgen is required for --auto-multiplicity. Install with: pip install pymatgen")
    n_electrons = sum(Element(s).Z for s in symbols) - charge
    if n_electrons < 0:
        raise ValueError(f"Total electrons would be negative (charge={charge} too large for composition)")
    return 1 if n_electrons % 2 == 0 else 2


def write_xyz_with_charge(
    output_path: Path,
    symbols: list[str],
    positions: list[tuple[float, float, float]],
    charge: int = 2,
    multiplicity: int = 1,
    comment: str | None = None,
) -> None:
    """
    Write XYZ file with charge and multiplicity in the comment line (line 2).

    Standard XYZ: line1 = natoms, line2 = comment, then "Symbol x y z" per atom.
    Comment format: "charge=2 multiplicity=1" plus optional comment text.
    """
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
    carboxyl: int = 0,
    hydroxyl: int = 0,
    charge: int = 2,
    multiplicity: int | None = 1,
) -> int:
    """
    Prepare only initial.xyz and final.xyz in output_dir (no optimization).

    If multiplicity is None, it is set for closed-shell when possible (even
    electrons → 1, odd electrons → 2). Returns 0 on success, 1 on error.
    """
    if NebCalculator is None:
        print("Error: zn2_adsorption is not installed. Install the package or run from repo root.", file=sys.stderr)
        return 1

    calc = NebCalculator()
    pristine = calc.surface_builder.build_pristine_surface(supercell_size=size)
    surface = calc.surface_builder.add_functional_groups(
        pristine,
        num_carboxyl=carboxyl,
        num_hydroxyl=hydroxyl,
    )
    initial_pmg = calc.create_initial_structure(surface.copy(), start_distance)
    final_pmg = calc.create_final_structure(surface.copy(), end_distance)
    initial_pmg = calc._ensure_atom_ordering(initial_pmg)
    final_pmg = calc._ensure_atom_ordering(final_pmg)

    symbols_initial = [site.specie.symbol for site in initial_pmg]
    positions_initial = [tuple(map(float, site.coords)) for site in initial_pmg]
    symbols_final = [site.specie.symbol for site in final_pmg]
    positions_final = [tuple(map(float, site.coords)) for site in final_pmg]

    if multiplicity is None:
        try:
            multiplicity = closed_shell_multiplicity(symbols_initial, charge)
        except (ValueError, RuntimeError) as e:
            print(f"Error: {e}", file=sys.stderr)
            return 1
        print(f"Auto multiplicity (closed-shell when possible): {multiplicity}")

    write_xyz_with_charge(
        output_dir / "initial.xyz",
        symbols_initial,
        positions_initial,
        charge=charge,
        multiplicity=multiplicity,
        comment="initial endpoint (reactant)",
    )
    write_xyz_with_charge(
        output_dir / "final.xyz",
        symbols_final,
        positions_final,
        charge=charge,
        multiplicity=multiplicity,
        comment="final endpoint (product)",
    )

    print(f"Wrote {output_dir / 'initial.xyz'} (charge={charge}, multiplicity={multiplicity})")
    print(f"Wrote {output_dir / 'final.xyz'}   (charge={charge}, multiplicity={multiplicity})")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Prepare initial and final NEB endpoint XYZ files (no optimization)."
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
        help="Set multiplicity for closed-shell when possible (even electrons → 1, odd → 2).",
    )
    args = parser.parse_args()

    try:
        nx, ny = map(int, args.size.lower().split("x"))
        size = (nx, ny)
    except Exception:
        parser.error("--size must be of the form NxM (e.g. 2x2, 4x4)")

    multiplicity: int | None = None if args.auto_multiplicity else args.multiplicity
    return run(
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
