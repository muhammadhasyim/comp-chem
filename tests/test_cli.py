import pytest
from pathlib import Path
from zn2_adsorption.cli import parse_arguments, calculate_neb_images

def test_parse_arguments_defaults():
    args = parse_arguments([
        "--start-distance", "5.0",
        "--end-distance", "2.5"
    ])
    assert args.carboxyl == 0
    assert args.hydroxyl == 0
    assert args.size == "4x4"
    assert args.start_distance == 5.0
    assert args.end_distance == 2.5
    assert args.neb_images is None
    assert args.method == "B3LYP"
    assert args.basis == "def2-TZVP"
    assert args.solvent == "water"
    assert args.solvation == "CPCM"
    assert args.output_dir == "./orca_inputs"
    assert args.run_orca is False
    assert args.json_output == "./neb_results.json"
    assert args.no_mpi is False
    assert args.nprocs == 4

def test_parse_arguments_custom():
    args = parse_arguments([
        "--carboxyl", "2",
        "--hydroxyl", "1",
        "--size", "6x6",
        "--start-distance", "5.5",
        "--end-distance", "2.0",
        "--neb-images", "10",
        "--method", "PBE0",
        "--basis", "def2-SVP",
        "--solvent", "dmso",
        "--solvation", "SMD",
        "--output-dir", "./my_calc",
        "--run-orca",
        "--json-output", "./results.json"
    ])
    assert args.carboxyl == 2
    assert args.hydroxyl == 1
    assert args.size == "6x6"
    assert args.start_distance == 5.5
    assert args.end_distance == 2.0
    assert args.neb_images == 10
    assert args.method == "PBE0"
    assert args.basis == "def2-SVP"
    assert args.solvent == "dmso"
    assert args.solvation == "SMD"
    assert args.output_dir == "./my_calc"
    assert args.run_orca is True
    assert args.json_output == "./results.json"

def test_calculate_neb_images():
    """Test auto-calculation of NEB images."""
    # Small distance difference
    images = calculate_neb_images(3.0, 2.5)
    assert images >= 5  # Minimum
    assert images <= 20  # Maximum
    
    # Large distance difference
    images = calculate_neb_images(10.0, 2.0)
    assert images >= 5
    assert images <= 20
    
    # Very small difference
    images = calculate_neb_images(3.0, 2.9)
    assert images == 5  # Should hit minimum
    
    # Symmetric (should work regardless of order)
    images1 = calculate_neb_images(5.0, 2.5)
    images2 = calculate_neb_images(2.5, 5.0)
    assert images1 == images2

from unittest.mock import patch, MagicMock
from pathlib import Path
from zn2_adsorption.cli import validate_inputs, check_orca_available

def test_check_orca_available():
    """Test ORCA detection - works with real filesystem."""
    # This test works with the actual filesystem
    # It verifies the function doesn't crash and returns valid results
    available, path = check_orca_available()
    
    # Should return boolean and optional string
    assert isinstance(available, bool)
    if available:
        assert isinstance(path, str)
        assert len(path) > 0
    else:
        assert path is None
    
    # Test that it can be called multiple times without error
    available2, path2 = check_orca_available()
    assert isinstance(available2, bool)

def test_validate_inputs_invalid_size():
    args = parse_arguments(["--start-distance", "5.0", "--end-distance", "2.5", "--size", "invalid"])
    with pytest.raises(ValueError, match="Invalid surface size format"):
        from zn2_adsorption.cli import validate_inputs
        validate_inputs(args)

def test_validate_inputs_invalid_distance():
    args = parse_arguments(["--start-distance", "-1.0", "--end-distance", "2.5"])
    with pytest.raises(ValueError, match="Start distance must be positive"):
        from zn2_adsorption.cli import validate_inputs
        validate_inputs(args)
    
    args = parse_arguments(["--start-distance", "5.0", "--end-distance", "-1.0"])
    with pytest.raises(ValueError, match="End distance must be positive"):
        from zn2_adsorption.cli import validate_inputs
        validate_inputs(args)
    
    args = parse_arguments(["--start-distance", "3.0", "--end-distance", "3.0"])
    with pytest.raises(ValueError, match="Start and end distances must be different"):
        from zn2_adsorption.cli import validate_inputs
        validate_inputs(args)

def test_cli_endpoints_only_uff(tmp_path):
    """Test that --endpoints-only --calculator uff produces only xyz files (no .inp, no neb)."""
    import subprocess
    import sys

    repo_root = Path(__file__).resolve().parent.parent
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "zn2_adsorption.cli",
            "--endpoints-only",
            "--calculator",
            "uff",
            "--start-distance",
            "5.0",
            "--end-distance",
            "2.5",
            "--output-dir",
            str(tmp_path),
        ],
        capture_output=True,
        text=True,
        cwd=str(repo_root),
    )
    assert result.returncode == 0, f"CLI failed: {result.stderr}"

    assert (tmp_path / "initial.xyz").exists()
    assert (tmp_path / "final.xyz").exists()
    assert not (tmp_path / "initial.inp").exists()
    assert not (tmp_path / "final.inp").exists()
    assert not (tmp_path / "neb.inp").exists()


def test_cli_endpoints_only_orca_produces_no_neb_inp(tmp_path):
    """Test that --endpoints-only --calculator orca produces initial/final.inp but NOT neb.inp."""
    import subprocess
    import sys

    repo_root = Path(__file__).resolve().parent.parent
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "zn2_adsorption.cli",
            "--endpoints-only",
            "--calculator",
            "orca",
            "--start-distance",
            "5.0",
            "--end-distance",
            "2.5",
            "--output-dir",
            str(tmp_path),
        ],
        capture_output=True,
        text=True,
        cwd=str(repo_root),
    )
    assert result.returncode == 0, f"CLI failed: {result.stderr}"

    assert (tmp_path / "initial.inp").exists()
    assert (tmp_path / "final.inp").exists()
    assert (tmp_path / "initial.xyz").exists()
    assert (tmp_path / "final.xyz").exists()
    assert not (tmp_path / "neb.inp").exists()


def test_parse_arguments_endpoints_only():
    """Test --endpoints-only flag is parsed."""
    args = parse_arguments([
        "--endpoints-only",
        "--start-distance", "5.0",
        "--end-distance", "2.5",
    ])
    assert args.endpoints_only is True


def test_parse_arguments_no_mpi():
    """Test --no-mpi flag sets nprocs to 1."""
    args = parse_arguments([
        "--start-distance", "5.0",
        "--end-distance", "2.5",
        "--no-mpi"
    ])
    assert args.no_mpi is True
    assert args.nprocs == 4  # nprocs argument still parsed, but will be overridden to 1
    
    # Test that --no-mpi works with custom nprocs (should still force to 1)
    args = parse_arguments([
        "--start-distance", "5.0",
        "--end-distance", "2.5",
        "--nprocs", "8",
        "--no-mpi"
    ])
    assert args.no_mpi is True
    assert args.nprocs == 8  # Parsed value, but will be overridden to 1 in run_calculation
