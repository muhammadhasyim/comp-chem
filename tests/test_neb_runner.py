"""Tests for NEB runner module."""

import pytest

try:
    from ase import Atoms
    from zn2_adsorption.neb_runner import run_neb_uff
    ASE_AVAILABLE = True
except ImportError:
    ASE_AVAILABLE = False
    pytestmark = pytest.mark.skip("ASE not available")

try:
    from openbabel import openbabel
    OPENBABEL_AVAILABLE = True
except ImportError:
    OPENBABEL_AVAILABLE = False


@pytest.mark.skipif(not OPENBABEL_AVAILABLE, reason="Open Babel not installed")
def test_run_neb_uff_basic():
    """Test NEB runner produces images and energies."""
    initial = Atoms("H2", positions=[[0, 0, 0], [0, 0, 0.74]])
    final = Atoms("H2", positions=[[0, 0, 0], [0, 0, 1.2]])
    result = run_neb_uff(
        initial,
        final,
        n_images=4,
        fmax=0.5,
        max_steps=10,
    )
    assert "images" in result
    assert "energies" in result
    assert "ts_index" in result
    assert "converged" in result
    assert len(result["images"]) == 4
    assert len(result["energies"]) == 4
    assert 0 <= result["ts_index"] < 4


@pytest.mark.skipif(not OPENBABEL_AVAILABLE, reason="Open Babel not installed")
def test_run_neb_uff_output_dir(tmp_path):
    """Test NEB runner writes XYZ files when output_dir given."""
    initial = Atoms("H2", positions=[[0, 0, 0], [0, 0, 0.74]])
    final = Atoms("H2", positions=[[0, 0, 0], [0, 0, 1.2]])
    run_neb_uff(
        initial,
        final,
        n_images=3,
        fmax=0.5,
        max_steps=5,
        output_dir=tmp_path,
    )
    assert (tmp_path / "neb_000.xyz").exists()
    assert (tmp_path / "neb_001.xyz").exists()
    assert (tmp_path / "neb_002.xyz").exists()
