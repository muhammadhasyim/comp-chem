"""
Tests for NEB calculator module.
"""

import pytest
from pathlib import Path
import tempfile
import shutil

try:
    from ase import Atoms
    ASE_AVAILABLE = True
except ImportError:
    ASE_AVAILABLE = False
    pytestmark = pytest.mark.skip("ASE not available")

from zn2_adsorption.neb_calculator import NebCalculator
from zn2_adsorption.orca_generator import OrcaInputGenerator
from zn2_adsorption.surface_builder import FunctionalizedGrapheneBuilder


@pytest.fixture
def temp_dir():
    """Create temporary directory for tests."""
    tmpdir = tempfile.mkdtemp()
    yield Path(tmpdir)
    shutil.rmtree(tmpdir)


@pytest.fixture
def neb_calculator():
    """Create NebCalculator instance for testing."""
    orca_gen = OrcaInputGenerator(
        method="B3LYP",
        basis_set="def2-TZVP",
        solvent="water",
        solvation_model="CPCM"
    )
    return NebCalculator(orca_generator=orca_gen)


def test_neb_calculator_initialization():
    """Test NebCalculator initialization."""
    calc = NebCalculator()
    assert calc.orca_gen is not None
    assert calc.surface_builder is not None


def test_create_initial_structure(neb_calculator):
    """Test creation of initial structure."""
    builder = neb_calculator.surface_builder
    pristine = builder.build_pristine_surface(supercell_size=(2, 2))
    surface, _ = builder.add_functional_groups(pristine, num_carboxyl=1, num_hydroxyl=0)
    
    initial = neb_calculator.create_initial_structure(surface, distance=5.0)
    
    # Check that Zn²⁺ was added
    zn_symbols = [site.specie.symbol for site in initial if site.specie.symbol == "Zn"]
    assert len(zn_symbols) == 1
    
    # Check that Zn²⁺ is at approximately correct distance
    zn_site = [site for site in initial if site.specie.symbol == "Zn"][0]
    zn_z = zn_site.coords[2]
    
    # Surface should be around z=0, Zn²⁺ should be around z=5.0
    # Allow some tolerance
    assert abs(zn_z - 5.0) < 1.0


def test_create_final_structure(neb_calculator):
    """Test creation of final structure."""
    builder = neb_calculator.surface_builder
    pristine = builder.build_pristine_surface(supercell_size=(2, 2))
    surface, _ = builder.add_functional_groups(pristine, num_carboxyl=1, num_hydroxyl=0)
    
    final = neb_calculator.create_final_structure(surface, distance=2.5)
    
    # Check that Zn²⁺ was added
    zn_symbols = [site.specie.symbol for site in final if site.specie.symbol == "Zn"]
    assert len(zn_symbols) == 1
    
    # Check that Zn²⁺ is at approximately correct distance
    zn_site = [site for site in final if site.specie.symbol == "Zn"][0]
    zn_z = zn_site.coords[2]
    
    # Surface should be around z=0, Zn²⁺ should be around z=2.5
    assert abs(zn_z - 2.5) < 1.0


def test_ensure_atom_ordering(neb_calculator):
    """Test that atom ordering ensures Zn²⁺ is last."""
    builder = neb_calculator.surface_builder
    pristine = builder.build_pristine_surface(supercell_size=(2, 2))
    surface, _ = builder.add_functional_groups(pristine, num_carboxyl=1, num_hydroxyl=0)
    
    # Add Zn²⁺ (should already be last, but test reordering anyway)
    structure_with_zn = builder.add_zn2_ion(surface, distance=3.0)
    
    # Ensure ordering
    ordered = neb_calculator._ensure_atom_ordering(structure_with_zn)
    
    # Check that Zn²⁺ is last
    last_site = ordered[-1]
    assert last_site.specie.symbol == "Zn"
    
    # Check that all atoms are still present
    assert len(ordered) == len(structure_with_zn)


def test_prepare_scan_base(neb_calculator):
    """Test prepare_scan_base returns surface without Zn and metadata."""
    base = neb_calculator.prepare_scan_base(
        num_carboxyl=1,
        num_hydroxyl=0,
        surface_size=(2, 2),
    )
    assert "surface_structure" in base
    assert "n_graphene" in base
    assert "fg_bond_constraints" in base
    # No Zn in surface
    symbols = [s.specie.symbol for s in base["surface_structure"]]
    assert "Zn" not in symbols
    assert base["n_graphene"] > 0
    assert isinstance(base["fg_bond_constraints"], list)


def test_prepare_neb_calculation(neb_calculator, temp_dir):
    """Test complete NEB calculation preparation."""
    results = neb_calculator.prepare_neb_calculation(
        num_carboxyl=1,
        num_hydroxyl=0,
        surface_size=(2, 2),
        start_distance=5.0,
        end_distance=2.5,
        neb_images=7,
        output_dir=temp_dir
    )
    
    # Check that all required files were created
    assert "initial_input" in results
    assert "final_input" in results
    assert "neb_input" in results
    assert results["initial_xyz"] is not None
    assert results["final_xyz"] is not None
    
    # Check that files exist
    assert (temp_dir / "initial.inp").exists()
    assert (temp_dir / "final.inp").exists()
    assert (temp_dir / "neb.inp").exists()
    assert (temp_dir / "initial.xyz").exists()
    assert (temp_dir / "final.xyz").exists()
    
    # Check that NEB input contains required keywords
    neb_content = (temp_dir / "neb.inp").read_text()
    assert "NEB-TS" in neb_content
    assert "NEB_END_XYZFILE" in neb_content
    assert "PREOPT_ENDS" in neb_content
    assert "TRUE" in neb_content or "true" in neb_content.lower()

    # With constrain_endpoints=True (default): initial/final.inp should have Constraints
    initial_content = (temp_dir / "initial.inp").read_text()
    assert "Constraints" in initial_content
    assert "{ C 0 C }" in initial_content  # Cartesian freeze for graphene carbons (0-based)
    assert "{ B " in initial_content  # Bond constraint for Zn-surface distance


def test_prepare_endpoints_only(neb_calculator, temp_dir):
    """Test prepare_endpoints_only produces initial.inp and final.inp only (no neb.inp)."""
    results = neb_calculator.prepare_endpoints_only(
        num_carboxyl=1,
        num_hydroxyl=0,
        surface_size=(2, 2),
        start_distance=5.0,
        end_distance=2.5,
        output_dir=temp_dir,
    )
    assert "initial_input" in results
    assert "final_input" in results
    assert "initial_xyz" in results
    assert "final_xyz" in results
    assert "neb_input" not in results

    assert (temp_dir / "initial.inp").exists()
    assert (temp_dir / "final.inp").exists()
    assert (temp_dir / "initial.xyz").exists()
    assert (temp_dir / "final.xyz").exists()
    assert not (temp_dir / "neb.inp").exists()

    # Constrained: surface + Zn frozen, FG-surface bonds fixed (1 carboxyl)
    initial_content = (temp_dir / "initial.inp").read_text()
    assert "Constraints" in initial_content
    assert "{ C 0 C }" in initial_content
    assert "{ B " in initial_content  # FG-surface bond constraint(s)


def test_prepare_endpoints_only_uff_mode(neb_calculator, temp_dir):
    """Test prepare_endpoints_only with generate_orca_inputs=False (UFF mode)."""
    results = neb_calculator.prepare_endpoints_only(
        num_carboxyl=0,
        num_hydroxyl=0,
        surface_size=(2, 2),
        start_distance=5.0,
        end_distance=2.5,
        output_dir=temp_dir,
        generate_orca_inputs=False,
    )
    assert "initial_xyz" in results
    assert "final_xyz" in results
    assert "initial_ase" in results
    assert "final_ase" in results
    assert "n_graphene" in results
    assert "initial_input" not in results
    assert "final_input" not in results

    assert (temp_dir / "initial.xyz").exists()
    assert (temp_dir / "final.xyz").exists()
    assert not (temp_dir / "initial.inp").exists()
    assert not (temp_dir / "final.inp").exists()


def test_prepare_endpoints_only_unconstrained(neb_calculator, temp_dir):
    """Test prepare_endpoints_only with constrain_endpoints=False."""
    neb_calculator.prepare_endpoints_only(
        num_carboxyl=0,
        num_hydroxyl=0,
        surface_size=(2, 2),
        start_distance=5.0,
        end_distance=2.5,
        output_dir=temp_dir,
        constrain_endpoints=False,
    )
    assert not (temp_dir / "neb.inp").exists()
    initial_content = (temp_dir / "initial.inp").read_text()
    assert "Constraints" not in initial_content


def test_prepare_neb_calculation_unconstrained(neb_calculator, temp_dir):
    """Test NEB preparation with constrain_endpoints=False."""
    results = neb_calculator.prepare_neb_calculation(
        num_carboxyl=1,
        num_hydroxyl=0,
        surface_size=(2, 2),
        start_distance=5.0,
        end_distance=2.5,
        output_dir=temp_dir,
        constrain_endpoints=False,
    )
    # Unconstrained: no %geom Constraints block in initial.inp
    initial_content = (temp_dir / "initial.inp").read_text()
    assert "Constraints" not in initial_content


def test_calculate_multiplicity(neb_calculator):
    """Test multiplicity calculation."""
    # Test with even number of electrons
    atoms_even = Atoms(symbols=['C', 'H', 'H', 'H', 'H'], positions=[[0,0,0]]*5)
    mult = neb_calculator._calculate_multiplicity(atoms_even, charge=0)
    assert mult == 1  # Even electrons → singlet
    
    # Test with odd number of electrons
    atoms_odd = Atoms(symbols=['C', 'H', 'H', 'H'], positions=[[0,0,0]]*4)
    mult = neb_calculator._calculate_multiplicity(atoms_odd, charge=0)
    assert mult == 2  # Odd electrons → doublet
    
    # Test with charge
    atoms_charged = Atoms(symbols=['Zn'], positions=[[0,0,0]])
    mult = neb_calculator._calculate_multiplicity(atoms_charged, charge=2)
    # Zn has 30 electrons, -2 charge = 28 electrons (even) → singlet
    assert mult == 1


def test_export_xyz(neb_calculator, temp_dir):
    """Test XYZ file export."""
    atoms = Atoms(symbols=['C', 'H', 'H', 'H', 'H'], positions=[[0,0,0]]*5)
    output_file = temp_dir / "test.xyz"
    
    neb_calculator._export_xyz(atoms, output_file)
    
    assert output_file.exists()
    content = output_file.read_text()
    assert "5" in content  # Number of atoms
    assert "C" in content
    assert "H" in content


def test_prepare_neb_calculation_atom_ordering(neb_calculator, temp_dir):
    """Test that initial and final structures have consistent atom ordering."""
    results = neb_calculator.prepare_neb_calculation(
        num_carboxyl=0,
        num_hydroxyl=0,
        surface_size=(2, 2),
        start_distance=5.0,
        end_distance=2.5,
        output_dir=temp_dir
    )
    
    # Read XYZ files and check atom ordering
    from ase.io import read
    
    initial_xyz = read(str(temp_dir / "initial.xyz"))
    final_xyz = read(str(temp_dir / "final.xyz"))
    
    # Check that both have same number of atoms
    assert len(initial_xyz) == len(final_xyz)
    
    # Check that atom symbols match
    initial_symbols = initial_xyz.get_chemical_symbols()
    final_symbols = final_xyz.get_chemical_symbols()
    assert initial_symbols == final_symbols
    
    # Check that Zn²⁺ is last in both
    assert initial_symbols[-1] == "Zn"
    assert final_symbols[-1] == "Zn"
