import pytest
from pathlib import Path
from zn2_adsorption.adsorption_calculator import AdsorptionEnergyCalculator
from zn2_adsorption.surface_builder import FunctionalizedGrapheneBuilder

def test_calculator_init():
    calc = AdsorptionEnergyCalculator()
    assert calc.orca_gen is not None
    assert isinstance(calc.surface_builder, FunctionalizedGrapheneBuilder)

def test_calculator_workflow(tmp_path):
    calc = AdsorptionEnergyCalculator()
    output_dir = tmp_path / "orca_files"
    
    result = calc.prepare_calculations(
        num_carboxyl=1,
        num_hydroxyl=1,
        surface_size=(2, 2),
        zn2_distance=3.0,
        constrain_distance=True,
        output_dir=output_dir
    )
    
    assert (output_dir / "zn2_surface.inp").exists()
    assert (output_dir / "surface.inp").exists()
    assert (output_dir / "zn2_ion.inp").exists()
    assert "zn2_surface_input" in result
    assert "%geom" in result["zn2_surface_input"]
