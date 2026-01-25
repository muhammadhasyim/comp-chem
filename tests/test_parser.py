import pytest
from pathlib import Path
from zn2_adsorption.cli import extract_energy_from_orca

def test_extract_energy_from_orca(tmp_path):
    output_content = """
    ...
    -------------------------
    FINAL SINGLE POINT ENERGY      -123.456789012345
    -------------------------
    ...
    """
    out_file = tmp_path / "test.out"
    out_file.write_text(output_content)
    
    energy = extract_energy_from_orca(out_file)
    assert energy == -123.456789012345

def test_extract_energy_not_found(tmp_path):
    out_file = tmp_path / "empty.out"
    out_file.write_text("No energy here")
    
    energy = extract_energy_from_orca(out_file)
    assert energy is None
