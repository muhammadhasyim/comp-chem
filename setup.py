from setuptools import setup, find_packages

setup(
    name="zn2-adsorption",
    version="0.1.0",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    install_requires=[
        "pymatgen>=2024.1.20",
        "ase>=3.22.1",
        "numpy>=1.24.0",
        "scipy>=1.10.0",
        "matplotlib>=3.7.0",
        "ccinput>=1.0.0",
        "rdkit>=2022.9.5",
    ],
    entry_points={
        "console_scripts": [
            "run-adsorption=zn2_adsorption.cli:main",
        ],
    },
    author="Scientific Software Engineer",
    description="Zn2+ Adsorption Energy Calculator CLI",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    python_requires=">=3.9",
)
