import pytest
import random
import numpy as np
from pymatgen.core import Structure
from zn2_adsorption.surface_builder import FunctionalizedGrapheneBuilder


def test_builder_pristine():
    builder = FunctionalizedGrapheneBuilder()
    structure = builder.build_pristine_surface(supercell_size=(2, 2))
    assert isinstance(structure, Structure)
    assert len(structure) == 8 # 2 atoms per unit cell * 2 * 2 supercell
    assert structure.lattice is not None

def test_builder_functionalized():
    builder = FunctionalizedGrapheneBuilder()
    structure = builder.build_pristine_surface(supercell_size=(4, 4))
    functionalized, fg_bonds = builder.add_functional_groups(
        structure,
        num_carboxyl=1,
        num_hydroxyl=1
    )
    assert len(functionalized) > len(structure)
    symbols = [s.symbol for s in functionalized.species]
    assert 'O' in symbols
    assert 'H' in symbols

def test_add_zn2():
    builder = FunctionalizedGrapheneBuilder()
    structure = builder.build_pristine_surface(supercell_size=(4, 4))
    zn_structure = builder.add_zn2_ion(structure, distance=3.0)
    assert 'Zn' in [s.symbol for s in zn_structure.species]
    
    # Check distance
    zn_idx = [i for i, s in enumerate(zn_structure.species) if s.symbol == 'Zn'][0]
    zn_pos = zn_structure.cart_coords[zn_idx]
    
    # Find surface z
    c_coords = np.array([s.coords for s in zn_structure if s.specie.symbol == "C"])
    surface_z = np.mean(c_coords[:, 2])
    
    assert abs(zn_pos[2] - surface_z - 3.0) < 0.1


def test_functional_groups_distinct_sites():
    """Verify hydroxyl and carboxyl occupy distinct carbon sites (no overlap)."""
    builder = FunctionalizedGrapheneBuilder()
    structure = builder.build_pristine_surface(supercell_size=(4, 4))
    n_graphene = len(structure)

    # With fixed seed, add 2 carboxylate + 2 hydroxyl (4 distinct sites)
    random.seed(42)
    functionalized, _ = builder.add_functional_groups(
        structure,
        num_carboxyl=2,
        num_hydroxyl=2,
    )

    # Expected atom count: 4 groups on distinct sites
    # Each -OH adds 1 O + 1 H; each -COO⁻ adds 1 C + 2 O (no H)
    expected_extra = 2 * (1 + 1) + 2 * (1 + 2)  # 4 + 6 = 10
    assert len(functionalized) == n_graphene + expected_extra

    # Run with multiple seeds to exercise no-overlap logic
    for seed in [0, 1, 99, 12345]:
        random.seed(seed)
        f, _ = builder.add_functional_groups(
            structure.copy(),
            num_carboxyl=2,
            num_hydroxyl=2,
        )
        assert len(f) == n_graphene + expected_extra


def test_functional_groups_no_nearest_neighbors():
    """Verify FGs are not placed on nearest-neighbor carbons (uniform placement)."""
    builder = FunctionalizedGrapheneBuilder()
    structure = builder.build_pristine_surface(supercell_size=(4, 4))
    n_graphene = len(structure)
    C_C = 1.5  # NN cutoff (slightly above 1.42)

    # Add several FGs with fixed seed
    random.seed(123)
    functionalized, fg_bonds = builder.add_functional_groups(
        structure, num_carboxyl=3, num_hydroxyl=2
    )

    # fg_bonds: (anchor_idx, fg_first_idx, bond_len) - anchor_idx are the surface C indices
    anchor_indices = [b[0] for b in fg_bonds]
    assert len(anchor_indices) == 5

    # Verify no pair of anchors are nearest neighbors
    for i, a in enumerate(anchor_indices):
        for b in anchor_indices[i + 1 :]:
            dist = structure.get_distance(a, b)
            assert dist >= C_C, (
                f"FGs on sites {a} and {b} are nearest neighbors (d={dist:.3f} < {C_C})"
            )


def test_functional_groups_clustered_around_zn():
    """Verify FGs are clustered around Zn center, not directly below it."""
    builder = FunctionalizedGrapheneBuilder()
    structure = builder.build_pristine_surface(supercell_size=(4, 4))
    n_graphene = len(structure)
    # Zn center = surface center (same as add_zn2_ion)
    c_coords = np.array([s.coords for s in structure if s.specie.symbol == "C"])
    zn_center = (float(np.mean(c_coords[:, 0])), float(np.mean(c_coords[:, 1])))

    functionalized, fg_bonds = builder.add_functional_groups(
        structure,
        num_carboxyl=1,
        num_hydroxyl=2,
        zn_center_xy=zn_center,
    )

    r_exclude = builder._r_exclude
    r_cluster = builder._r_cluster
    cx, cy = zn_center
    for anchor_idx, _, _ in fg_bonds:
        pos = structure.cart_coords[anchor_idx]
        r = np.sqrt((pos[0] - cx) ** 2 + (pos[1] - cy) ** 2)
        assert r >= r_exclude, f"Anchor {anchor_idx} in exclusion zone (r={r:.3f} < {r_exclude})"
        assert r <= r_cluster, f"Anchor {anchor_idx} outside cluster zone (r={r:.3f} > {r_cluster})"


def test_find_central_carbon_away_from_fg():
    """Central carbon is not an FG anchor or its neighbor."""
    builder = FunctionalizedGrapheneBuilder()
    structure = builder.build_pristine_surface(supercell_size=(4, 4))
    n_graphene = len(structure)
    c_coords = np.array([s.coords for s in structure if s.specie.symbol == "C"])
    center_xy = tuple(np.mean(c_coords[:, :2], axis=0))

    random.seed(42)
    functionalized, fg_bonds = builder.add_functional_groups(
        structure, num_carboxyl=2, num_hydroxyl=1,
    )
    fg_anchors = [i for i, _, _ in fg_bonds]
    forbidden = set(fg_anchors)
    for a in fg_anchors:
        forbidden |= builder._get_neighbor_indices(
            functionalized, {a}, n_graphene
        )

    target = builder.find_central_carbon_away_from_fg(
        functionalized, center_xy, n_graphene, fg_anchor_indices=fg_anchors,
    )
    assert target not in forbidden


def test_functional_groups_insufficient_sites_raises():
    """Requesting more groups than graphene carbons raises ValueError."""
    builder = FunctionalizedGrapheneBuilder()
    structure = builder.build_pristine_surface(supercell_size=(4, 4))
    n_carbons = sum(1 for s in structure.species if s.symbol == "C")

    with pytest.raises(ValueError, match="Not enough distinct carbon sites"):
        builder.add_functional_groups(
            structure,
            num_carboxyl=n_carbons,  # uses all sites
            num_hydroxyl=1,          # no site left
        )
