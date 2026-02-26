"""
Plotting utilities for NEB and MEP results.
"""

from __future__ import annotations

from pathlib import Path
from typing import List, Optional

import numpy as np


def plot_neb_energy_profile(
    energies_eV: List[float],
    ts_index: int,
    output_path: Path,
    title: Optional[str] = "NEB Energy Profile",
    ylabel: str = "Energy (eV)",
) -> None:
    """
    Plot NEB energy pathway: image index vs energy.

    Args:
        energies_eV: Energy (eV) for each image.
        ts_index: Index of transition state (max energy).
        output_path: Path to save figure (e.g. path.png).
        title: Plot title.
        ylabel: Y-axis label.
    """
    import matplotlib.pyplot as plt

    n = len(energies_eV)
    x = np.arange(n)
    e = np.array(energies_eV, dtype=float)

    fig, ax = plt.subplots()
    ax.plot(x, e, "o-", color="C0", markersize=6, label="NEB path")
    ax.plot(ts_index, e[ts_index], "s", color="C3", markersize=10, label="TS guess")
    ax.set_xlabel("Image index")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=150)
    plt.close(fig)


def plot_scan_profile(
    distances: List[float],
    energies: List[float],
    output_path: Path,
    title: str = "Distance Scan (UFF)",
    ylabel: str = "Energy (eV)",
    xlabel: str = "Zn-surface distance (Å)",
) -> None:
    """
    Plot E(d) for distance scan: energy vs Zn-surface distance.

    Args:
        distances: Zn-surface distances (Å).
        energies: Optimized energies (eV) at each distance.
        output_path: Path to save figure (e.g. path.png).
        title: Plot title.
        ylabel: Y-axis label.
        xlabel: X-axis label.
    """
    import matplotlib.pyplot as plt

    d = np.array(distances, dtype=float)
    e = np.array(energies, dtype=float)

    fig, ax = plt.subplots()
    ax.plot(d, e, "o-", color="C0", markersize=6)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=150)
    plt.close(fig)
