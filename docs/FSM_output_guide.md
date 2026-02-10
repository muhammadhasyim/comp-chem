# FSM (Freezing String Method) Output Guide

When you run FSM with ORCA (or any calculator), you get two kinds of output:

1. **ORCA (or other calculator) files** – many `.inp`, `.out`, `.gbw`, etc. These are temporary single-point runs for each node. ASE runs the calculator (e.g. in the current working directory or in subdirs) and you can ignore or archive them after the run.
2. **FSM result files** – the actual output of the program. These live in a **single output directory** under your reaction directory.

---

## Where FSM writes its results

Path pattern:

```
<reaction_dir>/fsm_interp_<interp>_method_<method>_maxls_<maxls>_maxiter_<maxiter>_nnodesmin_<nnodes>_<calculator>_<suffix>/
```

Example for your run:

```
orca_inputs/fsm_reaction/fsm_interp_ric_method_L-BFGS-B_maxls_2_maxiter_2_nnodesmin_6_orca_run1/
```

---

## FSM output files

| File | Description |
|------|-------------|
| **vfile_01.xyz**, **vfile_02.xyz**, … | One multi-frame XYZ per FSM iteration. Each frame = one node along the path. **Comment line (line 2) of each frame:** `s  E_rel` = arc length (Å) and relative energy (eV). |
| **ngrad.txt** | Written when FSM has finished (string stopped growing). Single number = total gradient (energy+force) evaluations. |

---

## How to read the result

- **vfile_XX.xyz**  
  - Frame 1: reactant (s ≈ 0, E_rel ≈ 0).  
  - Last frame: product (largest s, E_rel often ≈ 0 again).  
  - **The frame with the largest E_rel is the transition-state (TS) guess.**

- **Comment line format**  
  For each structure the second line is e.g. `0.41491 0.000` → arc length 0.41 Å, relative energy 0.000 eV. The TS is the structure whose second line has the **maximum** second number.

- **Console**  
  At the end of the run you should see something like:  
  `Gradient calls: <N>`  
  and each iteration logs:  
  `INFO:root:ITERATION: k DIST: d ENERGY: [e0 e1 e2 ...]`  
  (relative energies in eV along the string).

---

## Summary

- The “bunch of ORCA outputs” are **temporary** calculator runs; the **real** FSM output is in the **fsm_interp_...** directory.
- Look there for **vfile_*.xyz** (path + energies) and **ngrad.txt** (when finished).
- The **TS guess** is the geometry in the vfile frame whose comment line has the **highest relative energy (eV)**.
