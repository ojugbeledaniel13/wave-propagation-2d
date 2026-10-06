# wave-propagation-2d
Vectorized finite difference solver for 2D acoustic wave propagation on a staggered Yee lattice, solving the first-order hyperbolic Friedrichs system.

# 2D Acoustic Wave Propagation Simulator

This repository contains a high-performance, fully vectorized finite difference solver for simulating acoustic wave propagation in a two-dimensional unbounded domain. It resolves the first-order hyperbolic Friedrichs system using a spatial and temporal staggered Yee-lattice arrangement.

## Project Structure

```text
wave-propagation-2d/
├── src/
│   └── wave_solver.py       # Core simulation script utilizing NumPy array slicing
├── notebooks/
│   └── colab_runtime.ipynb  # Interactive Google Colab notebook with plotting cells
├── report/
│   ├── figures/             # Saved .png wave snapshots for report documentation
│   └── main.tex             # Project report source LaTeX file
├── README.md                # Repository documentation manual
└── .gitignore               # Standard Python ignore file to exclude cached data
```

## Physical and Numerical Specifications

The current implementation runs under a homogeneous Dirichlet boundary condition (p = 0) on the outer artificial boundary (Sigma). The baseline simulation uses the following parameters:

* **Half-box size (L):** 1.0
* **Grid resolution (N):** 120 cells per side
* **Spatial step size (h):** 2.0 * L / N
* **CFL parameter (nu):** 0.6 (satisfies the analytical stability limit nu <= 0.707)
* **Total runtime (T):** 5.0 seconds
* **Acoustic Source:** Off-center Ricker wavelet localized at position (0.1, -0.05) with sigma = 0.08 and base frequency f0 = 1.5.

## Staggered Grid Allocation Logic

To ensure centered finite differences without introducing numerical artificial dissipation, variables are offset in memory as follows:
* **Pressure (p):** Centered at standard grid nodes (i, j) at integer time levels t^n.
* **Horizontal Velocity (vx):** Placed on horizontal cell edges (i + 1/2, j) at half-integer time levels t^(n+1/2).
* **Vertical Velocity (vy):** Placed on vertical cell edges (i, j + 1/2) at half-integer time levels t^(n+1/2).

## How to Run

Ensure dependencies are installed:
```bash
pip install numpy matplotlib
```

Execute the primary solver:
```bash
python src/wave_solver.py
```
