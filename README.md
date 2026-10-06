# wave-propagation-2d
Vectorized finite difference solver for 2D acoustic wave propagation on a staggered Yee lattice, solving the first-order hyperbolic Friedrichs system.
# 2D Acoustic Wave Propagation Simulator

This repository contains a high-performance, fully vectorized finite difference solver for simulating acoustic wave propagation in a two-dimensional unbounded domain. It resolves the first-order hyperbolic Friedrichs system using a spatial and temporal staggered Yee-lattice arrangement.

## Project Structure

* `src/wave_solver.py`: The production-grade Python script containing the core simulation loop utilizing NumPy SIMD array slicing.
* `notebooks/colab_runtime.ipynb`: Interactive Jupyter/Google Colab notebook with built-in live plotting cells.
* `report/`: Source files (`.tex`) and exported simulation figures used for the project documentation.

## Physical and Numerical Specifications

The current implementation runs under a homogeneous Dirichlet boundary condition (\(p = 0\)) on the outer artificial boundary (\(\Sigma\)). The baseline simulation uses the following parameters:

* **Half-box size (\(L\)):** 1.0
* **Grid resolution (\(N\)):** 120 cells per side
* **Spatial step size (\(h\)):** 2.0 * L / N
* **CFL parameter (\(\nu\)):** 0.6 (satisfies the analytical stability limit \(\nu \leq 1/\sqrt{2}\))
* **Total runtime (\(T\)):** 5.0 seconds
* **Acoustic Source:** Off-center Ricker wavelet localized at \(x_s = (0.1, -0.05)\) with \(\sigma = 0.08\) and base frequency \(f_0 = 1.5\).

## Staggered Grid Allocation Logic

To ensure centered finite differences without introducing numerical artificial dissipation, variables are offset in memory as follows:
* Pressure (\(p\)): Centered at standard grid nodes \((i, j)\) at integer time levels \(t^n\).
* Horizontal Velocity (\(v_x\)): Placed on horizontal cell edges \((i + 1/2, j)\) at half-integer time levels \(t^{n+1/2}\).
* Vertical Velocity (\(v_y\)): Placed on vertical cell edges \((i, j + 1/2)\) at half-integer time levels \(t^{n+1/2}\).

## How to Run

1. Ensure dependencies are installed:
   ```bash
   pip install numpy matplotlib
   ```
2. Execute the primary solver:
   ```bash
   python src/wave_solver.py
   ```
