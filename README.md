# QuDSim: 3D Finite Element Simulation Tool for Nanoscale Devices

## Overview

QuDSim is an open-source, computationally efficient 3D simulator designed for simulating the electrical characteristics of nanoscale semiconductor devices, particularly at advanced technology nodes. Built upon the Finite Element Method (FEM), QuDSim solves coupled differential equations across complex device geometries. It is capable of simulating both 2D and 3D device structures, including 3D-to-2D equivalent configurations due to geometric symmetries.

QuDSim integrates a range of advanced techniques such as solving the Schrodinger-Poisson equations, simulating tunneling currents, and evaluating quasi-bound states. The tool also supports MPI-based parallelization for large-scale simulations and features a user-friendly graphical interface.

## Key Features

- **Self-consistent simulation** of Schrodinger-Poisson equations for electrostatic potential and charge density.
- **Polynomial eigenvalue solvers** to compute tunneling current and lifetimes of quasi-bound states.
- **Graphene Quantum Dot Solver** using Dirac-like equations for graphene.
- **Gate leakage analysis** and **quantum tunneling simulations** for nanowire FETs.
- **3D-to-2D equivalence** to simulate complex geometries with reduced computational cost.
- **MPI Parallelization** for scalable simulations on multi-core architectures, achieving significant speedup (e.g., 9.42x on 32 cores).
- **Graphical User Interface (GUI)** for simplifying simulation workflows, including mesh generation, visualization, and result processing.

## Computational Framework

QuDSim is implemented in C++ and integrates several modular open-source software packages for enhanced flexibility and performance:
- **DUNE** (Distributed and Unified Numerics Environment) for mesh generation and grid-based solutions.
- **PETSc** for parallel sparse matrix operations and iterative solvers.
- **SLEPc** for scalable eigenvalue problem solutions.
- **GMSH** for mesh generation.
- **LAPACK**, **SuperLU** for linear algebra operations.
- **ParaView** for result visualization.

## Performance and Scalability

QuDSim has been optimized for parallel computing, and its performance is demonstrated through several benchmarks:
- **Speedup**: Achieved up to 9.42x speedup on 32 cores for self-consistent simulations.
- **Efficiency**: The parallel efficiency for the Schrodinger solver is consistently close to the theoretical maximum as per Amdahl’s Law.

## Examples

### 1. **NC Flash Memory**
Simulated the program, erase, and retention operations in floating-gate (FG) flash memory devices using quantum tunneling models. The framework accurately resolves tunneling currents and quantum confinement effects within the nanocrystal floating gate.

### 2. **Gate-All-Around Nanowire FET**
Computed leakage currents in gate-all-around (GAA) FETs using a multi-dimensional FEM approach. Simulated the tunneling current across oxide layers with varying thickness and demonstrated the impact of material properties on leakage.

### 3. **Graphene Quantum Dot**
Implemented a Dirac equation solver for graphene quantum dots, calculating eigenstates and tunneling current. QuDSim integrates a polynomial eigenvalue solver to obtain tunneling lifetimes for quasi-bound states.

## Getting Started

### Prerequisites

- C++ compiler
- CMake
- DUNE, PETSc, SLEPc, LAPACK, SuperLU libraries
- MPI support for parallel execution
- GMSH for mesh generation
- Python (for GUI and visualization support)

### Installation

Clone the repository:
```bash
git clone https://github.com/qudsimcnrnano/QuDSimprojects.git
cd QuDSimprojects

## For Ubuntu/Linux
sudo apt install build-essential cmake g++ libopenmpi-dev

mkdir build
cd build
cmake ..
make


### running
/opt/local/dune_210/dune-common-2.10.dev20221009/bin/dunecontrol --opts=config321.opts --only=dune-projectname all
