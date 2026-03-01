
# QuDSim — Open-Source 3D Finite Element Simulation Tool for Nanoscale Devices

**QuDSim** is an open-source, MPI-parallel device simulation framework based on the finite element method (FEM) for modeling the electrical characteristics of nanoscale semiconductor devices. It provides self-consistent Schrödinger–Poisson solvers, polynomial eigenvalue–based tunneling lifetime extraction, and a Dirac equation solver for graphene quantum dots — all within a single extensible architecture.

> **Reference:** A. S, R. Nimje, R. Solanki, A. K. Jain, R. K. Pandey, and A. Mahajan, "QuDSim: An Open-Source 3D Finite Element Simulation Tool for Nanoscale Devices," *IEEE Transactions on Computer-Aided Design of Integrated Circuits and Systems* (under review), 2025.

---

## Table of Contents

- [Key Features](#key-features)
- [Architecture](#architecture)
- [Repository Structure](#repository-structure)
- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Running Simulations](#running-simulations)
- [Application Examples](#application-examples)
  - [Self-Consistent Schrödinger–Poisson Solver](#1-self-consistent-schrödinger-poisson-solver)
  - [Nanocrystal Flash Memory](#2-nanocrystal-flash-memory)
  - [Gate-All-Around Nanowire FET](#3-gate-all-around-nanowire-fet)
  - [Graphene Quantum Dot](#4-graphene-quantum-dot)
- [MPI Parallelization & Performance](#mpi-parallelization--performance)
- [Graphical User Interface](#graphical-user-interface)
- [Citation](#citation)
- [Contributing](#contributing)
- [Contact](#contact)
- [Acknowledgments](#acknowledgments)

---

## Key Features

- **Self-consistent Schrödinger–Poisson solver** with multi-valley effective mass support for obtaining electrostatic potential, quantized subbands, and charge density distributions.
- **Polynomial eigenvalue problem (PEP) formulation** for extracting tunneling lifetimes from quasi-bound states — avoids WKB approximations used by commercial TCAD tools.
- **Quantum tunneling current models**: direct tunneling, Fowler–Nordheim tunneling, trap-assisted tunneling, and quantum well-to-well tunneling via Bardeen's transfer Hamiltonian.
- **Dirac equation solver** for graphene quantum dots, formulated as a polynomial eigenvalue problem.
- **MPI-parallel execution** using domain decomposition with demonstrated 9.42× speedup on 32 cores.
- **Arbitrary geometry support** via unstructured FEM meshes (GMSH), including 2D, 3D, and symmetry-reduced configurations.
- **Lightweight GUI** (Python/Tkinter) for mesh editing, simulation configuration, and result visualization.
- **Extensible design**: users can incorporate custom governing equations, material models, and boundary conditions.

---

## Architecture

QuDSim is implemented in C++ and integrates the following open-source libraries:

| Component | Library | Role |
|-----------|---------|------|
| Mesh & grid management | [DUNE](https://www.dune-project.org/) + ALUGrid | Unstructured FEM grid, parallel load balancing |
| Sparse linear solvers | [PETSc](https://petsc.org/) | Distributed matrices, KSP iterative solvers |
| Eigenvalue solvers | [SLEPc](https://slepc.upv.es/) | Standard, generalized, and polynomial eigenvalue problems |
| Mesh generation | [GMSH](https://gmsh.info/) | 2D/3D unstructured mesh generation |
| Dense linear algebra | LAPACK, SuperLU | LU factorization, direct solvers |
| Visualization | [ParaView](https://www.paraview.org/) | VTK-based post-processing and 3D visualization |

All device physics formulations, FEM operator assembly, open-boundary condition implementation, self-consistent algorithms, and MPI integration are **original code** developed by the authors. The external libraries are used as unmodified computational backends.

---


## Prerequisites

### System Requirements

- Linux (tested on Ubuntu 20.04/22.04, CentOS 7/8, Red Hat)
- C++17-compatible compiler (GCC ≥ 9 or Clang ≥ 10)
- CMake ≥ 3.16
- MPI implementation (OpenMPI or MPICH)

### Required Libraries

| Library | Minimum Version | Installation |
|---------|----------------|-------------|
| DUNE | 2.10 | [dune-project.org](https://www.dune-project.org/) |
| ALUGrid | 2.10 | Included with DUNE |
| PETSc | 3.18 | [petsc.org](https://petsc.org/) |
| SLEPc | 3.18 | [slepc.upv.es](https://slepc.upv.es/) |
| GMSH | 4.11 | `sudo apt install gmsh` or [gmsh.info](https://gmsh.info/) |
| LAPACK | — | `sudo apt install liblapack-dev` |
| SuperLU | 5.3 | `sudo apt install libsuperlu-dev` |

### Optional

- **ParaView** for visualization of `.vtu` output files
- **Python 3.8+** with Tkinter for the GUI

---

## Installation

### 1. Install system dependencies (Ubuntu/Debian)

```bash
sudo apt update
sudo apt install build-essential cmake g++ gfortran \
    libopenmpi-dev openmpi-bin \
    liblapack-dev libblas-dev libsuperlu-dev \
    gmsh python3-tk
```

### 2. Install PETSc

```bash
git clone -b release https://gitlab.com/petsc/petsc.git petsc
cd petsc
./configure --with-cc=mpicc --with-cxx=mpicxx --with-fc=mpif90 \
    --with-debugging=0 --with-shared-libraries=1 \
    --download-superlu_dist
make all
export PETSC_DIR=$(pwd)
export PETSC_ARCH=arch-linux-c-opt
```

### 3. Install SLEPc

```bash
git clone -b release https://gitlab.com/slepc/slepc.git slepc
cd slepc
export SLEPC_DIR=$(pwd)
./configure
make all
```

### 4. Install DUNE (with ALUGrid)

Follow the instructions at [dune-project.org/doc/installation](https://www.dune-project.org/doc/installation/) to install the DUNE core modules and ALUGrid.

### 5. Clone and build QuDSim

```bash
git clone https://github.com/qudsimcnrnano/QuDSimprojects.git
cd QuDSimprojects
```

Build a specific module using `dunecontrol`:

```bash
/path/to/dune-common/bin/dunecontrol \
    --opts=config321.opts \
    --only=dune-selfconsistentsolver all
```

Replace `dune-selfconsistentsolver` with the desired module name (`dune-nanowire`, `dune-graphene`, or `dune-poisson3d`).

---

## Running Simulations

### Self-consistent solver (MOSCAP example)

```bash
cd dune-selfconsistentsolver/build-cmake/src

# Serial execution
./dune-selfconsistentsolver input.geo

# Parallel execution with 8 MPI processes
mpirun -np 8 ./dune-selfconsistentsolver input.geo
```

### Multiple bias points in parallel

QuDSim supports distributing different bias points across MPI process groups:

```bash
# Run 5 gate voltages × 8 cores each = 40 total cores
mpirun -np 40 ./dune-selfconsistentsolver input.geo --bias-parallel 5
```

### Output

Simulation results are exported as `.vtu` files that can be visualized with ParaView:

```bash
paraview output/potential.vtu
```

---

## Application Examples

### 1. Self-Consistent Schrödinger–Poisson Solver

**Module:** `dune-selfconsistentsolver`

Solves the coupled Schrödinger and Poisson equations iteratively until convergence, computing quantized subbands, charge density, and electrostatic potential for MOS structures. Supports multi-valley effective mass models for silicon (separate longitudinal and transverse masses) and isotropic mass for III–V semiconductors.

### 2. Nanocrystal Flash Memory

**Module:** `dune-selfconsistentsolver` + `dune-nanowire`

Models program, erase, and retention operations in nanocrystal-based flash memory with high-κ interpoly dielectrics (HfO₂, Al₂O₃, HfAlO). Features include:
- Quasi-bound state lifetime extraction via polynomial eigenvalue formulation
- Bardeen's transfer Hamiltonian for charge capture/emission
- Single-trap and multi-trap percolation decay models
- Validated against experimental data from Molas et al. (IEDM 2007)

### 3. Gate-All-Around Nanowire FET

**Module:** `dune-nanowire`

Computes gate leakage currents in cylindrical GAA nanowire MOSFETs by solving the open-boundary Schrödinger equation as a quadratic eigenvalue problem. Features include:
- Radial tunneling lifetime extraction using Hankel function boundary conditions
- Multi-material support: Si/SiO₂, InAs/Al₂O₃, In₀.₅₃Ga₀.₄₇As/Al₂O₃, GaN/Al₂O₃
- Composite oxide stack analysis (Al₂O₃/HfO₂)
- Strain-dependent leakage via deformation potential theory
- Validated against experimental data from Spinelli et al. (IEEE TED, 2012)

### 4. Graphene Quantum Dot

**Module:** `dune-graphene`

Solves the Dirac-like equation for massless fermions in graphene quantum dots using a weak formulation that decouples into two Schrödinger-like polynomial eigenvalue problems. Computes ground and excited state wavefunctions under vanishing boundary conditions.

---

## MPI Parallelization & Performance

QuDSim uses a two-level parallelization strategy:

1. **Domain decomposition** via DUNE-ALUGrid's space-filling curve (SFC) partitioning
2. **Task-level parallelism** across bias points



## Graphical User Interface

The GUI is built with Python's Tkinter library and provides:

- **Geometry editing**: load and modify `.geo` files for GMSH
- **Mesh visualization**: view `.msh` files with quality indicators
- **Simulation control**: configure solver parameters and launch runs
- **Result visualization**: process `.vtu` files via ParaView integration

```bash
cd GUI
python3 main.py
```

---

## Citation

If you use QuDSim in your research, please cite:

```bibtex
@article{qudsim2025,
  author  = {Athira, S. and Nimje, Rohit and Solanki, Ravi and Jain, Aakash Kumar 
             and Pandey, Rajan Kumar and Mahajan, Ashutosh},
  title   = {{QuDSim}: An Open-Source {3D} Finite Element Simulation Tool 
             for Nanoscale Devices},
  journal = {IEEE Transactions on Computer-Aided Design of Integrated Circuits 
             and Systems},
  year    = {2025},
  note    = {Under review}
}
```

---

## Contributing

We welcome contributions from the community. To contribute:

1. Fork the repository
2. Create a feature branch: `git checkout -b feature/your-feature`
3. Commit your changes: `git commit -m 'Add new feature'`
4. Push to the branch: `git push origin feature/your-feature`
5. Open a Pull Request

For bug reports and feature requests, please use the [Issues](https://github.com/qudsimcnrnano/QuDSimprojects/issues) tab.

---

## Contact

**Maintainer:** Dr. Ashutosh Mahajan

Centre for Nanotechnology Research (CNR)  
Vellore Institute of Technology (VIT), Vellore – 632014, India

- **Email:** ashutosh.mahajan@vit.ac.in
- **GitHub:** [qudsimcnrnano](https://github.com/qudsimcnrnano)

---

## Acknowledgments

This work was supported by the **Science and Engineering Research Board (SERB)**, Government of India, under Grant No. CRG/2022/009394.



---

## License

This project is released as open-source software for academic and research use. See the repository for license details.
