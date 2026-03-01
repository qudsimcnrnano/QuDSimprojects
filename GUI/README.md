
# QuDSim GUI — Graphical User Interface for QuDSim 3D FEM Tool

A lightweight Python/Tkinter-based graphical user interface for the QuDSim open-source 3D finite element simulation tool for nanoscale devices.

## Overview

This GUI provides an integrated environment for managing QuDSim simulations, from geometry creation through mesh generation, solver execution, and results visualization. It bridges the following computational modules into a unified interface:

| Module   | Purpose                              |
|----------|--------------------------------------|
| **DUNE** | Finite element grid management       |
| **GMSH** | 2D/3D mesh generation                |
| **PETSc**| Parallel sparse linear solvers       |
| **SLEPc**| Eigenvalue and PEP solvers           |
| **LAPACK**| Dense linear algebra                |
| **SuperLU**| Sparse direct LU factorization     |

## Features

- **Geometry Editor**: Edit `.geo` files with syntax highlighting
- **Mesh Generation**: One-click mesh generation via GMSH
- **Solver Execution**: Run MPI-parallel QuDSim solvers with configurable process counts
- **Results Visualization**: Open `.vtu` results in ParaView
- **Project Templates**: Pre-built templates for MOSCAP, GAA nanowire, graphene quantum dot, and NC flash memory structures
- **Material Database**: Built-in material parameter viewer
- **Input Parameter Editor**: Edit solver configuration files
- **Session Persistence**: Automatic save/restore of project state

## Installation

### Prerequisites

- **Python 3.6+** with Tkinter (included in standard Python installations)
- **GMSH** (for mesh generation): https://gmsh.info
- **ParaView** (for visualization): https://www.paraview.org
- **MPI** (for parallel solver execution): OpenMPI or MPICH
- **QuDSim** solver executables (compiled from the main repository)

### Running the GUI

```bash
cd QuDSim_GUI/
python3 main.py
```

No additional Python packages are required — the GUI uses only the standard library (Tkinter).

## File Structure

```
QuDSim_GUI/
├── main.py          # Main application entry point
├── load_file.py     # File loading (.geo, .msh, .vtu)
├── save_file.py     # Project saving and export
├── new_file.py      # New project creation with templates
└── README.md        # This file
```

## Usage

### Creating a New Project

1. **File → New Project** (or Ctrl+N)
2. Enter a project name and select a device template
3. The geometry file opens automatically in the built-in editor
4. Edit the geometry as needed

### Simulation Workflow

1. **Generate Mesh**: Click the GMSH button or use Simulation → Generate Mesh
2. **Configure**: Set MPI processes via Tools → Configure MPI
3. **Run Solver**: Click Solve or use Simulation → Run Solver
4. **Visualize**: Click ParaView or use Simulation → View Results

### Supported File Types

| Extension | Description | Action |
|-----------|-------------|--------|
| `.geo`    | GMSH geometry definition | Edit in built-in text editor |
| `.msh`    | GMSH mesh file | View info, open in GMSH |
| `.vtu`    | VTK simulation results | Open in ParaView |
| `.qdsproj`| QuDSim project file | Save/load project state |

## Module Check

Click any module button (DUNE, GMSH, PETSc, SLEPc, LAPACK, SuperLU) to verify if it is installed on your system. The GUI checks for:
- Command-line tools (`gmsh`, `dunecontrol`) in PATH
- Environment variables (`PETSC_DIR`, `SLEPC_DIR`, etc.)

## Device Templates

The GUI includes geometry templates for the three examples presented in the QuDSim paper:

1. **MOSCAP** — MOS capacitor for Schrödinger–Poisson simulations
2. **GAA Nanowire** — Gate-all-around cross-section for gate leakage analysis
3. **Graphene Quantum Dot** — Circular dot for Dirac equation solver
4. **NC Flash Memory** — Nanocrystal flash memory unit cell

## Citation



## License

This GUI is released as part of the QuDSim open-source project.

## Contact

For questions or feedback:
- Repository: https://github.com/qudsimcnrnano/QuDSimprojects
- Email: ashutosh.mahajan@vit.ac.in
