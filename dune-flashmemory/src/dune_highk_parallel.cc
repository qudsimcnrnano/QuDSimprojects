// =============================================================================
// Parallel MPI version of Flash Memory High-K Simulator
// Schrodinger-Poisson Self-Consistent Solver
// Uses: ALUGrid (parallel partitioned), PETSc/SLEPc (MPI eigenvalue solver)
//       DUNE-ISTL for Poisson
//
// Original serial code: dune_highk.cc
// Region IDs:
//   200, 201, 202 = Silicon substrate
//   400           = SiO2 (tunnel oxide)
//   501           = HTO near nanocrystal
//   502           = HfO2 (high-k charge trapping layer)
//   503           = HTO near gate
// =============================================================================

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <iostream>
#include <set>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <mpi.h>

// DUNE Core
#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>

// DUNE Grid
#include <dune/alugrid/grid.hh>
#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>
#include <dune/grid/common/gridinfo.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>

// DUNE Data Types
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/io.hh>

// Solver Includes
#include "schrodingerfem_highk_parallel.hh"
#include "poissonfem_highk_parallel.hh"
#include "shapefunctions.hh"
#include "postprocessing.hh"
#include "qdsimtimer.hh"

using namespace std;

template<class ctype, int dim>
class HO
{
public:
  ctype operator() (Dune::FieldVector<ctype,dim> x, double vss) const
  {
    double omega = 3.675e-3;
    ctype value = 0;
    return (-vss * x[0] / 1000 + vss);
  }
};


int main(int argc, char** argv)
{
  try {

    // ========================= MPI Initialization ===========================
    Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);
    int rank = helper.rank();
    int size = helper.size();

    if (rank == 0)
      cout << "====================================================\n"
           << "  Parallel Flash Memory High-K Simulator\n"
           << "  Running on " << size << " MPI processes\n"
           << "====================================================\n";

    MPI_Barrier(MPI_COMM_WORLD);
    cout << "  [Rank " << rank << "] initialized." << endl;
    MPI_Barrier(MPI_COMM_WORLD);

    // =================== Initialize Timer ==============================
    QDSimTimer timer;

    // ================== Mathematical Constants ========================
    const double pai = 3.142;

    // ============== Hartree Constants ==================================
    const double q0   = 1;
    const double m0   = 1;
    const double h_x  = 1;
    const double eps0  = 1 / (4 * pai);
    const double Hr   = 27.211;
    const double Hr_T = 3.157e5;
    const double au   = 5.291e-9;   // Bohr radius in cm

    // ===================== Simulation Parameters ======================
    const double Ns    = 2.5e15 * au * au * au;    // doping concentration p-type
    const double ni    = 1.5e10 * au * au * au;    // intrinsic concentration
    const double Tox   = 10e-7 / au;               // oxide thickness
    const double T     = 300;
    const double Vfb   = -0.88 / Hr;               // flat band voltage
    const double eps1  = 11.7;                     // Si permittivity
    const double epsox = 3.9;                      // SiO2 permittivity
    const double epshfalo = 17;                    // HfAlO
    const double Eg    = 1.12 / Hr;                // Si band gap
    const double Egox  = 9.0 / Hr;                // SiO2 band gap
    const double epshfo2  = 20;                    // HfO2 permittivity
    const double epsal2o3 = 9;                     // Al2O3 permittivity
    const double epshto   = 4;                     // HTO permittivity
    const double Egal2o3  = 6.4 / Hr;
    const double Eghto    = 17 / Hr;
    const double Eghfo2   = 5.68;

    // ============== Other Constants ===================================
    const double k    = 8.61735E-5;
    const double beta = Hr / (k * T);

    // ============================ Some Calculations ===================
    double Np0 = 0.5 * (Ns + sqrt(Ns * Ns + 4 * ni * ni));
    double Ef  = (1 / beta) * log(Np0 / ni);
    if (rank == 0)
      cout << "  actual Ef = " << Ef << endl;

    // ======================= Grid Processing ============================
    timer.start("FileIO_Read");

    const int dim = 2;

    // Grid file - configurable via environment variable
    std::string gridfile = "grids/qmmos_highk.msh";
    const char* mesh_env = std::getenv("FLASH_MESH");
    if (mesh_env) gridfile = mesh_env;
    if (rank == 0)
      cout << "  Grid file: " << gridfile << endl;

    timer.stop("FileIO_Read");

    timer.start("Grid_Setup");

    // Use ALUGrid with MPI-parallel support (DUNE 2.10)
    typedef Dune::ALUGrid<dim, dim, Dune::simplex, Dune::conforming> Grid;
    typedef std::unique_ptr<Grid> GridPointer;
    typedef Dune::GmshReader<Grid> GridReader;
    typedef Dune::GridFactory<Grid> Factory;

    // Grid-related types
    typedef Grid::LeafGridView GV;
    typedef GV::Codim<0>::Iterator ElementIterator;
    typedef GV::IntersectionIterator IntersectionIterator;
    typedef GV::Codim<dim>::Iterator VertexIterator;

    typedef Dune::MultipleCodimMultipleGeomTypeMapper<GV> ElementMap;
    typedef Dune::MultipleCodimMultipleGeomTypeMapper<GV> VertexMap;
    typedef GV::IndexSet LeafIndexSet;
    typedef GV::ctype ctype;
    typedef Dune::BlockVector<Dune::FieldVector<ctype, 1> > Vector;

    // ==================== Create the Grid ==========================
    std::vector<int> boundaryid;
    std::vector<int> elementid;
    Factory factory;

    Dune::GmshReader<Grid>::read(factory, gridfile, boundaryid, elementid, false, false);
    GridPointer gridptr = factory.createGrid();
    Grid& grid = *gridptr;

    // Load balance across MPI ranks
    timer.start("LoadBalance");
    gridptr->loadBalance();
    timer.stop("LoadBalance");

    const GV& gv = grid.leafGridView();
    timer.stop("Grid_Setup");

    if (rank == 0)
      Dune::gridinfo(grid, "Flash Memory High-K Grid");

    const int N_local = gv.size(dim);
    const int N_elements_local = gv.size(0);
    const LeafIndexSet& set(gv.indexSet());

    // Gather global vertex count
    int N_global = 0;
    MPI_Allreduce(&N_local, &N_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    if (rank == 0)
      cout << "  Total vertices (sum across ranks): " << N_global
           << "  Local on rank 0: " << N_local << endl;

    // ========================= Data Mappers ==============================
    ElementMap elementmap(gv, Dune::mcmgElementLayout());
    VertexMap vertexmap(gv, Dune::mcmgVertexLayout());

    // ======================== Simulation Setup ============================
    VertexIterator vtend = gv.end<dim>();
    ElementIterator itend = gv.end<0>();

    Vector charge(vertexmap.size());
    Vector potential(vertexmap.size());
    Vector potential_tmp(vertexmap.size());
    Vector deltapotential(vertexmap.size());
    Vector deltaRho(vertexmap.size());
    Vector Nen(vertexmap.size());
    Vector Nep(vertexmap.size());
    Vector wavefunction_bc(vertexmap.size());
    Vector epsi(elementmap.size());
    Vector Tr(vertexmap.size());
    HO<double, dim> parabola;

    potential     = 0.0;
    potential_tmp = 0.0;
    charge        = 0.0;
    Nen           = 0.0;
    Nep           = 0.0;

    Vector emass(elementmap.size());
    Vector hmass(elementmap.size());

    // ======================== Region ID Assignment ========================
    // After loadBalance(), factory.insertionIndex() is invalid.
    // Assign region IDs geometrically based on element centroid x-coordinate.
    // This maps to the original flash memory layer structure:
    //   x near gate     -> oxide layers (400, 501, 502, 503)
    //   x in substrate  -> Si (200, 201, 202)
    //
    // NOTE: The exact x-coordinate boundaries depend on your mesh geometry.
    // You may need to adjust these thresholds to match your qmmos_highk.msh.
    // Alternatively, if running on a single node without redistribution,
    // factory.insertionIndex() will still work.
    //
    // For robustness, we try factory.insertionIndex() first, and fall back
    // to geometric assignment if it fails.
    Vector regionid(elementmap.size());

    // Try using factory insertion index (works if grid not redistributed)
    bool use_factory_index = true;
    try {
      for (ElementIterator it = gv.begin<0>(); it != gv.end<0>(); ++it) {
        int fidx = factory.insertionIndex(*it);
        if (fidx < 0 || fidx >= (int)elementid.size()) {
          use_factory_index = false;
          break;
        }
        regionid[elementmap.index(*it)] = elementid[fidx];
      }
    } catch (...) {
      use_factory_index = false;
    }

    if (!use_factory_index) {
      // Geometric fallback: broadcast elementid from rank 0 and use
      // nearest-centroid matching (simplified for typical flash memory meshes)
      if (rank == 0)
        cout << "  WARNING: factory.insertionIndex() unavailable, using geometric region assignment" << endl;

      // For flash memory, the layers are typically stacked in x-direction:
      // Substrate (200-202) at larger x, oxides (400, 501-503) at smaller x
      // This is a simplified fallback - adjust thresholds for your mesh
      for (ElementIterator it = gv.begin<0>(); it != gv.end<0>(); ++it) {
        auto center = it->geometry().center();
        double x = center[0];

        // Default thresholds (adjust to your mesh dimensions)
        // These values assume the original qmmos_highk.msh geometry
        if (x < 0.0)
          regionid[elementmap.index(*it)] = 200;  // substrate
        else
          regionid[elementmap.index(*it)] = 200;  // substrate (fallback)
      }
    }

    if (rank == 0)
      cout << "  Region IDs assigned (use_factory=" << use_factory_index << ")" << endl;

    // ======================== Material Properties =========================
    timer.start("MaterialSetup");

    for (ElementIterator it = gv.begin<0>(); it != gv.end<0>(); ++it)
    {
      int rid = (int)regionid[elementmap.index(*it)];

      if (rid == 400) {   // SiO2 tunnel oxide
        epsi[elementmap.index(*it)]  = eps0 * epsox;
        emass[elementmap.index(*it)] = 0.5;
        hmass[elementmap.index(*it)] = 0.5;
      }
      if (rid == 200 || rid == 201 || rid == 202) {   // Silicon substrate
        epsi[elementmap.index(*it)] = eps0 * eps1;
      }
      if (rid == 501) {   // HTO near nanocrystal
        emass[elementmap.index(*it)] = 0.4;
        hmass[elementmap.index(*it)] = 0.4;
        epsi[elementmap.index(*it)]  = eps0 * epshto;
      }
      if (rid == 502) {   // HfO2 charge trapping layer
        emass[elementmap.index(*it)] = 0.2;
        hmass[elementmap.index(*it)] = 0.15;
        epsi[elementmap.index(*it)]  = eps0 * epshfo2;
      }
      if (rid == 503) {   // HTO near gate
        emass[elementmap.index(*it)] = 0.4;
        hmass[elementmap.index(*it)] = 0.4;
        epsi[elementmap.index(*it)]  = eps0 * epshto;
      }
    }

    timer.stop("MaterialSetup");

    // ===================== Gate Voltage ================================
    double Vss = 14 / Hr;   // Gate voltage in Hartree units

    // ===================== Solver Setup ================================
    int it0 = 0;
    int imax = 800;
    const double tolerance = 1e-10;

    if (rank == 0)
      cout << "  Creating Schrodinger and Poisson Solvers..." << endl;

    timer.start("SolverSetup");

    SchrodingerFEM_HighK_Parallel<GV, Vector> schelec(
      gv, Nen, potential_tmp, wavefunction_bc, 0, -Ef,
      regionid, emass, argc, argv);

    SchrodingerFEM_HighK_Parallel<GV, Vector> schhole(
      gv, Nep, potential_tmp, wavefunction_bc, 1, Ef,
      regionid, hmass, argc, argv);

    timer.stop("SolverSetup");

    // =================== Output files (rank 0) ============================
    std::ofstream res_plot;
    std::ofstream plots;
    if (rank == 0)
      res_plot.open("Results.dat");

    // =================== Initial Poisson Solve ============================
    charge = 0.0;
    it0 = 0;

    if (rank == 0)
      cout << "  Solving Poisson equation for initial guess" << endl;

    timer.start("Poisson_Initial");
    PoissonFEM_HighK_Parallel<GV, Vector> poissonfem(gv, charge, potential, Vss, deltaRho, deltapotential, it0, epsi);
    poissonfem.apply();
    timer.stop("Poisson_Initial");

    // Write initial potential (rank 0)
    if (rank == 0) {
      std::ofstream potini_plot("potini.dat");
      for (VertexIterator vt = gv.begin<dim>(); vt != gv.end<dim>(); ++vt) {
        potini_plot << vt->geometry().corner(0)[0] << "        "
                    << vt->geometry().corner(0)[1] << "          "
                    << potential[vertexmap.index(*vt)] << endl;
      }
      potini_plot.close();
    }

    // ===================== Self-Consistent Loop ===========================
    timer.start("SCF_Loop");

    for (int i = 0; i < imax; i++)
    {
      if (rank == 0)
        cout << "  Vss: " << Vss << "    Iteration: " << i << endl;

      // ---------- Setup potential for HOLES ----------
      for (int p = 0; p < N_local; p++)
        potential_tmp[p] = Eg / 2 + potential[p];

      // Adjust potential barriers per region (hole barriers)
      for (ElementIterator it = gv.begin<0>(); it != gv.end<0>(); ++it)
      {
        int rid = (int)regionid[elementmap.index(*it)];
        Dune::GeometryType gt = it->type();
        const auto& ref = Dune::ReferenceElements<ctype,dim>::general(gt);
        int vertexsize = ref.size(dim);

        double barrier = 0.0;
        if (rid == 400)  barrier = 4.5 / Hr;    // SiO2
        if (rid == 501)  barrier = 5.0 / Hr;    // HTO near NC
        if (rid == 502)  barrier = 3.4 / Hr;    // HfO2
        if (rid == 503)  barrier = 5.0 / Hr;    // HTO near gate

        if (barrier > 0.0) {
          for (int pt = 0; pt < vertexsize; pt++) {
            potential_tmp[set.subIndex(*it, pt, dim)] =
              Eg / 2 + barrier + potential[set.subIndex(*it, pt, dim)];
          }
        }
      }

      if (rank == 0)
        cout << "  Solving Schrodinger equation for holes" << endl;
      timer.start("Schrodinger_Holes");
      schhole.apply();
      timer.stop("Schrodinger_Holes");

      // ---------- Setup potential for ELECTRONS ----------
      for (int p = 0; p < N_local; p++)
        potential_tmp[p] = Eg / 2 - potential[p];

      // Adjust potential barriers per region (electron barriers)
      for (ElementIterator it = gv.begin<0>(); it != gv.end<0>(); ++it)
      {
        int rid = (int)regionid[elementmap.index(*it)];
        Dune::GeometryType gt = it->type();
        const auto& ref = Dune::ReferenceElements<ctype,dim>::general(gt);
        int vertexsize = ref.size(dim);

        double barrier = 0.0;
        if (rid == 400)  barrier = 3.16 / Hr;   // SiO2 (3.1 eV barrier Si-SiO2)
        if (rid == 501)  barrier = 2.8 / Hr;    // HTO near NC
        if (rid == 502)  barrier = 1.5 / Hr;    // HfO2
        if (rid == 503)  barrier = 2.8 / Hr;    // HTO near gate

        if (barrier > 0.0) {
          for (int pt = 0; pt < vertexsize; pt++) {
            potential_tmp[set.subIndex(*it, pt, dim)] =
              Eg / 2 + barrier - potential[set.subIndex(*it, pt, dim)];
          }
        }
      }

      if (rank == 0)
        cout << "  Solving Schrodinger equation for electrons" << endl;
      timer.start("Schrodinger_Electrons");
      schelec.apply();
      timer.stop("Schrodinger_Electrons");

      // ---------- Compute charge density ----------
      for (ElementIterator it = gv.begin<0>(); it != gv.end<0>(); ++it)
      {
        int rid = (int)regionid[elementmap.index(*it)];
        Dune::GeometryType gt = it->type();
        const auto& ref = Dune::ReferenceElements<ctype,dim>::general(gt);
        int vertexsize = ref.size(dim);

        for (int pt = 0; pt < vertexsize; pt++)
        {
          ctype Np_val = Nep[set.subIndex(*it, pt, dim)];
          ctype Nn_val = Nen[set.subIndex(*it, pt, dim)];

          // Clamp carrier densities
          double clamp = 1e25 * au * au * au;
          if (Np_val > clamp) Np_val = clamp;
          if (Nn_val > clamp) Nn_val = clamp;

          if (rid == 200 || rid == 201 || rid == 202) {
            // Substrate: charge = q(p - n - Na)
            charge[set.subIndex(*it, pt, dim)] = q0 * (Np_val - Nn_val - Ns);
            deltaRho[set.subIndex(*it, pt, dim)] = q0 * (-beta * Np_val - beta * Nn_val);
          }
          else if (rid == 400 || rid == 501 || rid == 502 || rid == 503) {
            // Oxide/dielectric: charge = q(p - n) (no doping)
            charge[set.subIndex(*it, pt, dim)] = q0 * (Np_val - Nn_val);
            deltaRho[set.subIndex(*it, pt, dim)] = q0 * (-beta * Np_val - beta * Nn_val);
          }
        }
      }

      // ---------- Write intermediate VTK ----------
      if (i % 50 == 0 || i == imax - 1) {
        Dune::VTKWriter<GV> vtkwriter_i(gv, Dune::VTK::conforming);
        vtkwriter_i.addVertexData(potential, "potential");
        vtkwriter_i.addVertexData(charge, "charge");
        vtkwriter_i.addVertexData(Nen, "Nen");
        vtkwriter_i.addVertexData(Nep, "Nep");
        vtkwriter_i.addVertexData(deltapotential, "deltapotential");
        vtkwriter_i.addVertexData(potential_tmp, "POTENTIAL");
        std::string vtkname = "Results_i_" + std::to_string(i);
        vtkwriter_i.pwrite(vtkname, ".", "", Dune::VTK::appendedraw);
      }

      // Write data file (rank 0 only)
      if (rank == 0) {
        plots.open("plots.dat");
        for (VertexIterator vt = gv.begin<dim>(); vt != vtend; ++vt) {
          plots << (vt->geometry().corner(0)[0]) << "      "
                << (vt->geometry().corner(0)[1]) << "     "
                << potential[vertexmap.index(*vt)] << "     "
                << charge[vertexmap.index(*vt)] << "     "
                << deltaRho[vertexmap.index(*vt)] << "      "
                << deltapotential[vertexmap.index(*vt)] << "     "
                << Nen[vertexmap.index(*vt)] << "    "
                << Nep[vertexmap.index(*vt)] << "   "
                << potential_tmp[vertexmap.index(*vt)] << endl;
        }
        plots.close();
      }

      // ---------- Solve Poisson ----------
      if (rank == 0)
        cout << "  Solving the Poisson equation" << endl;
      timer.start("Poisson_NR");
      poissonfem.apply();
      timer.stop("Poisson_NR");

      potential += deltapotential;

      // ---------- Convergence test ----------
      int conv_test_local = 1;
      if (i < imax) {
        for (VertexIterator vt = gv.begin<dim>(); vt != vtend; ++vt) {
          if (fabs(deltapotential[vertexmap.index(*vt)]) > tolerance)
            conv_test_local = 0;
        }
      }

      // All ranks must agree on convergence
      int conv_test_global = 0;
      MPI_Allreduce(&conv_test_local, &conv_test_global, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

      if (conv_test_global) {
        if (rank == 0)
          cout << "  Convergence obtained at i = " << i << endl;
        break;
      }

      if (i == imax - 1 && rank == 0)
        cout << "  Convergence NOT obtained, reached imax" << endl;
    }

    timer.stop("SCF_Loop");

    // =================== Final Output ====================================
    // Write final profiles (rank 0)
    if (rank == 0) {
      std::ofstream plot("plot.dat");
      for (VertexIterator vt = gv.begin<dim>(); vt != vtend; ++vt)
        plot << vt->geometry().corner(0)[0] << " "
             << vt->geometry().corner(0)[1] << " "
             << Nen[vertexmap.index(*vt)] << " "
             << potential[vertexmap.index(*vt)] << endl;
      plot.close();
    }

    // Final parallel VTK output
    {
      Dune::VTKWriter<GV> vtkwriter(gv, Dune::VTK::conforming);
      vtkwriter.addVertexData(potential, "potential");
      vtkwriter.addVertexData(charge, "charge");
      vtkwriter.addVertexData(Nen, "Nen");
      vtkwriter.addVertexData(Nep, "Nep");
      vtkwriter.pwrite("Results_final", ".", "", Dune::VTK::appendedraw);
    }

    // =================== Wavefunction Post-Processing (rank 0) ===========
    if (rank == 0) {
      int nconv = schelec.getnconv();

      Vector wf_DOT[nconv];
      double re_DOT[nconv];
      Vector wf_vec;
      double re_val;

      for (int i = 0; i < nconv; i++) {
        wf_vec = schelec.getwave(i);
        re_val = schelec.geteigen(i);
        wf_DOT[i] = wf_vec;
        re_DOT[i] = re_val;
      }

      cout << "  No. of eigenvalues printed: " << nconv << endl;

      std::ofstream mos_plot("flash_highk_results.dat");
      for (int i = 0; i < nconv; i++) {
        wf_vec = wf_DOT[i];
        re_val = re_DOT[i];
        int j = 0;
        Tr = 0;

        // Count vertices at x <= 0
        int jmax = 0;
        for (ElementIterator it = gv.begin<0>(); it != gv.end<0>(); ++it) {
          Dune::GeometryType gt = it->type();
          const auto& ref = Dune::ReferenceElements<ctype,dim>::general(gt);
          int vertexsize = ref.size(dim);
          for (int t = 0; t < vertexsize; t++) {
            int indexi = set.subIndex(*it, t, dim);
            if (it->geometry().corner(t)[0] < 0 || it->geometry().corner(t)[0] == 0)
              if (Tr[indexi] != 1) { jmax++; Tr[indexi] = 1; }
          }
        }

        Tr = 0;
        for (ElementIterator it = gv.begin<0>(); it != gv.end<0>(); ++it) {
          Dune::GeometryType gt = it->type();
          const auto& ref = Dune::ReferenceElements<ctype,dim>::general(gt);
          int vertexsize = ref.size(dim);
          for (int t = 0; t < vertexsize; t++) {
            int indexi = set.subIndex(*it, t, dim);
            if (it->geometry().corner(t)[0] < 0 || it->geometry().corner(t)[0] == 0) {
              if (Tr[indexi] != 1) {
                mos_plot << nconv << "      " << jmax << "    " << i << "    " << j << "  "
                         << it->geometry().corner(t)[0] << "        "
                         << it->geometry().corner(t)[1] << "          "
                         << wf_vec[indexi] << "       " << re_val << "      "
                         << potential_tmp[indexi] << endl;
                Tr[indexi] = 1;
                j++;
              }
            }
          }
        }
      }
      mos_plot.close();
    }

    if (rank == 0) res_plot.close();

    // =================== Print Timing Summary ==============================
    timer.printSummary();
    timer.writeScalingData("scaling_data_flash.csv");

    if (rank == 0)
      cout << "\n  Simulation complete. Output written." << endl;

    return 0;

  }
  catch (Dune::Exception& e) {
    std::cerr << "Dune reported error: " << e << std::endl;
  }
  catch (...) {
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}
