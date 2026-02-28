// =============================================================================
// Parallel MPI version of GAA Nanowire FET Simulator
// Schrodinger-Poisson Self-Consistent Solver
// Uses: ALUGrid (parallel partitioned), PETSc/SLEPc (MPI eigenvalue solver)
//       DUNE-ISTL with Overlapping Schwarz for Poisson
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
//#include <dune/common/mpihelper.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/parallel/indexset.hh>
#include <dune/common/parallel/communicator.hh>

// DUNE Grid
//#include <dune/grid/alugrid.hh>
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
#include "schrodingerfem_gaa_parallel.hh"
#include "poissonfem_NR_parallel.hh"
#include "shapefunctions.hh"
#include "postprocessing.hh"
#include "qdsimtimer.hh"

// ======================== Utility: Parallel-Safe I/O ==========================
// Only rank 0 reads files; data is broadcast to all ranks
// =============================================================================

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
      std::cout << "====================================================\n"
                << "  Parallel GAA Nanowire FET Simulator\n"
                << "  Running on " << size << " MPI processes\n"
                << "====================================================\n";

    MPI_Barrier(MPI_COMM_WORLD);
    std::cout << "  [Rank " << rank << "] initialized." << std::endl;
    MPI_Barrier(MPI_COMM_WORLD);

    // =================== Initialize Timer ==============================
    QDSimTimer timer;

    // ================== Mathematical Constants ========================
    const double pai = 3.142;

    // ============== Hartree Constants ==================================
    const double q0   = 1;
    const double m0   = 1;
    const double h_x  = 1;
    const double eps0 = 1 / (4 * pai);
    const double Hr   = 27.211;
    const double Hr_T = 3.157e5;
    const double au   = 5.291e-9;

    // ===================== Simulation Parameters ======================
    const double T = 300;

    // Material parameters - rank 0 reads, then broadcasts
    double Ns = 0, ni = 0, gv1 = 0, gv2 = 0, ml = 0, mt = 0;
    double Eg = 0, eps = 0, Egox1 = 0, meox1 = 0, mhox1 = 0, epsox1 = 0;
    double Egox2 = 0, meox2 = 0, mhox2 = 0, epsox2 = 0;
    double delEc1 = 0, delEv1 = 0, delEc2 = 0, delEv2 = 0;

    timer.start("FileIO_Read");

    if (rank == 0) {
      std::ifstream mat_para;
      mat_para.open("mat_para.dat");
      mat_para >> Ns >> ni >> Eg >> eps >> Egox1 >> meox1 >> mhox1 >> epsox1
               >> Egox2 >> meox2 >> mhox2 >> epsox2 >> delEc1 >> delEv1
               >> delEc2 >> delEv2;
      mat_para.close();
      std::cout << "  [Rank 0] Read material parameters." << std::endl;
    }

    // Broadcast material parameters to all ranks
    {
      double params[16] = {Ns, ni, Eg, eps, Egox1, meox1, mhox1, epsox1,
                           Egox2, meox2, mhox2, epsox2, delEc1, delEv1,
                           delEc2, delEv2};
      MPI_Bcast(params, 16, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      Ns = params[0]; ni = params[1]; Eg = params[2]; eps = params[3];
      Egox1 = params[4]; meox1 = params[5]; mhox1 = params[6]; epsox1 = params[7];
      Egox2 = params[8]; meox2 = params[9]; mhox2 = params[10]; epsox2 = params[11];
      delEc1 = params[12]; delEv1 = params[13]; delEc2 = params[14]; delEv2 = params[15];
    }

    // Unit conversions
    Ns    = Ns * au * au * au;
    ni    = ni * au * au * au;
    Eg    = Eg / Hr;
    Egox1 = Egox1 / Hr;
    Egox2 = Egox2 / Hr;

    const double solver = 1;
    const double k    = 8.61735E-5;
    const double beta = Hr / (k * T);

    typedef std::vector<float> vec;
    double Ef = 0.0;

    // --------------------- Reading grid file name --------------------------------
    // Mesh file is in the src directory
    // Mesh file - can be set via NANOWIRE_MESH env var
    std::string str = "/home/athira/dune-projects/dune-nanowire/src/gaa_25k.msh";
    const char* mesh_env = std::getenv("NANOWIRE_MESH");
    if (mesh_env) str = mesh_env;
    if (rank == 0)
      std::cout << "  Name of Grid file included = " << str << std::endl;

    // ======================= Grid Processing ============================
    // ALUGrid supports parallel partitioning natively
    // All ranks read the mesh; ALUGrid partitions automatically
    timer.stop("FileIO_Read");

    timer.start("Grid_Setup");

    const int dim = 2;
    const std::string gridfile = str;

    // Use ALUSimplexGrid with MPI-parallel support
    typedef Dune::ALUGrid<dim, dim, Dune::simplex, Dune::conforming> Grid;
    typedef std::unique_ptr<Grid> GridPointer;

    typedef Dune::GmshReader<Grid> GridReader;
    typedef Dune::GridFactory<Grid> Factory;

    // Grid-related types
    typedef Grid::LeafGridView GV;
    typedef Grid::LeafIntersectionIterator GridIntersectionIterator;

    typedef GV::Codim<0>::Iterator ElementIterator;
    typedef GV::IntersectionIterator IntersectionIterator;
    typedef GV::Codim<dim>::Iterator VertexIterator;

    typedef Dune::MultipleCodimMultipleGeomTypeMapper<GV> ElementMap;
    
    typedef Dune::MultipleCodimMultipleGeomTypeMapper<GV> VertexMap;
    typedef GV::IndexSet LeafIndexSet;

    typedef GV::ctype ctype;
    typedef Dune::BlockVector<Dune::FieldVector<ctype, 1> > Vector;

    

    // ==================== Create the Grid ==========================
    // All ranks participate in grid creation for parallel ALUGrid
    std::vector<int> boundaryid;
    std::vector<int> elementid;
    Factory factory;

    Dune::GmshReader<Grid>::read(factory, gridfile, boundaryid, elementid, false, false);
    GridPointer gridptr = factory.createGrid();
    Grid& grid = *gridptr;

    // Load balance: distribute elements across MPI ranks
    // ALUGrid handles this automatically when MPI is active
    timer.start("LoadBalance");
    gridptr->loadBalance();
    timer.stop("LoadBalance");

    // Define the grid view (each rank sees its local partition + ghost/overlap)
    const GV& gv_local = grid.leafGridView();

    timer.stop("Grid_Setup");

    if (rank == 0)
      Dune::gridinfo(grid, "GAA Nanowire Parallel Grid");

    const int N_local = gv_local.size(dim);
    const int N_elements_local = gv_local.size(0);
    const LeafIndexSet& set(gv_local.indexSet());

    // Gather global vertex count
    int N_global = 0;
    MPI_Allreduce(&N_local, &N_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    if (rank == 0)
      std::cout << "  Total vertices (sum across ranks): " << N_global
                << "  Local on rank 0: " << N_local << std::endl;

    // ========================= Data Mappers ==============================
    ElementMap elementmap(gv_local, Dune::mcmgElementLayout());
    VertexMap vertexmap(gv_local, Dune::mcmgVertexLayout());

    // ======================== Simulation Setup ============================
    VertexIterator vtend = gv_local.end<dim>();
    ElementIterator itend = gv_local.end<0>();

    Vector Nen(vertexmap.size());
    Vector Nep(vertexmap.size());
    Vector wavefunction_bc(vertexmap.size());
    Vector epsi(elementmap.size());
    Vector emass(elementmap.size());
    Vector hmass(elementmap.size());
    Vector em(elementmap.size());
    Vector Tr(vertexmap.size());

    HO<double, dim> parabola;

    // Region IDs from element centroid radius
    // After loadBalance(), factory.insertionIndex() is invalid.
    // Instead, assign region IDs geometrically based on centroid distance from origin.
    // This matches the GAA nanowire concentric-circle geometry:
    //   PhysID 201: channel (r <= R_nw)
    //   PhysID 200: oxide 1 / IL (R_nw < r <= R_nw + tox1)
    //   PhysID 199: oxide 2 / gate oxide (r > R_nw + tox1)
    double R_nw_geom, tox1_geom, tox2_geom;
    if (rank == 0) {
      std::ifstream dp("dev_para.dat");
      dp >> R_nw_geom >> tox1_geom >> tox2_geom;
      dp.close();
    }
    MPI_Bcast(&R_nw_geom, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&tox1_geom, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&tox2_geom, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    Vector regionid(elementmap.size());
    for (ElementIterator it = gv_local.begin<0>(); it != gv_local.end<0>(); ++it)
    {
      // Compute element centroid radius
      auto center = it->geometry().center();
      double r_cent = std::sqrt(center[0]*center[0] + center[1]*center[1]);

      int phys_id;
      if (r_cent <= R_nw_geom)
        phys_id = 201;  // channel
      else if (r_cent <= R_nw_geom + tox1_geom)
        phys_id = 200;  // oxide 1 (IL)
      else
        phys_id = 199;  // oxide 2 (gate oxide)

      regionid[elementmap.index(*it)] = phys_id;
    }

    Nen = 0.0;

    // =================== Initial Potential ================================
    // Initialize potential to zero (flat band). No file needed for initial run.
    // For self-consistent iteration, potential will be updated by Poisson solver.
    Vector potential(vertexmap.size());
    potential = 0.0;

    if (rank == 0)
      std::cout << "  Initialized flat-band potential (V=0)." << std::endl;

    // =================== Assign Material Properties =======================
    timer.start("MaterialSetup");
    for (ElementIterator it = gv_local.begin<0>(); it != gv_local.end<0>(); ++it)
    {
      int rid = (int)regionid[elementmap.index(*it)];

      if (rid == 199 || rid == 299) {  // oxide region 2
        epsi[elementmap.index(*it)]  = eps0 * epsox2;
        emass[elementmap.index(*it)] = meox2;
        hmass[elementmap.index(*it)] = mhox2;
      }

      if (rid == 200 || rid == 300) {  // oxide region 1 (IL)
        epsi[elementmap.index(*it)]  = eps0 * epsox1;
        emass[elementmap.index(*it)] = meox1;
        hmass[elementmap.index(*it)] = mhox1;
      }

      if (rid == 201 || rid == 301) {  // channel (substrate)
        epsi[elementmap.index(*it)] = eps0 * eps;
      }
    }

    // =================== Solve Schrodinger Equation ========================
    timer.stop("MaterialSetup");

    timer.start("Schrodinger_Total");
    if (rank == 0)
      std::cout << "\n  Solving Schrodinger equation for electrons (parallel)..." << std::endl;

    MPI_Barrier(MPI_COMM_WORLD);

    SchrodingerFEM_Parallel<GV, Vector> schelec(
      gv_local, Nen, potential, wavefunction_bc, 0, -Ef,
      emass, regionid, argc, argv);
    schelec.apply();

    MPI_Barrier(MPI_COMM_WORLD);

    if (rank == 0)
      std::cout << "  Schrodinger solve completed." << std::endl;

    timer.stop("Schrodinger_Total");

    // =================== Write Coordinates (rank 0 only) ====================
    timer.start("FileIO_Write");
    if (rank == 0) {
      std::ofstream plots;
      plots.open("coordinates.dat");
      for (VertexIterator vt = gv_local.begin<dim>(); vt != vtend; ++vt) {
        plots << vt->geometry().corner(0)[0] << "        "
              << vt->geometry().corner(0)[1] << "        "
              << vertexmap.index(*vt) << std::endl;
      }
      plots.close();
    }

    // =================== Store Potential (parallel VTK) =====================
    // Use DUNE's parallel VTK writer
    {
      Dune::VTKWriter<GV> vtkwriter(gv_local, Dune::VTK::conforming);
      vtkwriter.addVertexData(potential, "potential");

      // pwrite creates rank-specific files + .pvtu master file
      vtkwriter.pwrite("potential_output", ".", "",
                       Dune::VTK::appendedraw);
    }

    // Rank 0 writes the serial-style output for backward compatibility
    if (rank == 0) {
      std::ofstream plots;
      plots.open("plots1.dat");
      for (VertexIterator vt = gv_local.begin<dim>(); vt != vtend; ++vt) {
        plots << (vt->geometry().corner(0)[0]) << ","
              << (vt->geometry().corner(0)[1]) << ","
              << 0.0 << ","  // z=0 for 2D
              << potential[vertexmap.index(*vt)] << std::endl;
      }
      plots.close();
    }

    MPI_Barrier(MPI_COMM_WORLD);
    timer.stop("FileIO_Write");

    // =================== Print Timing Summary ==============================
    timer.printSummary();
    timer.writeScalingData("scaling_data.csv");

    if (rank == 0)
      std::cout << "\n  Simulation complete. Output written." << std::endl;

    return 0;

  }
  catch (Dune::Exception& e) {
    std::cerr << "Dune reported error: " << e << std::endl;
  }
  catch (...) {
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}
