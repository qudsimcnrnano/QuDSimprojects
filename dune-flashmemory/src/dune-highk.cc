// =============================================================================
// Parallel MPI version of Flash Memory High-K Simulator
// Schrodinger-Poisson Self-Consistent Solver
// Uses: ALUGrid (parallel partitioned), PETSc/SLEPc (MPI eigenvalue solver)
//       DUNE-ISTL for Poisson
//
// Original serial code: dune_highk.cc
//
// Device structure (x-direction, gate to substrate):
//   x = 0          : Gate contact (Dirichlet BC)
//   x = [0, 5]     : HTO near gate     (region 503, eps=4)
//   x = [5, 12]    : HfO2 trapping     (region 502, eps=20)
//   x = [12, 15]   : HTO near NC       (region 501, eps=4)
//   x = [15, 19]   : SiO2 tunnel oxide (region 400, eps=3.9)
//   x = [19, 24]   : Si transition     (region 200, eps=11.7)
//   x = [24, 39]   : Si channel        (region 201, eps=11.7)
//   x = [39, 50]   : Si bulk           (region 202, eps=11.7)
//   x = 50         : Substrate contact (Dirichlet BC)
//
// All coordinates in nanometers. Mesh: grids/qmmos_highk.msh
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
    // Read from flash_mat_para.dat (rank 0 reads, broadcast to all)

    double Ns_raw = 2.5e15, ni_raw = 1.5e10, Eg_raw = 1.12, eps1 = 11.7;
    double Vfb_raw = -0.88, T = 300;
    double epsox = 3.9, Egox_raw = 9.0, meox = 0.5, mhox = 0.5;
    double delEc_sio2 = 3.16, delEv_sio2 = 4.5;
    double epshto = 4.0, Eghto_raw = 17.0, mehto = 0.4, mhhto = 0.4;
    double delEc_hto = 2.8, delEv_hto = 5.0;
    double epshfo2 = 20.0, Eghfo2_raw = 5.68, mehfo2 = 0.2, mhhfo2 = 0.15;
    double delEc_hfo2 = 1.5, delEv_hfo2 = 3.4;
    double Vss_raw = 14.0;
    double t_hto_gate = 5.0, t_hfo2 = 7.0, t_hto_nc = 3.0, t_sio2 = 4.0;
    double t_si_trans = 5.0, t_si_chan = 15.0, t_si_bulk = 11.0;
    double y_width = 50.0;

    if (rank == 0) {
      std::ifstream fpara("flash_mat_para.dat");
      if (fpara.is_open()) {
        std::string key;
        while (fpara >> key) {
          if (key[0] == '#') { std::getline(fpara, key); continue; }
          if      (key == "Ns")          fpara >> Ns_raw;
          else if (key == "ni")          fpara >> ni_raw;
          else if (key == "Eg")          fpara >> Eg_raw;
          else if (key == "eps_si")      fpara >> eps1;
          else if (key == "Vfb")         fpara >> Vfb_raw;
          else if (key == "T")           fpara >> T;
          else if (key == "epsox")       fpara >> epsox;
          else if (key == "Egox")        fpara >> Egox_raw;
          else if (key == "meox")        fpara >> meox;
          else if (key == "mhox")        fpara >> mhox;
          else if (key == "delEc_sio2")  fpara >> delEc_sio2;
          else if (key == "delEv_sio2")  fpara >> delEv_sio2;
          else if (key == "epshto")      fpara >> epshto;
          else if (key == "Eghto")       fpara >> Eghto_raw;
          else if (key == "mehto")       fpara >> mehto;
          else if (key == "mhhto")       fpara >> mhhto;
          else if (key == "delEc_hto")   fpara >> delEc_hto;
          else if (key == "delEv_hto")   fpara >> delEv_hto;
          else if (key == "epshfo2")     fpara >> epshfo2;
          else if (key == "Eghfo2")      fpara >> Eghfo2_raw;
          else if (key == "mehfo2")      fpara >> mehfo2;
          else if (key == "mhhfo2")      fpara >> mhhfo2;
          else if (key == "delEc_hfo2")  fpara >> delEc_hfo2;
          else if (key == "delEv_hfo2")  fpara >> delEv_hfo2;
          else if (key == "Vgs")         fpara >> Vss_raw;
          else if (key == "t_hto_gate")  fpara >> t_hto_gate;
          else if (key == "t_hfo2")      fpara >> t_hfo2;
          else if (key == "t_hto_nc")    fpara >> t_hto_nc;
          else if (key == "t_sio2")      fpara >> t_sio2;
          else if (key == "t_si_trans")  fpara >> t_si_trans;
          else if (key == "t_si_chan")   fpara >> t_si_chan;
          else if (key == "t_si_bulk")   fpara >> t_si_bulk;
          else if (key == "y_width")     fpara >> y_width;
          else { std::getline(fpara, key); }  // skip unknown keys
        }
        fpara.close();
        cout << "  [Rank 0] Read flash_mat_para.dat" << endl;
      } else {
        cout << "  WARNING: flash_mat_para.dat not found, using defaults" << endl;
      }
    }

    // Broadcast all parameters (pack into array)
    {
      double p[30] = { Ns_raw, ni_raw, Eg_raw, eps1, Vfb_raw, T,
                        epsox, Egox_raw, meox, mhox, delEc_sio2, delEv_sio2,
                        epshto, Eghto_raw, mehto, mhhto, delEc_hto, delEv_hto,
                        epshfo2, Eghfo2_raw, mehfo2, mhhfo2, delEc_hfo2, delEv_hfo2,
                        Vss_raw, t_hto_gate, t_hfo2, t_hto_nc, t_sio2, t_si_trans };
      MPI_Bcast(p, 30, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      Ns_raw=p[0]; ni_raw=p[1]; Eg_raw=p[2]; eps1=p[3]; Vfb_raw=p[4]; T=p[5];
      epsox=p[6]; Egox_raw=p[7]; meox=p[8]; mhox=p[9]; delEc_sio2=p[10]; delEv_sio2=p[11];
      epshto=p[12]; Eghto_raw=p[13]; mehto=p[14]; mhhto=p[15]; delEc_hto=p[16]; delEv_hto=p[17];
      epshfo2=p[18]; Eghfo2_raw=p[19]; mehfo2=p[20]; mhhfo2=p[21]; delEc_hfo2=p[22]; delEv_hfo2=p[23];
      Vss_raw=p[24]; t_hto_gate=p[25]; t_hfo2=p[26]; t_hto_nc=p[27]; t_sio2=p[28]; t_si_trans=p[29];

      double p2[4] = { t_si_chan, t_si_bulk, y_width, 0.0 };
      MPI_Bcast(p2, 4, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      t_si_chan=p2[0]; t_si_bulk=p2[1]; y_width=p2[2];
    }

    // Unit conversions (to Hartree atomic units)
    const double Ns  = Ns_raw * au * au * au;
    const double ni  = ni_raw * au * au * au;
    const double Eg  = Eg_raw / Hr;
    const double Vfb = Vfb_raw / Hr;

    // ============== Other Constants ===================================
    const double k    = 8.61735E-5;
    const double beta = Hr / (k * T);

    // ============================ Some Calculations ===================
    double Np0 = 0.5 * (Ns + sqrt(Ns * Ns + 4 * ni * ni));
    double Ef  = (1 / beta) * log(Np0 / ni);
    if (rank == 0) {
      cout << "  actual Ef = " << Ef << endl;
      cout << "  Parameters: Ns=" << Ns_raw << " ni=" << ni_raw << " Eg=" << Eg_raw
           << " eps_Si=" << eps1 << endl;
      cout << "  SiO2:  eps=" << epsox << " me*=" << meox << " delEc=" << delEc_sio2 << endl;
      cout << "  HTO:   eps=" << epshto << " me*=" << mehto << " delEc=" << delEc_hto << endl;
      cout << "  HfO2:  eps=" << epshfo2 << " me*=" << mehfo2 << " delEc=" << delEc_hfo2 << endl;
      cout << "  Vgs=" << Vss_raw << " V" << endl;
    }

    // ======================= Device geometry constants ================
    // Computed from layer thicknesses (read from file or defaults)
    const double x_gate     = 0.0;
    const double x_hto1_end = t_hto_gate;
    const double x_hfo2_end = t_hto_gate + t_hfo2;
    const double x_hto2_end = t_hto_gate + t_hfo2 + t_hto_nc;
    const double x_sio2_end = t_hto_gate + t_hfo2 + t_hto_nc + t_sio2;
    const double x_sub1_end = x_sio2_end + t_si_trans;
    const double x_sub2_end = x_sub1_end + t_si_chan;
    const double x_sub      = x_sub2_end + t_si_bulk;

    if (rank == 0) {
      cout << "  Device geometry (nm):" << endl;
      cout << "    HTO(gate) [503]: x=[0, " << x_hto1_end << "]" << endl;
      cout << "    HfO2     [502]: x=[" << x_hto1_end << ", " << x_hfo2_end << "]" << endl;
      cout << "    HTO(NC)  [501]: x=[" << x_hfo2_end << ", " << x_hto2_end << "]" << endl;
      cout << "    SiO2     [400]: x=[" << x_hto2_end << ", " << x_sio2_end << "]" << endl;
      cout << "    Si(trans)[200]: x=[" << x_sio2_end << ", " << x_sub1_end << "]" << endl;
      cout << "    Si(chan) [201]: x=[" << x_sub1_end << ", " << x_sub2_end << "]" << endl;
      cout << "    Si(bulk) [202]: x=[" << x_sub2_end << ", " << x_sub << "]" << endl;
    }

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
    // After loadBalance(), factory.insertionIndex() may be invalid.
    // Use geometric assignment based on element centroid x-coordinate.
    // This is robust for the structured flash memory mesh.
    Vector regionid(elementmap.size());

    for (ElementIterator it = gv.begin<0>(); it != gv.end<0>(); ++it)
    {
      auto center = it->geometry().center();
      double x_cent = center[0];

      int phys_id;
      if (x_cent < x_hto1_end)
        phys_id = 503;   // HTO near gate
      else if (x_cent < x_hfo2_end)
        phys_id = 502;   // HfO2 trapping
      else if (x_cent < x_hto2_end)
        phys_id = 501;   // HTO near NC
      else if (x_cent < x_sio2_end)
        phys_id = 400;   // SiO2 tunnel oxide
      else if (x_cent < x_sub1_end)
        phys_id = 200;   // Si transition
      else if (x_cent < x_sub2_end)
        phys_id = 201;   // Si channel
      else
        phys_id = 202;   // Si bulk

      regionid[elementmap.index(*it)] = phys_id;
    }

    if (rank == 0) {
      // Count elements per region for verification
      int cnt[8] = {0};
      int rids[7] = {503, 502, 501, 400, 200, 201, 202};
      for (ElementIterator it = gv.begin<0>(); it != gv.end<0>(); ++it) {
        int rid = (int)regionid[elementmap.index(*it)];
        for (int k = 0; k < 7; k++) if (rid == rids[k]) cnt[k]++;
      }
      cout << "  Region element counts (rank 0 partition):" << endl;
      const char* names[] = {"HTO(gate)", "HfO2", "HTO(NC)", "SiO2", "Si(trans)", "Si(chan)", "Si(bulk)"};
      for (int k = 0; k < 7; k++)
        cout << "    " << rids[k] << " [" << names[k] << "]: " << cnt[k] << " elements" << endl;
    }

    // ======================== Material Properties =========================
    timer.start("MaterialSetup");

    for (ElementIterator it = gv.begin<0>(); it != gv.end<0>(); ++it)
    {
      int rid = (int)regionid[elementmap.index(*it)];

      if (rid == 400) {   // SiO2 tunnel oxide
        epsi[elementmap.index(*it)]  = eps0 * epsox;
        emass[elementmap.index(*it)] = meox;
        hmass[elementmap.index(*it)] = mhox;
      }
      if (rid == 200 || rid == 201 || rid == 202) {   // Silicon substrate
        epsi[elementmap.index(*it)] = eps0 * eps1;
      }
      if (rid == 501) {   // HTO near nanocrystal
        emass[elementmap.index(*it)] = mehto;
        hmass[elementmap.index(*it)] = mhhto;
        epsi[elementmap.index(*it)]  = eps0 * epshto;
      }
      if (rid == 502) {   // HfO2 charge trapping layer
        emass[elementmap.index(*it)] = mehfo2;
        hmass[elementmap.index(*it)] = mhhfo2;
        epsi[elementmap.index(*it)]  = eps0 * epshfo2;
      }
      if (rid == 503) {   // HTO near gate
        emass[elementmap.index(*it)] = mehto;
        hmass[elementmap.index(*it)] = mhhto;
        epsi[elementmap.index(*it)]  = eps0 * epshto;
      }
    }

    timer.stop("MaterialSetup");

    // ===================== Gate Voltage ================================
    double Vss = Vss_raw / Hr;   // Gate voltage in Hartree units

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
        if (rid == 400)  barrier = delEv_sio2 / Hr;   // SiO2
        if (rid == 501)  barrier = delEv_hto / Hr;     // HTO near NC
        if (rid == 502)  barrier = delEv_hfo2 / Hr;    // HfO2
        if (rid == 503)  barrier = delEv_hto / Hr;     // HTO near gate

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
        if (rid == 400)  barrier = delEc_sio2 / Hr;   // SiO2
        if (rid == 501)  barrier = delEc_hto / Hr;     // HTO near NC
        if (rid == 502)  barrier = delEc_hfo2 / Hr;    // HfO2
        if (rid == 503)  barrier = delEc_hto / Hr;     // HTO near gate

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

      // ---------- Write intermediate VTK (every 50 iterations) ----------
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

        // Count vertices in oxide+gate region (x <= x_sio2_end)
        int jmax = 0;
        for (ElementIterator it = gv.begin<0>(); it != gv.end<0>(); ++it) {
          Dune::GeometryType gt = it->type();
          const auto& ref = Dune::ReferenceElements<ctype,dim>::general(gt);
          int vertexsize = ref.size(dim);
          for (int t = 0; t < vertexsize; t++) {
            int indexi = set.subIndex(*it, t, dim);
            if (it->geometry().corner(t)[0] <= x_sio2_end)
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
            if (it->geometry().corner(t)[0] <= x_sio2_end) {
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
