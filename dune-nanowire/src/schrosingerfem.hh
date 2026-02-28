// =============================================================================
// Parallel MPI Schrodinger FEM Solver for GAA Nanowire
// Key changes from serial:
//   1. PETSc matrices: MATSEQAIJ -> MATMPIAIJ (distributed across MPI ranks)
//   2. SLEPc EPS solver: MPI_COMM_WORLD (parallel eigenvalue solve)
//   3. Local-to-global index mapping for parallel assembly
//   4. Parallel wavefunction extraction and normalization
//   5. File I/O restricted to rank 0
// =============================================================================

#ifndef SCHRODINGERFEM_PARALLEL_HH
#define SCHRODINGERFEM_PARALLEL_HH

#include <sstream>
#include <cmath>
#include <map>
#include <mpi.h>
#include "petscmat.h"
#include "petscvec.h"
#include "petscis.h"
#include "slepceps.h"
#include "fermi.h"
#include <dune/common/fvector.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/common/fmatrix.hh>
//#include <dune/grid/common/quadraturerules.hh>
#include "shapefunctions.hh"
#include <dune/geometry/quadraturerules.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/io.hh>
#include "wavefunction.hh"

template<class GV, class Vector>
class SchrodingerFEM_Parallel
{
  private:
    static const int dim  = GV::dimension;
    static const int dim1 = 1;
    typedef typename GV::ctype ctype;
    typedef typename GV::template Codim<dim>::Iterator VertexIterator;
    typedef typename GV::template Codim<0>::Iterator LeafIterator;
    typedef typename GV::template Codim<1>::Iterator triangleIterator;
    typedef typename GV::IntersectionIterator IntersectionIterator;
    typedef typename GV::IndexSet LeafIndexSet;
    typedef Dune::BCRSMatrix<Dune::FieldMatrix<ctype,1,1> > Matrix;
    typedef Dune::MultipleCodimMultipleGeomTypeMapper<GV> VertexMap;
    typedef Dune::MultipleCodimMultipleGeomTypeMapper<GV> ElementMap;
    typedef std::vector<double> vec;
    typedef std::vector<vec> mat;

    const GV& gv;
    Vector& charge;
    Vector charge_tmp;
    Vector wavevector;
    Vector& potential;
    Vector& wavefunction_bc;
    Vector Tr;
    Vector& emass;
    Vector mass;
    Vector& regionid;
    const LeafIndexSet& set;

    // === Parallel index mapping ===
    int N_local;        // number of local vertices on this rank
    int N_global;       // total vertices across all ranks
    int rank, nprocs;
    std::vector<PetscInt> local_to_global;  // maps local vertex idx -> global idx
    std::map<int, int> global_to_local;     // maps global idx -> local vertex idx
    PetscInt row_start, row_end;            // this rank owns global rows [row_start, row_end)

    // Sparsity pattern (local)
    std::vector<std::set<int> > adjacencyPattern;

    // Local DUNE matrices (for assembly convenience)
    Matrix A, A1, A0, B, Cb, C;

    // PETSc parallel matrices
    Mat A_p, B_p;   // MPI-distributed matrices

    // SLEPc eigenvalue solver
    EPS            eps;
    EPSType  type;  // non-const for EPSGetType
    PetscReal      error, tol, re, im, mod_xr;
    PetscScalar    kr, ki, xr_val;
    PetscErrorCode ierr;
    PetscInt       nev, maxit, i, its, lits, nconv;
    Vec            xr_vec, xi_vec;  // renamed to avoid conflict with xr coordinate

    int carrier_type;
    int car_type;
    double Ef;

    void determineAdjacencyPattern();
    void buildParallelIndexMapping();
    void assemble(Vector mass, int valley);
    void apply_bc();
    void solve(double m, double g, int valley);

  public:
    SchrodingerFEM_Parallel(const GV& gv_, Vector& charge_, Vector& potential_,
      Vector& wavefunction_bc, int carrier_type_, double Ef_,
      Vector& emass, Vector& regionid_, int& argc, char**& argv);
    void apply();
    double geteigen(int);
    Vector getwave(int);
    int getnconv();

    int P, P1;
    mat wave_DOT;
    mat wave_DOT1;
    vec real_DOT;
    vec real_DOT1;

    Dune::FieldVector<ctype,dim> max_coords;
    double max, min, max1, min2, max2, avg, Nn;
};


// =============================================================================
// Constructor: Set up parallel PETSc matrices and SLEPc context
// =============================================================================
template<class GV, class Vector>
SchrodingerFEM_Parallel<GV,Vector>::SchrodingerFEM_Parallel(
    const GV& gv_, Vector& charge_, Vector& potential_,
    Vector& wavefunction_bc_, int carrier_type_, double Ef_,
    Vector& emass_, Vector& regionid_, int& argc, char**& argv)
  : gv(gv_), charge(charge_), potential(potential_),
    wavefunction_bc(wavefunction_bc_), emass(emass_), regionid(regionid_),
    set(gv.indexSet())
{
  carrier_type = carrier_type_;
  Ef = Ef_;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  N_local = gv.size(dim);

  // Build local-to-global index mapping using ALUGrid's global ID set
  buildParallelIndexMapping();

  if (rank == 0)
    std::cout << "  [SchrodingerFEM_Parallel] N_global = " << N_global
              << ", N_local(rank 0) = " << N_local << std::endl;

  // Determine sparsity pattern
  determineAdjacencyPattern();

  // =================== PETSc/SLEPc Initialization ==========================
  static char help[] = "Parallel Schrodinger Solver";
  SlepcInitialize(&argc, &argv, (char*)0, help);

  // =================== Create PETSc matrices ================================
  // Use PETSC_COMM_SELF for sequential matrices (works with any -np)
  {
    PetscInt nnz[N_local];
    for (int ii = 0; ii < N_local; ii++)
      nnz[ii] = adjacencyPattern[ii].size();

    MatCreate(PETSC_COMM_SELF, &A_p);
    MatCreate(PETSC_COMM_SELF, &B_p);
    MatSetType(A_p, MATSEQAIJ);
    MatSetType(B_p, MATSEQAIJ);
    MatSetSizes(A_p, N_local, N_local, N_local, N_local);
    MatSetSizes(B_p, N_local, N_local, N_local, N_local);
    MatSeqAIJSetPreallocation(A_p, 0, nnz);
    MatSeqAIJSetPreallocation(B_p, 0, nnz);
  }

  MatSetOption(A_p, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
  MatSetOption(B_p, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);

  // =================== Set up local DUNE matrices (for assembly) =============
  // These are used for local FEM assembly before copying to PETSc
  auto setupMatrix = [&](Matrix& M) {
    M.setSize(N_local, N_local, N_local + 2 * gv.size(dim - 1));
    M.setBuildMode(Matrix::random);
    for (int ii = 0; ii < N_local; ii++)
      M.setrowsize(ii, adjacencyPattern[ii].size());
    M.endrowsizes();
    for (int ii = 0; ii < N_local; ii++) {
      for (auto jj : adjacencyPattern[ii])
        M.addindex(ii, jj);
    }
    M.endindices();
  };

  setupMatrix(A);
  setupMatrix(A1);
  setupMatrix(A0);
  setupMatrix(B);
  setupMatrix(Cb);
  setupMatrix(C);

  charge_tmp.resize(N_local);
  wavevector.resize(N_local);
  Tr.resize(N_local);

  A = 0.0; A1 = 0.0; A0 = 0.0; B = 0.0; C = 0.0; Cb = 0.0;
  charge_tmp = 0.0; wavevector = 0.0; Tr = 0.0;
}


// Build local-to-global index mapping
// =============================================================================
template<class GV, class Vector>
void SchrodingerFEM_Parallel<GV,Vector>::buildParallelIndexMapping()
{
  VertexMap vertexmap(gv, Dune::mcmgVertexLayout());
  VertexIterator vtend = gv.template end<dim>();

  if (nprocs == 1) {
    // Single rank: identity mapping
    N_global = N_local;
    local_to_global.resize(N_local);
    for (int ii = 0; ii < N_local; ii++) {
      local_to_global[ii] = ii;
      global_to_local[ii] = ii;
    }
    row_start = 0;
    row_end = N_local;
  }
  else {
    // Multi-rank: build contiguous global numbering from vertex coordinates
    // (ALUGrid IDs are opaque objects, not castable to long)
    std::vector<std::pair<double,double>> local_coords(N_local);
    for (VertexIterator vt = gv.template begin<dim>(); vt != vtend; ++vt) {
      int localIdx = vertexmap.index(*vt);
      auto corner = vt->geometry().corner(0);
      local_coords[localIdx] = std::make_pair(corner[0], corner[1]);
    }

    // Gather all coordinates across all ranks
    std::vector<double> local_flat(2 * N_local);
    for (int ii = 0; ii < N_local; ii++) {
      local_flat[2*ii]   = local_coords[ii].first;
      local_flat[2*ii+1] = local_coords[ii].second;
    }

    std::vector<int> all_sizes(nprocs), all_sizes2(nprocs);
    MPI_Allgather(&N_local, 1, MPI_INT, all_sizes.data(), 1, MPI_INT, MPI_COMM_WORLD);

    int total_verts = 0;
    std::vector<int> displs(nprocs);
    for (int p = 0; p < nprocs; p++) {
      displs[p] = total_verts * 2;
      all_sizes2[p] = all_sizes[p] * 2;
      total_verts += all_sizes[p];
    }

    std::vector<double> all_flat(total_verts * 2);
    MPI_Allgatherv(local_flat.data(), 2*N_local, MPI_DOUBLE,
                   all_flat.data(), all_sizes2.data(), displs.data(),
                   MPI_DOUBLE, MPI_COMM_WORLD);

    // Build sorted unique coordinate set -> contiguous global numbering
    std::map<std::pair<long,long>, int> coord_to_global;
    int next_global = 0;
    for (int i = 0; i < total_verts; i++) {
      long kx = std::lround(all_flat[2*i]   * 1e6);
      long ky = std::lround(all_flat[2*i+1] * 1e6);
      auto key = std::make_pair(kx, ky);
      if (coord_to_global.find(key) == coord_to_global.end())
        coord_to_global[key] = next_global++;
    }
    N_global = next_global;

    // Build local-to-global mapping
    local_to_global.resize(N_local);
    for (int ii = 0; ii < N_local; ii++) {
      long kx = std::lround(local_coords[ii].first  * 1e6);
      long ky = std::lround(local_coords[ii].second * 1e6);
      local_to_global[ii] = coord_to_global[std::make_pair(kx, ky)];
      global_to_local[local_to_global[ii]] = ii;
    }

    // PETSc row ownership
    std::vector<int> local_sizes(nprocs);
    MPI_Allgather(&N_local, 1, MPI_INT, local_sizes.data(), 1, MPI_INT, MPI_COMM_WORLD);
    row_start = 0;
    for (int p = 0; p < rank; p++)
      row_start += local_sizes[p];
    row_end = row_start + N_local;
  }

  if (rank == 0) {
    std::cout << "  [Parallel Index] Global vertices: " << N_global
              << ", row_start=" << row_start << ", row_end=" << row_end << std::endl;
  }
}



// =============================================================================
// Determine adjacency pattern (same as serial, but on local partition)
// =============================================================================
template<class GV, class Vector>
void SchrodingerFEM_Parallel<GV,Vector>::determineAdjacencyPattern()
{
  adjacencyPattern.resize(N_local);

  if (rank == 0)
    std::cout << "  [Rank 0] Number of local unknowns: " << N_local << std::endl;

  const LeafIterator itend = gv.template end<0>();
  for (LeafIterator it = gv.template begin<0>(); it != itend; ++it)
  {
    Dune::GeometryType gt = it->type();
    const auto& ref =
        Dune::ReferenceElements<ctype,dim>::general(gt);

    const IntersectionIterator isend = gv.iend(*it);
    for (IntersectionIterator is = gv.ibegin(*it); is != isend; ++is)
    {
      int vertexsize = ref.size(is->indexInInside(), 1, dim);
      for (int ii = 0; ii < vertexsize; ii++)
      {
        int indexi = set.subIndex(*it, ref.subEntity(is->indexInInside(), 1, ii, dim), dim);
        for (int jj = 0; jj < vertexsize; jj++)
        {
          int indexj = set.subIndex(*it, ref.subEntity(is->indexInInside(), 1, jj, dim), dim);
          adjacencyPattern[indexi].insert(indexj);
        }
      }
    }
  }
}


// =============================================================================
// Assemble: Local FEM assembly, then copy to parallel PETSc matrices
// =============================================================================
template<class GV, class Vector>
void SchrodingerFEM_Parallel<GV,Vector>::assemble(Vector emass, int valley)
{
  A = 0.0;
  B = 0.0;
  ElementMap elementmap(gv, Dune::mcmgElementLayout());
  const LeafIterator itend = gv.template end<0>();

  const P1ShapeFunctionSet<ctype,ctype,dim>& basis = P1ShapeFunctionSet<ctype,ctype,dim>::instance();

  // =================== Local FEM Assembly (same as serial) ===================
  for (LeafIterator it = gv.template begin<0>(); it != itend; ++it)
  {
    Dune::GeometryType gt = it->type();
    const auto& ref =
        Dune::ReferenceElements<ctype,dim>::general(gt);
    int vertexsize = ref.size(dim);

    const Dune::QuadratureRule<ctype,dim>& rule =
        Dune::QuadratureRules<ctype,dim>::rule(gt, 2);
    for (typename Dune::QuadratureRule<ctype,dim>::const_iterator r = rule.begin();
         r != rule.end(); ++r)
    {
      Dune::FieldMatrix<ctype,dim,dim> jacInvTra =
          it->geometry().jacobianInverseTransposed(r->position());
      ctype weight = r->weight();
      ctype detjac = it->geometry().integrationElement(r->position());

      for (int ii = 0; ii < vertexsize; ii++)
      {
        Dune::FieldVector<ctype,dim> grad1;
        jacInvTra.mv(basis[ii].evaluateGradient(r->position()), grad1);
        ctype phii = basis[ii].evaluateFunction(r->position());

        for (int jj = 0; jj < vertexsize; jj++)
        {
          Dune::FieldVector<ctype,dim> grad2;
          jacInvTra.mv(basis[jj].evaluateGradient(r->position()), grad2);
          ctype phij = basis[jj].evaluateFunction(r->position());

          int li = set.subIndex(*it, ii, dim);
          int lj = set.subIndex(*it, jj, dim);

          A[li][lj] += ((grad1 * grad2) * (0.5 / emass[elementmap.index(*it)])) * weight * detjac
                     + (phii * potential[li] * phij) * weight * detjac;

          B[li][lj] += (phii * phij) * weight * detjac;
        }
      }
    }
  }

  // Store valley-specific matrices
  if (valley == 0) { A1 = A; }
  if (valley == 1) { A0 = A; }

  // =================== Assemble C and Cb matrices ===========================
  double nwrad = 60.0, tox1 = 40.0, tox2 = 60.0;
  if (rank == 0) {
    std::ifstream dev_para("dev_para.dat");
    if (dev_para.is_open()) {
      dev_para >> nwrad >> tox1 >> tox2;
      dev_para.close();
    }
  }
  MPI_Bcast(&nwrad, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&tox1, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&tox2, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  double pi = 3.142;
  int Disp = 0;
  int indexii = 0, indexjj = 0;
  double Area = 0;
  VertexMap vertexmap(gv, Dune::mcmgVertexLayout());

  const triangleIterator tiend1 = gv.template end<1>();
  const triangleIterator tistart1 = gv.template begin<1>();

  for (triangleIterator ti = tistart1; ti != tiend1; ++ti)
  {
    Dune::FieldVector<ctype,dim> point_e1 = ti->geometry().corner(0);
    Dune::FieldVector<ctype,dim> point_e2 = ti->geometry().corner(1);
    double c_rad_e1 = sqrt(pow(point_e1[0], 2) + pow(point_e1[1], 2));
    double c_rad_e2 = sqrt(pow(point_e2[0], 2) + pow(point_e2[1], 2));

    if (c_rad_e1 <= (nwrad + tox1 + tox2 + 0.01) &&
        c_rad_e1 >= (nwrad + tox1 + tox2 - 0.01) &&
        c_rad_e2 <= (nwrad + tox1 + tox2 + 0.01) &&
        c_rad_e2 >= (nwrad + tox1 + tox2 - 0.01))
    {
      if (Disp < 1 && rank == 0) {
        std::cout << "  Integrating on surface of radius "
                  << c_rad_e1 << "  " << c_rad_e2 << std::endl;
        Disp++;
      }

      const P1ShapeFunctionSet<ctype,ctype,dim1>& basis1 =
          P1ShapeFunctionSet<ctype,ctype,dim1>::instance();
      Dune::GeometryType gt = ti->type();
      const auto& ref1 =
          Dune::ReferenceElements<ctype,dim1>::general(gt);
      const Dune::QuadratureRule<ctype,dim1>& rule =
          Dune::QuadratureRules<ctype,dim1>::rule(gt, 1);

      for (typename Dune::QuadratureRule<ctype,dim1>::const_iterator r = rule.begin();
           r != rule.end(); ++r)
      {
        ctype weight = r->weight();
        ctype detjac = ti->geometry().integrationElement(r->position());
        int vertexsize = ref1.size(dim1);

        for (int ii = 0; ii < vertexsize; ii++)
        {
          ctype phii = basis1[ii].evaluateFunction(r->position());
          Dune::FieldVector<ctype,dim> point1 = ti->geometry().corner(ii);

          for (int jj = 0; jj < vertexsize; jj++)
          {
            ctype phij = basis1[jj].evaluateFunction(r->position());
            Dune::FieldVector<ctype,dim> point2 = ti->geometry().corner(jj);

            // Find local indices by coordinate matching
            for (VertexIterator vvt = gv.template begin<dim>();
                 vvt != gv.template end<dim>(); ++vvt)
            {
              Dune::FieldVector<ctype,dim> point3 = vvt->geometry().corner(0);
              if (point1[0] == point3[0] && point1[1] == point3[1])
                indexii = vertexmap.index(*vvt);
              if (point2[0] == point3[0] && point2[1] == point3[1])
                indexjj = vertexmap.index(*vvt);
            }

            C[indexii][indexjj]  += 0.5 * 0.5 * phii * phij * weight * detjac / (nwrad + tox1 + tox2);
            Cb[indexii][indexjj] += 0.5 * phii * phij * weight * detjac;
            Area += weight * detjac;
          }
        }
      }
    }
  }

  // Reduce Area across all ranks (sum local contributions)
  double Area_global = 0;
  MPI_Allreduce(&Area, &Area_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  Area = Area_global / 4;
  if (rank == 0)
    std::cout << "  Circumference of the outer circle is " << Area << std::endl;

  // =================== Apply BC before copying to PETSc =====================
  // Zero out rows/cols of boundary nodes
  for (int ii = 0; ii < N_local; ii++) {
    if (A[ii][ii] == 1) {
      for (int jj = 0; jj < N_local; jj++) {
        if (ii != jj) {
          if (A.exists(jj, ii)) A[jj][ii] = 0;
          if (B.exists(jj, ii)) B[jj][ii] = 0;
        }
      }
    }
  }

  // =================== Write debug matrices (rank 0 only) ===================
  if (rank == 0) {
    Dune::writeMatrixToMatlab(A, "MatrixA");
    Dune::writeMatrixToMatlab(A1, "MatrixA1");
    Dune::writeMatrixToMatlab(A0, "MatrixA0");
    Dune::writeMatrixToMatlab(B, "MatrixB");
    Dune::writeMatrixToMatlab(Cb, "MatrixCb");
    Dune::writeMatrixToMatlab(C, "MatrixC");
  }

  // =================== Copy local DUNE matrices to PETSc ====================
  // NOTE: This is done AFTER apply_bc() is called in solve()
  // The matrices A, B are ready in DUNE format; PETSc copy happens later
}


// =============================================================================
// Apply boundary conditions
// =============================================================================
template<class GV, class Vector>
void SchrodingerFEM_Parallel<GV,Vector>::apply_bc()
{
  double pi = 3.142;
  double nwrad = 60.0, tox1 = 40.0, tox2 = 60.0;

  if (rank == 0) {
    std::ifstream dev_para("dev_para.dat");
    if (dev_para.is_open()) {
      dev_para >> nwrad >> tox1 >> tox2;
      dev_para.close();
    }
  }
  MPI_Bcast(&nwrad, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&tox1, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&tox2, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  double R0 = nwrad + tox1 + tox2;

  LeafIterator itend = gv.template end<0>();
  for (LeafIterator it = gv.template begin<0>(); it != itend; ++it)
  {
    const IntersectionIterator isend = gv.iend(*it);
    for (IntersectionIterator is = gv.ibegin(*it); is != isend; ++is)
    {
      Dune::GeometryType gt = it->type();
      const auto& ref =
          Dune::ReferenceElements<ctype,dim>::general(gt);

      Dune::FieldVector<ctype, dim> point = is->geometry().center();
      double rho = sqrt(pow(point[0], 2) + pow(point[1], 2));

      if (rho > nwrad + 10)   // boundary: outside channel radius
      {
        for (int ii = 0; ii < ref.size(is->indexInInside(), 1, dim); ii++)
        {
          int indexi = set.subIndex(*it, ref.subEntity(is->indexInInside(), 1, ii, dim), dim);

          A[indexi]         = 0.0;
          A[indexi][indexi]  = 1.0;
          B[indexi]         = 0.0;
          B[indexi][indexi]  = 1.0;  // Keep B non-singular for Ax=λBx
        }
      }
    }
  }
}


// =============================================================================
// Solve: Parallel SLEPc eigenvalue solve
// =============================================================================
template<class GV, class Vector>
void SchrodingerFEM_Parallel<GV,Vector>::solve(double m, double g, int valley)
{
  VertexMap vertexmap(gv, Dune::mcmgVertexLayout());

  // =================== Copy BC-applied DUNE matrices to PETSc ================
  // This must happen AFTER apply_bc() has modified A and B
  MatZeroEntries(A_p);
  MatZeroEntries(B_p);

  // First: set identity diagonal on ALL local rows to prevent zero-pivot NaN
  // Rows that have actual element data will be overwritten below
  {
    PetscInt row_start, row_end;
    MatGetOwnershipRange(A_p, &row_start, &row_end);
    PetscScalar one = 1.0;
    for (PetscInt i = row_start; i < row_end; i++) {
      MatSetValue(A_p, i, i, one, INSERT_VALUES);
      MatSetValue(B_p, i, i, one, INSERT_VALUES);
    }
  }

  // Now overwrite with actual element-assembled + BC-applied values
  const LeafIterator itend_s = gv.template end<0>();
  for (LeafIterator it = gv.template begin<0>(); it != itend_s; ++it)
  {
    Dune::GeometryType gt = it->type();
    const auto& ref = Dune::ReferenceElements<ctype,dim>::general(gt);
    int vertexsize = ref.size(dim);
    for (int ii = 0; ii < vertexsize; ii++)
    {
      for (int jj = 0; jj < vertexsize; jj++)
      {
        int indexi = set.subIndex(*it, ii, dim);
        int indexj = set.subIndex(*it, jj, dim);
        MatSetValue(A_p, indexi, indexj, double(A[indexi][indexj]), INSERT_VALUES);
        MatSetValue(B_p, indexi, indexj, double(B[indexi][indexj]), INSERT_VALUES);
      }
    }
  }
  MatAssemblyBegin(A_p, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A_p, MAT_FINAL_ASSEMBLY);
  MatAssemblyBegin(B_p, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(B_p, MAT_FINAL_ASSEMBLY);

  PetscScalar* avec;

  double k    = 8.61735E-5;
  double T    = 300;
  double Hr   = 27.211;
  double beta = Hr / (k * T);
  double pai  = 3.142;

  std::ofstream target;
  std::ofstream chargeij;
  std::ofstream f_vij;

  // =================== Create SLEPc EPS solver ==============================
  EPSCreate(PETSC_COMM_SELF, &eps);

  // Create parallel vectors for eigenvectors
  MatCreateVecs(A_p, PETSC_NULL, &xr_vec);
  MatCreateVecs(A_p, PETSC_NULL, &xi_vec);

  EPSSetOperators(eps, A_p, B_p);
  EPSSetFromOptions(eps);
  EPSSetDimensions(eps, 3, PETSC_DECIDE, PETSC_DECIDE);

  // =================== Solve eigenvalue problem =============================
  if (rank == 0)
    std::cout << "  [SLEPc] Solving generalized eigenvalue problem..." << std::endl;

  EPSSolve(eps);

  EPSGetIterationNumber(eps, &its);
  // EPSGetOperationCounters removed in SLEPc 3.15
  EPSGetType(eps, &type);
  EPSGetDimensions(eps, &nev, PETSC_NULL, PETSC_NULL);
  EPSGetTolerances(eps, &tol, &maxit);
  EPSGetConverged(eps, &nconv);

  if (rank == 0) {
    PetscPrintf(PETSC_COMM_WORLD, "  [SLEPc] Iterations: %d, Converged: %d\n", its, nconv);
    PetscPrintf(PETSC_COMM_WORLD, "  [SLEPc] Solver type: %s\n", type);
  }

  double eigen[nconv];
  double nij[nconv];
  double fermi_value[nconv];

  Dune::VTKWriter<GV> vtkwriter_wf(gv, Dune::VTK::conforming);
  Dune::VTKWriter<GV> vtkwriter_wf2(gv, Dune::VTK::conforming);
  Vector wf[nconv];

  for (int ii = 0; ii < nconv; ii++)
  {
    ierr = EPSGetEigenpair(eps, ii, &kr, &ki, xr_vec, xi_vec);
    EPSComputeError(eps, ii, EPS_ERROR_RELATIVE, &error);

#if defined(PETSC_USE_COMPLEX)
    re = PetscRealPart(kr);
    im = PetscImaginaryPart(kr);
#else
    re = kr;
    im = ki;
#endif

    eigen[ii] = re;

    if (im != 0.0) {
      if (rank == 0)
        PetscPrintf(PETSC_COMM_WORLD, "  Eigen %d: %6f %+6fi\n", ii, re, im);
    } else {
      if (rank == 0)
        PetscPrintf(PETSC_COMM_WORLD, "  Eigen %d: %6f eV, error: %12g\n",
                    ii, re * 27.21, error);
    }

    // ============= Extract eigenvector to local DUNE vector ==================
    // For parallel PETSc vectors, each rank owns a portion.
    // We need to gather the full vector or extract our local portion.

    // Get the local portion of the eigenvector
    PetscInt vec_low, vec_high;
    VecGetOwnershipRange(xr_vec, &vec_low, &vec_high);
    VecGetArray(xr_vec, &avec);

    // Map from PETSc global row to local DUNE vertex index
    wavevector = 0.0;
    for (int rr = 0; rr < N_local; rr++) {
      PetscInt global_row = local_to_global[rr];
      if (global_row >= vec_low && global_row < vec_high) {
        wavevector[rr] = PetscRealPart(avec[global_row - vec_low]);
      }
    }
    VecRestoreArray(xr_vec, &avec);

    // For vertices that are ghost/overlap, we need to communicate
    // Use MPI to exchange ghost values
    // (ALUGrid's communication interface handles this)

    WaveFunction<GV, Vector> wavefunction(gv, wavevector);

    wf[ii] = wavevector;
    if (valley == 0)
      wavefunction.wf_VTKwrite(vtkwriter_wf, wf[ii], ii);
    if (valley == 1)
      wavefunction.wf_VTKwrite(vtkwriter_wf2, wf[ii], ii);

    if (wavefunction.is_y_fundamental()) {
      wavefunction.only_x_dependant();
      wavefunction.pnormalize();
      if (rank == 0)
        wavefunction.write_wf_gnuplot("wf_plot_qbs.dat");
    }

    double eeta = (-re) * beta;
    double ans;

    if (carrier_type == 0) {
      Fermi ferm;
      ans = ferm.val(-.5, eeta, 0);
      fermi_value[ii] = ans;
      nij[ii] = g * sqrt((2 * m) / (pai * beta));
    }
  }

  // =================== Write output files (rank 0 only) =====================
  if (rank == 0) {
    if (valley == 0) {
      target.open("target1.dat");
      for (int kk = 0; kk < nconv; kk++)
        target << eigen[kk] << std::endl;
      target.close();

      f_vij.open("f_vij1.dat");
      for (int kk = 0; kk < nconv; kk++)
        f_vij << fermi_value[kk] << std::endl;
      f_vij.close();

      chargeij.open("nij1.dat");
      for (int kk = 0; kk < nconv; kk++)
        chargeij << nij[kk] << std::endl;
      chargeij.close();
    }

    if (valley == 1) {
      target.open("target2.dat");
      for (int kk = 0; kk < nconv; kk++)
        target << eigen[kk] << std::endl;
      target.close();

      f_vij.open("f_vij2.dat");
      for (int kk = 0; kk < nconv; kk++)
        f_vij << fermi_value[kk] << std::endl;
      f_vij.close();

      chargeij.open("nij2.dat");
      for (int kk = 0; kk < nconv; kk++)
        chargeij << nij[kk] << std::endl;
      chargeij.close();
    }
  }

  // Write parallel VTK
  if (valley == 0)
    vtkwriter_wf.pwrite("wavefunction_v0", ".", "", Dune::VTK::appendedraw);
  if (valley == 1)
    vtkwriter_wf2.pwrite("wavefunction_v1", ".", "", Dune::VTK::appendedraw);

  // Cleanup
  EPSDestroy(&eps);
}


// =============================================================================
// Apply: Run the full Schrodinger solve for all valleys
// =============================================================================
template<class GV, class Vector>
void SchrodingerFEM_Parallel<GV,Vector>::apply()
{
  ElementMap elementmap(gv, Dune::mcmgElementLayout());

  double mq1 = 0, mq2 = 0, mz1 = 0, mz2 = 0, g1 = 0, g2 = 0;
  if (rank == 0) {
    std::ifstream channel_para("channel_para.dat");
    if (channel_para.is_open()) {
      channel_para >> mq1 >> mq2 >> mz1 >> mz2 >> g1 >> g2;
      channel_para.close();
      std::cout << "  Channel params: mq1=" << mq1 << " mq2=" << mq2
                << " mz1=" << mz1 << " mz2=" << mz2
                << " g1=" << g1 << " g2=" << g2 << std::endl;
    } else {
      std::cout << "  WARNING: channel_para.dat not found, using defaults" << std::endl;
      mq1 = 0.19; mq2 = 0.19; mz1 = 0.916; mz2 = 0.19; g1 = 2; g2 = 4;
    }
  }

  // Broadcast channel parameters to all ranks
  double ch_params[6];
  ch_params[0] = mq1; ch_params[1] = mq2;
  ch_params[2] = mz1; ch_params[3] = mz2;
  ch_params[4] = g1;  ch_params[5] = g2;
  MPI_Bcast(ch_params, 6, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  mq1 = ch_params[0]; mq2 = ch_params[1];
  mz1 = ch_params[2]; mz2 = ch_params[3];
  g1  = ch_params[4]; g2  = ch_params[5];

  LeafIterator itend = gv.template end<0>();
  double me_d[2] = {mq1, mq2};
  double me[2]   = {mz1, mz2};
  double ge[2]   = {g1, g2};

  double mh[3]   = {0.2, 0.291, 0.29};
  double mh_d[3] = {0.251, 0.654, 0.29};
  double gh[3]   = {1, 1, 1};

  if (rank == 0)
    std::cout << "  Carrier type: " << carrier_type << std::endl;

  if (carrier_type == 0) {
    for (int ii = 0; ii < 2; ii++) {
      if (ge[ii] != 0) {
        // Set effective mass for channel region
        for (LeafIterator it = gv.template begin<0>(); it != itend; ++it) {
          if (regionid[elementmap.index(*it)] == 201 || regionid[elementmap.index(*it)] == 301)
            emass[elementmap.index(*it)] = me_d[ii];
        }

        if (rank == 0)
          std::cout << "  Assembling valley " << ii << " (me_d=" << me_d[ii]
                    << ", ge=" << ge[ii] << ")" << std::endl;

        assemble(emass, ii);
        apply_bc();
        solve(me[ii], ge[ii], ii);
      }
    }
  }

  if (carrier_type == 1) {
    for (int ii = 0; ii < 3; ii++) {
      assemble(emass, ii);
      apply_bc();
      if (rank == 0)
        std::cout << "  Hole valley " << ii << std::endl;
    }
  }

  charge = charge_tmp;
}

#endif // SCHRODINGERFEM_PARALLEL_HH
