// =============================================================================
// Parallel MPI Schrodinger FEM Solver for Flash Memory (High-K dielectric stack)
// Adapted from serial schrodingerfem_12.hh
// Key changes from serial:
//   1. PETSc matrices: MATSEQAIJ on PETSC_COMM_SELF per rank
//   2. SLEPc EPS solver: per-rank eigenvalue solve
//   3. Local-to-global index mapping for parallel assembly
//   4. MPI broadcast for file I/O (rank 0 reads)
//   5. Boundary conditions: x-based (gate side / substrate side)
//   6. DUNE 2.10 API: ReferenceElements, mcmgMapper, leafGridView()
//   7. B[i][i]=1 for boundary nodes (non-singular generalized eigenproblem)
//   8. PETSc copy AFTER apply_bc (correct ordering)
//   9. Identity pre-fill for ghost rows (no zero-pivot NaN)
// =============================================================================

#ifndef SCHRODINGERFEM_HIGHK_PARALLEL_HH
#define SCHRODINGERFEM_HIGHK_PARALLEL_HH

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
#include "shapefunctions.hh"
#include <dune/geometry/quadraturerules.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/io.hh>
#include "wavefunction.hh"

template<class GV, class Vector>
class SchrodingerFEM_HighK_Parallel
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
    int N_local;
    int N_global;
    int rank, nprocs;
    std::vector<PetscInt> local_to_global;
    std::map<int, int> global_to_local;
    PetscInt row_start, row_end;

    // Sparsity pattern
    std::vector<std::set<int> > adjacencyPattern;

    // Local DUNE matrices
    Matrix A, A1, A0, B, Cb, C;

    // PETSc matrices
    Mat A_p, B_p;

    // SLEPc eigenvalue solver
    EPS            eps;
    EPSType  type;
    PetscReal      error, tol, re, im, mod_xr;
    PetscScalar    kr, ki, xr_val;
    PetscErrorCode ierr;
    PetscInt       nev, maxit, i, its, lits, nconv;
    Vec            xr_vec, xi_vec;

    int carrier_type;
    int car_type;
    double Ef;

    void determineAdjacencyPattern();
    void buildParallelIndexMapping();
    void assemble(Vector mass, int valley);
    void apply_bc();
    void solve(double m, double g, int valley);

  public:
    SchrodingerFEM_HighK_Parallel(const GV& gv_, Vector& charge_, Vector& potential_,
      Vector& wavefunction_bc, int carrier_type_, double Ef_,
      Vector& regionid_, Vector& emass_, int& argc, char**& argv);
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
// Constructor
// =============================================================================
template<class GV, class Vector>
SchrodingerFEM_HighK_Parallel<GV,Vector>::SchrodingerFEM_HighK_Parallel(
    const GV& gv_, Vector& charge_, Vector& potential_,
    Vector& wavefunction_bc_, int carrier_type_, double Ef_,
    Vector& regionid_, Vector& emass_, int& argc, char**& argv)
  : gv(gv_), charge(charge_), potential(potential_),
    wavefunction_bc(wavefunction_bc_), emass(emass_), regionid(regionid_),
    set(gv.indexSet())
{
  carrier_type = carrier_type_;
  Ef = Ef_;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  N_local = gv.size(dim);

  buildParallelIndexMapping();

  if (rank == 0)
    std::cout << "  [SchrodingerFEM_HighK] N_global = " << N_global
              << ", N_local(rank 0) = " << N_local << std::endl;

  determineAdjacencyPattern();

  // =================== PETSc/SLEPc Initialization ==========================
  static char help[] = "Parallel Flash Memory Schrodinger Solver";
  SlepcInitialize(&argc, &argv, (char*)0, help);

  // Create PETSc matrices (PETSC_COMM_SELF per rank)
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

  // Setup DUNE matrices
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


// =============================================================================
// Build local-to-global index mapping
// =============================================================================
template<class GV, class Vector>
void SchrodingerFEM_HighK_Parallel<GV,Vector>::buildParallelIndexMapping()
{
  VertexMap vertexmap(gv, Dune::mcmgVertexLayout());

  // Gather vertex coordinates from all ranks
  std::vector<std::pair<double,double>> local_coords(N_local);
  for (auto vt = gv.template begin<dim>(); vt != gv.template end<dim>(); ++vt) {
    int li = vertexmap.index(*vt);
    local_coords[li] = std::make_pair(vt->geometry().corner(0)[0],
                                       vt->geometry().corner(0)[1]);
  }

  // Exchange all coordinates via MPI
  std::vector<int> all_sizes(nprocs);
  MPI_Allgather(&N_local, 1, MPI_INT, all_sizes.data(), 1, MPI_INT, MPI_COMM_WORLD);

  int total_verts = 0;
  std::vector<int> displs(nprocs);
  for (int p = 0; p < nprocs; p++) {
    displs[p] = total_verts * 2;
    total_verts += all_sizes[p];
  }

  std::vector<double> local_xy(N_local * 2);
  for (int ii = 0; ii < N_local; ii++) {
    local_xy[2*ii]   = local_coords[ii].first;
    local_xy[2*ii+1] = local_coords[ii].second;
  }

  std::vector<double> all_xy(total_verts * 2);
  std::vector<int> sendcounts(nprocs);
  for (int p = 0; p < nprocs; p++) sendcounts[p] = all_sizes[p] * 2;

  MPI_Allgatherv(local_xy.data(), N_local * 2, MPI_DOUBLE,
                 all_xy.data(), sendcounts.data(), displs.data(),
                 MPI_DOUBLE, MPI_COMM_WORLD);

  // Build unique coordinate-to-global-ID map
  std::map<std::pair<long,long>, int> coord_to_global;
  N_global = 0;
  for (int v = 0; v < total_verts; v++) {
    long kx = std::lround(all_xy[2*v]   * 1e6);
    long ky = std::lround(all_xy[2*v+1] * 1e6);
    auto key = std::make_pair(kx, ky);
    if (coord_to_global.find(key) == coord_to_global.end()) {
      coord_to_global[key] = N_global++;
    }
  }

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

  if (rank == 0) {
    std::cout << "  [Parallel Index] Global vertices: " << N_global
              << ", row_start=" << row_start << ", row_end=" << row_end << std::endl;
  }
}


// =============================================================================
// Determine adjacency pattern
// =============================================================================
template<class GV, class Vector>
void SchrodingerFEM_HighK_Parallel<GV,Vector>::determineAdjacencyPattern()
{
  adjacencyPattern.resize(N_local);

  if (rank == 0)
    std::cout << "  [Rank 0] Number of local unknowns: " << N_local << std::endl;

  const LeafIterator itend = gv.template end<0>();
  for (LeafIterator it = gv.template begin<0>(); it != itend; ++it)
  {
    Dune::GeometryType gt = it->type();
    const auto& ref = Dune::ReferenceElements<ctype,dim>::general(gt);

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
// Assemble: Local FEM assembly
// =============================================================================
template<class GV, class Vector>
void SchrodingerFEM_HighK_Parallel<GV,Vector>::assemble(Vector emass, int valley)
{
  A = 0.0;
  B = 0.0;
  ElementMap elementmap(gv, Dune::mcmgElementLayout());
  const LeafIterator itend = gv.template end<0>();

  const P1ShapeFunctionSet<ctype,ctype,dim>& basis = P1ShapeFunctionSet<ctype,ctype,dim>::instance();

  for (LeafIterator it = gv.template begin<0>(); it != itend; ++it)
  {
    Dune::GeometryType gt = it->type();
    const auto& ref = Dune::ReferenceElements<ctype,dim>::general(gt);
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

          // H = -hbar^2/(2m) nabla^2 + V(x)
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

  // Debug output (rank 0 only)
  if (rank == 0) {
    Dune::writeMatrixToMatlab(A, "MatrixA");
    Dune::writeMatrixToMatlab(B, "MatrixB");
  }
}


// =============================================================================
// Apply boundary conditions for flash memory
// Boundary: nodes on the outer boundary of the mesh (is->boundary())
// For flash memory: wavefunction = 0 at device boundaries
// =============================================================================
template<class GV, class Vector>
void SchrodingerFEM_HighK_Parallel<GV,Vector>::apply_bc()
{
  LeafIterator itend = gv.template end<0>();
  for (LeafIterator it = gv.template begin<0>(); it != itend; ++it)
  {
    const IntersectionIterator isend = gv.iend(*it);
    for (IntersectionIterator is = gv.ibegin(*it); is != isend; ++is)
    {
      if (is->boundary())
      {
        Dune::GeometryType gt = it->type();
        const auto& ref = Dune::ReferenceElements<ctype,dim>::general(gt);

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
// Solve: Copy DUNE matrices to PETSc, run SLEPc eigenvalue solve
// =============================================================================
template<class GV, class Vector>
void SchrodingerFEM_HighK_Parallel<GV,Vector>::solve(double m, double g, int valley)
{
  VertexMap vertexmap(gv, Dune::mcmgVertexLayout());

  // =================== Copy BC-applied DUNE matrices to PETSc ================
  MatZeroEntries(A_p);
  MatZeroEntries(B_p);

  // Pre-fill identity on ALL local rows (prevents zero-pivot NaN on ghost rows)
  {
    PetscInt row_s, row_e;
    MatGetOwnershipRange(A_p, &row_s, &row_e);
    PetscScalar one = 1.0;
    for (PetscInt i = row_s; i < row_e; i++) {
      MatSetValue(A_p, i, i, one, INSERT_VALUES);
      MatSetValue(B_p, i, i, one, INSERT_VALUES);
    }
  }

  // Overwrite with actual element-assembled + BC-applied values
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
    PetscInt vec_low, vec_high;
    VecGetOwnershipRange(xr_vec, &vec_low, &vec_high);
    VecGetArray(xr_vec, &avec);

    wavevector = 0.0;
    for (int rr = 0; rr < N_local; rr++) {
      PetscInt global_row = local_to_global[rr];
      if (global_row >= vec_low && global_row < vec_high) {
        wavevector[rr] = PetscRealPart(avec[global_row - vec_low]);
      }
    }
    VecRestoreArray(xr_vec, &avec);

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
        wavefunction.write_wf_gnuplot("wf_plot_flash.dat");
    }

    double eeta = (-re) * beta;
    double ans;

    if (carrier_type == 0) {
      Fermi ferm;
      ans = ferm.val(-.5, eeta, 0);
      fermi_value[ii] = ans;
      nij[ii] = g * sqrt((2 * m) / (pai * beta));
    }
    if (carrier_type == 1) {
      Fermi ferm;
      ans = ferm.val(-.5, eeta, 0);
      fermi_value[ii] = ans;
      nij[ii] = g * sqrt((2 * m) / (pai * beta));
    }
  }

  // =================== Compute charge density ================================
  charge_tmp = 0.0;
  for (int ii = 0; ii < nconv; ii++) {
    for (int rr = 0; rr < N_local; rr++) {
      if (carrier_type == 0)
        charge_tmp[rr] += nij[ii] * fermi_value[ii] * wf[ii][rr] * wf[ii][rr];
      if (carrier_type == 1)
        charge_tmp[rr] += nij[ii] * fermi_value[ii] * wf[ii][rr] * wf[ii][rr];
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
    vtkwriter_wf.pwrite("wavefunction_flash_v0", ".", "", Dune::VTK::appendedraw);
  if (valley == 1)
    vtkwriter_wf2.pwrite("wavefunction_flash_v1", ".", "", Dune::VTK::appendedraw);

  // Cleanup
  EPSDestroy(&eps);
}


// =============================================================================
// Apply: Run the full Schrodinger solve for all valleys
// Flash memory uses the same valley structure as serial code
// =============================================================================
template<class GV, class Vector>
void SchrodingerFEM_HighK_Parallel<GV,Vector>::apply()
{
  ElementMap elementmap(gv, Dune::mcmgElementLayout());

  // Read channel parameters from file (rank 0 reads, broadcast)
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
        for (LeafIterator it = gv.template begin<0>(); it != itend; ++it) {
          int rid = (int)regionid[elementmap.index(*it)];
          // Silicon substrate regions get the effective mass
          if (rid == 200 || rid == 201 || rid == 202)
            emass[elementmap.index(*it)] = me_d[ii];
        }

        if (rank == 0)
          std::cout << "  Assembling electron valley " << ii << " (me_d=" << me_d[ii]
                    << ", ge=" << ge[ii] << ")" << std::endl;

        assemble(emass, ii);
        apply_bc();
        solve(me[ii], ge[ii], ii);
      }
    }
  }

  if (carrier_type == 1) {
    for (int ii = 0; ii < 3; ii++) {
      for (LeafIterator it = gv.template begin<0>(); it != itend; ++it) {
        int rid = (int)regionid[elementmap.index(*it)];
        if (rid == 200 || rid == 201 || rid == 202)
          emass[elementmap.index(*it)] = mh_d[ii];
      }

      if (rank == 0)
        std::cout << "  Assembling hole valley " << ii << " (mh_d=" << mh_d[ii]
                  << ", gh=" << gh[ii] << ")" << std::endl;

      assemble(emass, ii);
      apply_bc();
      solve(mh[ii], gh[ii], ii);
    }
  }

  charge = charge_tmp;
}


// =============================================================================
// Accessor methods
// =============================================================================
template<class GV, class Vector>
double SchrodingerFEM_HighK_Parallel<GV,Vector>::geteigen(int i)
{
  return real_DOT[i];
}

template<class GV, class Vector>
Vector SchrodingerFEM_HighK_Parallel<GV,Vector>::getwave(int i)
{
  Vector temp(N_local);
  for (int j = 0; j < N_local; j++)
    temp[j] = wave_DOT[i][j];
  return temp;
}

template<class GV, class Vector>
int SchrodingerFEM_HighK_Parallel<GV,Vector>::getnconv()
{
  return nconv;
}

#endif // SCHRODINGERFEM_HIGHK_PARALLEL_HH
