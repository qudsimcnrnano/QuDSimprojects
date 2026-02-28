// poisson3d_petsc_slepc_mpi.cc
// ============================================================================
// 3D Poisson FEM (P1) solver on a tetrahedral Gmsh mesh using DUNE.
// Uses PETSc for parallel linear algebra and KSP solvers.
// SLEPc is initialized for potential eigenvalue extensions.
// Fully MPI-parallel with proper global index mapping via ISLocalToGlobalMapping.
//
// Solves:  -∇²u(x,y,z) = ρ(x,y,z)   in Ω = [0,16]³
//
// Boundary conditions:
//   Dirichlet: u = g(x,y,z) on ∂Ω
//     - u = 1.0   on the face z = 0   (grounded plate at potential 1V)
//     - u = -1.0  on the face z = 16  (grounded plate at potential -1V)
//     - u = 0.0   on all other faces   (homogeneous Dirichlet)
//
// Charge density (RHS):
//   ρ(x,y,z) = 100 · exp(-((x-8)² + (y-8)² + (z-8)²) / (2·σ²))
//   A Gaussian charge blob centered at (8,8,8) with σ = 2.0
//   This represents a localized point-like charge distribution.
//
// Compile (example):
//   mpicxx -O2 -o poisson3d poisson3d_petsc_slepc_mpi.cc \
//     $(pkg-config --cflags --libs dune-grid dune-common dune-geometry) \
//     $(pkg-config --cflags --libs PETSc slepc) -lm
//
// Run:
//   mpirun -np 4 ./poisson3d poisson3d.msh
// ============================================================================

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <sstream>
#include <iostream>
#include <algorithm>
#include <map>
#include <vector>
#include <set>
#include <cmath>
#include <cstdlib>
#include <string>
#include <numeric>

// ─── DUNE headers ───────────────────────────────────────────────────────────
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/yaspgrid.hh>

// ─── ALUGrid for 3D unstructured tetrahedral grids ──────────────────────────
#include <dune/alugrid/3d/alugrid.hh>

// ─── DUNE ISTL (optional, for adjacency patterns) ──────────────────────────
#include <dune/istl/bvector.hh>
#include <dune/istl/bcrsmatrix.hh>

// ─── PETSc / SLEPc headers ─────────────────────────────────────────────────
#include "petsc.h"
#include "petscmat.h"
#include "petscksp.h"
#include "petscis.h"
#include "petscvec.h"
#include "slepc.h"
#include "slepceps.h"
#include "slepcsys.h"

// ─── Local headers ──────────────────────────────────────────────────────────
#include "shapefunctions.hh"
#include "qudsim_timer.hh"


// ============================================================================
// Charge density function: Gaussian blob centered at (8,8,8)
// ρ(x,y,z) = amplitude * exp(-|x - x₀|² / (2σ²))
// ============================================================================
template<class ctype, int dim>
class GaussianCharge {
public:
  ctype operator()(const Dune::FieldVector<ctype, dim>& x) const
  {
    static_assert(dim == 3, "GaussianCharge is designed for 3D");
    const ctype amplitude = 100.0;
    const ctype sigma     = 2.0;
    const ctype x0 = 8.0, y0 = 8.0, z0 = 8.0;

    ctype r2 = (x[0] - x0) * (x[0] - x0)
             + (x[1] - y0) * (x[1] - y0)
             + (x[2] - z0) * (x[2] - z0);

    return amplitude * std::exp(-r2 / (2.0 * sigma * sigma));
  }
};


// ============================================================================
// Dirichlet boundary condition function
// g(x,y,z) = +1.0 on z=0 face, -1.0 on z=16 face, 0.0 elsewhere
// ============================================================================
template<class ctype, int dim>
class DirichletBC {
public:
  ctype operator()(const Dune::FieldVector<ctype, dim>& x) const
  {
    static_assert(dim == 3, "DirichletBC is designed for 3D");
    const ctype tol = 1e-10;
    const ctype zmin = 0.0;
    const ctype zmax = 16.0;

    if (std::abs(x[2] - zmin) < tol)
      return 1.0;    // bottom face: u = +1
    else if (std::abs(x[2] - zmax) < tol)
      return -1.0;   // top face:    u = -1
    else
      return 0.0;    // all other boundary faces
  }
};


// ============================================================================
// PETSc-based P1 FEM assembler/solver for 3D Poisson equation
// ============================================================================
template<class GV, class RHSFunc, class BCFunc>
class P1PoissonPetsc3D
{
public:
  static const int dim = GV::dimension;
  typedef typename GV::ctype ctype;

private:
  using LeafIterator           = typename GV::template Codim<0>::Iterator;
  using IntersectionIterator   = typename GV::IntersectionIterator;
  using LeafIndexSet           = typename GV::IndexSet;
  using JacobianInverseTransposed =
      typename GV::template Codim<0>::Geometry::JacobianInverseTransposed;

  const GV&      gv;
  const RHSFunc& rhsFunc;    // charge density / RHS function
  const BCFunc&  bcFunc;     // Dirichlet BC function

public:
  Mat A;                      // PETSc parallel stiffness matrix
  Vec b;                      // PETSc parallel RHS vector
  Vec u;                      // PETSc parallel solution vector
  KSP ksp;
  PetscErrorCode ierr;

  std::vector<std::set<int>> adjacencyPattern;
  std::vector<PetscInt> l2g;  // local-to-global index mapping
  PetscInt globalN;           // total number of unique global vertices

  P1PoissonPetsc3D(const GV& gv_, const RHSFunc& rhs_, const BCFunc& bc_)
    : gv(gv_), rhsFunc(rhs_), bcFunc(bc_), globalN(0),
      A(nullptr), b(nullptr), u(nullptr) {}

  void buildGlobalIndexMap(int mpiRank, int mpiSize);
  void determineAdjacencyPattern();
  void assemble(int& argc, char**& argv);
  void solve(int& argc, char**& argv);

  void destroyPetsc()
  {
    if (A) MatDestroy(&A);
    if (b) VecDestroy(&b);
    if (u) VecDestroy(&u);
    A = nullptr; b = nullptr; u = nullptr;
  }
};


// ============================================================================
// Build global index mapping (coordinate-based)
// ============================================================================
template<class GV, class RHSFunc, class BCFunc>
void P1PoissonPetsc3D<GV, RHSFunc, BCFunc>::buildGlobalIndexMap(int mpiRank, int mpiSize)
{
  const int N = gv.size(dim);
  const LeafIndexSet& set = gv.indexSet();

  // Step 1: Collect coordinates of all local vertices
  std::vector<double> myCoords(N * dim);
  for (auto it = gv.template begin<dim>(); it != gv.template end<dim>(); ++it)
  {
    int localIdx = set.index(*it);
    auto pos = it->geometry().center();
    for (int d = 0; d < dim; ++d)
      myCoords[localIdx * dim + d] = pos[d];
  }

  // Step 2: Gather all vertex coordinates from all ranks
  int myCount = N * dim;
  std::vector<int> allCounts(mpiSize);
  MPI_Allgather(&myCount, 1, MPI_INT, allCounts.data(), 1, MPI_INT, PETSC_COMM_WORLD);

  std::vector<int> displs(mpiSize, 0);
  for (int p = 1; p < mpiSize; ++p)
    displs[p] = displs[p - 1] + allCounts[p - 1];
  int totalDoubles = displs[mpiSize - 1] + allCounts[mpiSize - 1];

  std::vector<double> allCoords(totalDoubles);
  MPI_Allgatherv(myCoords.data(), myCount, MPI_DOUBLE,
                 allCoords.data(), allCounts.data(), displs.data(), MPI_DOUBLE,
                 PETSC_COMM_WORLD);

  // Step 3: Build sorted unique list of all vertex coordinate tuples
  int totalVerts = totalDoubles / dim;

  struct CoordTuple {
    double c[3];
    bool operator<(const CoordTuple& o) const {
      for (int d = 0; d < 3; ++d) {
        if (c[d] < o.c[d] - 1e-14) return true;
        if (c[d] > o.c[d] + 1e-14) return false;
      }
      return false;
    }
    bool operator==(const CoordTuple& o) const {
      for (int d = 0; d < 3; ++d)
        if (std::abs(c[d] - o.c[d]) > 1e-14) return false;
      return true;
    }
  };

  std::vector<CoordTuple> allVerts(totalVerts);
  for (int v = 0; v < totalVerts; ++v) {
    allVerts[v].c[0] = allVerts[v].c[1] = allVerts[v].c[2] = 0.0;
    for (int d = 0; d < dim; ++d)
      allVerts[v].c[d] = allCoords[v * dim + d];
  }

  std::sort(allVerts.begin(), allVerts.end());
  allVerts.erase(std::unique(allVerts.begin(), allVerts.end()), allVerts.end());
  globalN = (PetscInt)allVerts.size();

  // Step 4: Build map from coordinate -> contiguous global index
  std::map<CoordTuple, PetscInt> coordToGlobal;
  for (PetscInt i = 0; i < globalN; ++i)
    coordToGlobal[allVerts[i]] = i;

  // Step 5: Build local -> global mapping
  l2g.resize(N);
  for (int i = 0; i < N; ++i) {
    CoordTuple ct;
    ct.c[0] = ct.c[1] = ct.c[2] = 0.0;
    for (int d = 0; d < dim; ++d)
      ct.c[d] = myCoords[i * dim + d];

    auto findIt = coordToGlobal.find(ct);
    if (findIt != coordToGlobal.end()) {
      l2g[i] = findIt->second;
    } else {
      // Fallback: linear search with tolerance
      bool found = false;
      for (PetscInt g = 0; g < globalN; ++g) {
        bool match = true;
        for (int d = 0; d < dim; ++d) {
          if (std::abs(ct.c[d] - allVerts[g].c[d]) > 1e-12) { match = false; break; }
        }
        if (match) { l2g[i] = g; found = true; break; }
      }
      if (!found) {
        std::cerr << "ERROR: vertex " << i << " on rank " << mpiRank
                  << " not found in global list!" << std::endl;
        MPI_Abort(PETSC_COMM_WORLD, 1);
      }
    }
  }

  if (mpiRank == 0)
    std::cout << "Global index map built: " << globalN
              << " unique global vertices" << std::endl;
}


// ============================================================================
// Determine adjacency (sparsity) pattern from mesh connectivity
// ============================================================================
template<class GV, class RHSFunc, class BCFunc>
void P1PoissonPetsc3D<GV, RHSFunc, BCFunc>::determineAdjacencyPattern()
{
  const int N = gv.size(dim);
  adjacencyPattern.clear();
  adjacencyPattern.resize(N);

  const LeafIndexSet& set = gv.indexSet();

  // For 3D tetrahedra: iterate over elements and connect all vertex pairs
  for (auto it = gv.template begin<0>(); it != gv.template end<0>(); ++it)
  {
    auto geo = it->geometry();
    auto ref = Dune::referenceElement(geo);
    int nVerts = ref.size(dim);  // should be 4 for tetrahedra

    for (int i = 0; i < nVerts; ++i)
    {
      int indexi = set.subIndex(*it, i, dim);
      for (int j = 0; j < nVerts; ++j)
      {
        int indexj = set.subIndex(*it, j, dim);
        adjacencyPattern[indexi].insert(indexj);
      }
    }
  }
}


// ============================================================================
// Assemble: stiffness matrix + RHS vector + Dirichlet BCs
// ============================================================================
template<class GV, class RHSFunc, class BCFunc>
void P1PoissonPetsc3D<GV, RHSFunc, BCFunc>::assemble(int& argc, char**& argv)
{
  const int N = gv.size(dim);  // local vertex count (including ghosts)
  const LeafIndexSet& set = gv.indexSet();

  static char help[] = "3D Poisson FEM Solver";

  PetscFunctionBeginUser;
  ierr = PetscInitialize(&argc, &argv, (char*)0, help);

  PetscMPIInt rank, size;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  MPI_Comm_size(PETSC_COMM_WORLD, &size);

  std::cout << "Rank " << rank << ": N=" << N << " local vertices, "
            << gv.size(0) << " local elements (tetrahedra)" << std::endl;

  // ── Build global index mapping ──
  buildGlobalIndexMap(rank, size);

  // ── Determine adjacency pattern ──
  determineAdjacencyPattern();

  // ── P1 basis for 3D tetrahedra (4 shape functions) ──
  P1ShapeFunctionSet<ctype, ctype, dim> basis =
      P1ShapeFunctionSet<ctype, ctype, dim>::instance();

  // ══════════════════════════════════════════════════════════════════════════
  // Use MatSetLocalToGlobalMapping + MatSetValuesLocal
  // Insert values using LOCAL indices; PETSc maps to global via l2g.
  // ══════════════════════════════════════════════════════════════════════════

  // ── Create ISLocalToGlobalMapping ──
  ISLocalToGlobalMapping lgmap;
  ISLocalToGlobalMappingCreate(PETSC_COMM_WORLD, 1, N, l2g.data(),
                               PETSC_COPY_VALUES, &lgmap);

  // ── Create PETSc vectors ──
  VecCreate(PETSC_COMM_WORLD, &b);
  VecSetSizes(b, PETSC_DECIDE, globalN);
  VecSetFromOptions(b);
  VecSetLocalToGlobalMapping(b, lgmap);
  VecDuplicate(b, &u);
  VecSet(b, 0.0);

  if (rank == 0)
    std::cout << "PETSc vectors created (globalN=" << globalN << ")" << std::endl;

  // ── Create PETSc matrix ──
  MatCreate(PETSC_COMM_WORLD, &A);
  MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, globalN, globalN);
  MatSetFromOptions(A);
  MatSetLocalToGlobalMapping(A, lgmap, lgmap);

  // Preallocation: estimate max nonzeros per row
  {
    int maxNnz = 0;
    for (int i = 0; i < N; ++i) {
      int sz = (int)adjacencyPattern[i].size();
      if (sz > maxNnz) maxNnz = sz;
    }
    // In 3D, typical P1 tet vertex has ~15-25 neighbors; add safety margin
    maxNnz = std::max(maxNnz, 30);
    MatMPIAIJSetPreallocation(A, maxNnz, NULL, maxNnz, NULL);
    MatSeqAIJSetPreallocation(A, maxNnz, NULL);
  }

  // Allow extra nonzeros without error (safety net)
  MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);

  if (rank == 0)
    std::cout << "PETSc matrix created with local-to-global mapping." << std::endl;

  // ══════════════════════════════════════════════════════════════════════════
  // Element loop: assemble stiffness matrix K and RHS vector f
  // For each tetrahedron, compute 4×4 element stiffness and 4×1 load vector
  // ══════════════════════════════════════════════════════════════════════════
  for (auto it = gv.template begin<0>(); it != gv.template end<0>(); ++it)
  {
    auto geo = it->geometry();
    auto ref = Dune::referenceElement(geo);
    int vertexsize = ref.size(dim);  // 4 for tetrahedra

    // Collect local vertex indices for this element
    std::vector<PetscInt> localRows(vertexsize);
    for (int i = 0; i < vertexsize; ++i)
      localRows[i] = (PetscInt)set.subIndex(*it, i, dim);

    // ── Element stiffness matrix: Ke[i][j] = ∫ ∇φ_i · ∇φ_j dΩ ──
    std::vector<PetscScalar> Ke(vertexsize * vertexsize, 0.0);

    // Quadrature order 1 is exact for linear gradients on tetrahedra
    const auto& rule = Dune::QuadratureRules<ctype, dim>::rule(it->type(), 1);
    for (auto r = rule.begin(); r != rule.end(); ++r)
    {
      JacobianInverseTransposed JinvT =
          geo.jacobianInverseTransposed(r->position());

      ctype w    = r->weight();
      ctype detJ = geo.integrationElement(r->position());

      for (int i = 0; i < vertexsize; ++i)
      {
        Dune::FieldVector<ctype, dim> grad_i;
        JinvT.mv(basis[i].evaluateGradient(r->position()), grad_i);

        for (int j = 0; j < vertexsize; ++j)
        {
          Dune::FieldVector<ctype, dim> grad_j;
          JinvT.mv(basis[j].evaluateGradient(r->position()), grad_j);

          Ke[i * vertexsize + j] += (PetscScalar)((grad_i * grad_j) * w * detJ);
        }
      }
    }

    // Insert element stiffness into PETSc using LOCAL indices
    MatSetValuesLocal(A, vertexsize, localRows.data(),
                         vertexsize, localRows.data(),
                         Ke.data(), ADD_VALUES);

    // ── Element RHS vector: be[i] = ∫ φ_i · ρ(x) dΩ ──
    std::vector<PetscScalar> be(vertexsize, 0.0);
    // Use quadrature order 2 for the RHS (charge density is non-polynomial)
    const auto& rule2 = Dune::QuadratureRules<ctype, dim>::rule(it->type(), 2);
    for (auto r = rule2.begin(); r != rule2.end(); ++r)
    {
      ctype w    = r->weight();
      ctype detJ = geo.integrationElement(r->position());
      auto xGlobal = geo.global(r->position());

      for (int i = 0; i < vertexsize; ++i)
      {
        ctype phi_i = basis[i].evaluateFunction(r->position());
        be[i] += (PetscScalar)(phi_i * rhsFunc(xGlobal) * w * detJ);
      }
    }

    VecSetValuesLocal(b, vertexsize, localRows.data(), be.data(), ADD_VALUES);
  }

  if (rank == 0)
    std::cout << "Element loop complete. Assembling PETSc objects..." << std::endl;

  // ── Assemble PETSc matrix and vector ──
  MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
  VecAssemblyBegin(b);
  VecAssemblyEnd(b);

  if (rank == 0)
    std::cout << "Matrix and vector assembled." << std::endl;

  // ══════════════════════════════════════════════════════════════════════════
  // Apply inhomogeneous Dirichlet boundary conditions
  //
  // For inhomogeneous BCs (u = g on ∂Ω where g ≠ 0):
  //   1. Identify boundary DOFs and their prescribed values
  //   2. Zero the BC rows in A, put 1.0 on diagonal
  //   3. Set b[i] = g(x_i) at boundary DOFs
  //
  // Note: For a more accurate "lift" approach, one would modify the RHS
  // to account for the contribution of boundary values to interior DOFs.
  // Here we use the simpler row-zeroing approach which is standard.
  // ══════════════════════════════════════════════════════════════════════════

  // Collect boundary DOFs: global indices and their prescribed values
  std::vector<PetscInt> bcGlobalRows;
  std::map<PetscInt, PetscScalar> bcValues;  // global index -> prescribed value

  for (auto it = gv.template begin<0>(); it != gv.template end<0>(); ++it)
  {
    auto ref = Dune::referenceElement(it->geometry());

    for (auto is = gv.ibegin(*it); is != gv.iend(*it); ++is)
    {
      if (!is->boundary()) continue;

      // Get vertices on this boundary face
      int faceIdx = is->indexInInside();
      int nFaceVerts = ref.size(faceIdx, 1, dim);

      for (int k = 0; k < nFaceVerts; ++k)
      {
        int localVertInElem = ref.subEntity(faceIdx, 1, k, dim);
        int localIdx = set.subIndex(*it, localVertInElem, dim);
        PetscInt globalIdx = l2g[localIdx];

        // Get physical coordinates of this vertex
        auto vertexPos = it->geometry().corner(localVertInElem);
        PetscScalar gVal = (PetscScalar)bcFunc(vertexPos);

        bcGlobalRows.push_back(globalIdx);
        bcValues[globalIdx] = gVal;
      }
    }
  }

  // Remove duplicate global indices
  std::sort(bcGlobalRows.begin(), bcGlobalRows.end());
  bcGlobalRows.erase(std::unique(bcGlobalRows.begin(), bcGlobalRows.end()),
                     bcGlobalRows.end());

  // Zero the BC rows in A, putting 1.0 on the diagonal
  if (!bcGlobalRows.empty())
  {
    MatZeroRows(A, (PetscInt)bcGlobalRows.size(), bcGlobalRows.data(),
                1.0, NULL, NULL);
  }

  // Set b[i] = g(x_i) at boundary DOFs (use global indices)
  for (auto gi : bcGlobalRows)
  {
    PetscScalar val = 0.0;
    auto valIt = bcValues.find(gi);
    if (valIt != bcValues.end())
      val = valIt->second;
    VecSetValue(b, gi, val, INSERT_VALUES);
  }

  VecAssemblyBegin(b);
  VecAssemblyEnd(b);

  // Clean up the mapping
  ISLocalToGlobalMappingDestroy(&lgmap);

  if (rank == 0) {
    std::cout << "Inhomogeneous Dirichlet BCs applied." << std::endl;
    std::cout << "  z=0 face:  u = +1.0" << std::endl;
    std::cout << "  z=16 face: u = -1.0" << std::endl;
    std::cout << "  Other:     u =  0.0" << std::endl;
    std::cout << "Assembly complete." << std::endl;
  }
}


// ============================================================================
// Solve the linear system A u = b using PETSc KSP (CG solver)
// ============================================================================
template<class GV, class RHSFunc, class BCFunc>
void P1PoissonPetsc3D<GV, RHSFunc, BCFunc>::solve(int& argc, char**& argv)
{
  PetscMPIInt rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  KSP ksp;
  KSPCreate(PETSC_COMM_WORLD, &ksp);
  KSPSetOperators(ksp, A, A);

  // CG is appropriate for the symmetric positive-definite Poisson system
  KSPSetType(ksp, KSPCG);

  // Allow command-line override: -ksp_type gmres -pc_type asm etc.
  KSPSetFromOptions(ksp);
  KSPSetUp(ksp);

  if (rank == 0)
    std::cout << "KSP solver set up (CG). Solving 3D Poisson system..." << std::endl;

  KSPSolve(ksp, b, u);

  // Print convergence info
  PetscInt its;
  KSPConvergedReason reason;
  PetscReal rnorm;
  KSPGetIterationNumber(ksp, &its);
  KSPGetConvergedReason(ksp, &reason);
  KSPGetResidualNorm(ksp, &rnorm);

  if (rank == 0) {
    std::cout << "KSP converged in " << its << " iterations"
              << " (reason=" << reason << ", residual=" << rnorm << ")" << std::endl;
  }

  KSPDestroy(&ksp);
}


// ============================================================================
// MAIN
// ============================================================================
using namespace Dune;

int main(int argc, char** argv)
{
  try {
    // ── Initialize DUNE MPI + SLEPc (which also initializes PETSc) ──
    Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);
    SlepcInitialize(&argc, &argv, NULL, NULL);

    // ── Performance Timer ──
    QuDSimTimer timer(PETSC_COMM_WORLD);

    if (helper.rank() == 0) {
      std::cout << "======================================================" << std::endl;
      std::cout << "  3D Poisson FEM (P1) with DUNE + PETSc/SLEPc + MPI  " << std::endl;
      std::cout << "======================================================" << std::endl;
    }

    if (Dune::MPIHelper::isFake)
      std::cout << "This is a sequential program." << std::endl;
    else
      std::cout << "I am rank " << helper.rank() << " of " << helper.size()
                << " processes!" << std::endl;

    // ── Determine mesh file path ──
    // Default mesh path; can be overridden by first command-line argument
    std::string meshFile = "poisson3d.msh";
    if (argc > 1) {
      meshFile = argv[1];
    }

    // ── Grid setup: 3D ALUGrid (tetrahedral, conforming) ──
    timer.start("Mesh Loading");

    constexpr int dim = 3;
    using Grid = Dune::ALUGrid<dim, dim, Dune::simplex, Dune::conforming>;
    auto grid = Dune::GmshReader<Grid>::read("/home/athira/dune-projects/dune-p3d/src/poisson3d2.msh");
   // auto grid = Dune::GmshReader<Grid>::read(meshFile);
    grid->loadBalance();  // distribute mesh across MPI ranks

    auto gv = grid->leafGridView();

    typedef Grid                        GridType;
    typedef GridType::LeafGridView      GV;
    typedef GridType::ctype             ctype;
    typedef GaussianCharge<ctype, dim>  RHSFunc;
    typedef DirichletBC<ctype, dim>     BCFunc;

    if (helper.rank() == 0)
      std::cout << "3D Grid loaded and balanced from: " << meshFile << std::endl;
    timer.stop("Mesh Loading");

    std::cout << "Rank " << helper.rank() << ": "
              << gv.size(0) << " elements (tets), "
              << gv.size(dim) << " vertices (local)" << std::endl;

    // ── Create FEM solver ──
    RHSFunc rhs;
    BCFunc  bc;
    P1PoissonPetsc3D<GV, RHSFunc, BCFunc> fem(gv, rhs, bc);

    // ── Assemble ──
    timer.start("FEM Assembly");
    fem.assemble(argc, argv);
    timer.stop("FEM Assembly");

    // ── Solve ──
    timer.start("KSP Solve");
    fem.solve(argc, argv);
    timer.stop("KSP Solve");

    if (helper.rank() == 0)
      std::cout << "Solve done! Writing VTK output..." << std::endl;

    // ══════════════════════════════════════════════════════════════════════════
    // VTK output: scatter full solution to all ranks for parallel VTK writing
    // ══════════════════════════════════════════════════════════════════════════
    timer.start("VTK Output");

    Vec u_all = nullptr;
    VecScatter scatterAll = nullptr;

    VecScatterCreateToAll(fem.u, &scatterAll, &u_all);
    VecScatterBegin(scatterAll, fem.u, u_all, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(scatterAll, fem.u, u_all, INSERT_VALUES, SCATTER_FORWARD);

    // Each rank reads the full global solution and maps to its local vertices
    const PetscScalar* arr = nullptr;
    VecGetArrayRead(u_all, &arr);

    std::vector<double> u_vertex(gv.size(dim), 0.0);
    for (int i = 0; i < (int)gv.size(dim); ++i)
      u_vertex[i] = PetscRealPart(arr[fem.l2g[i]]);

    VecRestoreArrayRead(u_all, &arr);

    // Also compute and output the charge density at each vertex (for visualization)
    std::vector<double> rho_vertex(gv.size(dim), 0.0);
    for (auto it = gv.template begin<dim>(); it != gv.template end<dim>(); ++it)
    {
      int idx = gv.indexSet().index(*it);
      auto pos = it->geometry().center();
      rho_vertex[idx] = rhs(pos);
    }

    // Write parallel VTK
    Dune::VTKWriter<GV> vtkwriter(gv);
    vtkwriter.addVertexData(u_vertex, "potential_u");
    vtkwriter.addVertexData(rho_vertex, "charge_density_rho");
    vtkwriter.write("poisson3d_result", Dune::VTK::appendedraw);

    if (helper.rank() == 0)
      std::cout << "VTK file written: poisson3d_result" << std::endl;

    VecScatterDestroy(&scatterAll);
    VecDestroy(&u_all);
    timer.stop("VTK Output");

    // ── Compute and print solution statistics ──
    // NOTE: VecMin/VecMax are collective — ALL ranks must call them
    {
      PetscReal uMin, uMax;
      VecMin(fem.u, NULL, &uMin);
      VecMax(fem.u, NULL, &uMax);
      if (helper.rank() == 0) {
        std::cout << "\nSolution statistics:" << std::endl;
        std::cout << "  u_min = " << uMin << std::endl;
        std::cout << "  u_max = " << uMax << std::endl;
      }
    }

    // ── Print performance report ──
    timer.report();
    timer.reportCSV("timings.csv");
    timer.appendScalingLine("scaling_results.csv");

    // ── Cleanup ──
    fem.destroyPetsc();
    SlepcFinalize();

    return 0;
  }
  catch (Dune::Exception& e) {
    std::cerr << "DUNE error: " << e.what() << "\n";
  }
  catch (std::exception& e) {
    std::cerr << "std::exception: " << e.what() << "\n";
  }
  catch (...) {
    std::cerr << "Unknown exception!\n";
  }
  return 1;
}
