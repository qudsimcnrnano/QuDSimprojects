#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
#include <sstream>
#include <iostream>
#include <algorithm>
#include <map>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/exceptions.hh>
#include <stdlib.h>
#include <cmath>
#include <set>
#include "petscmat.h"
#include "slepceps.h"
#include <config.h>
#include <vector>
#include <set>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include</opt/local/dune_210/dune-grid-2.10.dev20221009/dune/grid/albertagrid.hh>
#include </opt/local/dune_210/dune-alugrid-2.10.dev20221009/dune/alugrid/3d/alugrid.hh>

#include "shapefunctions.hh"
#include "qudsim_timer.hh"
#include <dune/istl/bvector.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/istl/io.hh>
#include </opt/local/dune_210/dune-grid-2.10.dev20221009/dune/grid/common/scsgmapper.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <stdlib.h>
#include <cmath>

#include</opt/local/dune_210/dune-grid-2.10.dev20221009/dune/grid/albertagrid.hh>
#include </opt/local/dune_210/dune-grid-howto-master/unitcube_alugrid.hh>
#include </opt/local/dune_210/dune-grid-2.10.dev20221009/dune/grid/io/file/gmshreader.hh>
#include </opt/local/dune_210/dune-alugrid-2.10.dev20221009/dune/alugrid/3d/gridfactory.hh>
#include <dune/grid/yaspgrid.hh>

#include <dune/grid/yaspgrid.hh>
#include </opt/local/dune_210/dune-grid-howto-master/unitcube_alugrid.hh>
#include </opt/local/dune_210/dune-grid-2.10.dev20221009/dune/grid/common/scsgmapper.hh>
#include </opt/local/dune_210/dune-grid-2.10.dev20221009/dune/grid/io/file/gmshreader.hh>
#include </opt/local/dune_210/dune-alugrid-2.10.dev20221009/dune/alugrid/3d/gridfactory.hh>

#if HAVE_DUNE_ISTL
#include <dune/istl/bvector.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/ilu.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/io.hh>
#else
#include <dune/common/dynvector.hh>
#include <dune/common/dynmatrix.hh>
#endif

#include </opt/local/dune_210/dune-grid-2.10.dev20221009/dune/grid/uggrid.hh>

#include "shapefunctions.hh"
#include "petsc.h"
#include "slepc.h"
#include "petscmat.h"
#include "slepceps.h"
#include "slepcsys.h"

#include <petscdm.h>
#include <petscdmda.h>
#include <petscksp.h>
#include <petscsnes.h>
#include <petscis.h>


// -------------------- PETSc-based P1 FEM assembler/solver --------------------
template<class GV, class F>
class P1ElementsPetsc
{
public:
  static const int dim = GV::dimension;
  typedef typename GV::ctype ctype;
  typedef Dune::BCRSMatrix<Dune::FieldMatrix<ctype,1,1> > Matrix;

private:
  using LeafIterator = typename GV::template Codim<0>::Iterator;
  using IntersectionIterator = typename GV::IntersectionIterator;
  using LeafIndexSet = typename GV::IndexSet;
  using JacobianInverseTransposed =
      typename GV::template Codim<0>::Geometry::JacobianInverseTransposed;

  const GV& gv;
  const F&  f;

public:
  Mat A;                    // PETSc parallel matrix
  Vec b;                    // PETSc parallel RHS vector
  Vec u;                    // PETSc parallel solution vector
  KSP ksp;
  PetscErrorCode ierr;

  std::vector<std::set<int>> adjacencyPattern;
  std::vector<PetscInt> l2g;   // local-to-global index mapping
  PetscInt globalN;            // total number of unique global vertices

  P1ElementsPetsc(const GV& gv_, const F& f_) : gv(gv_), f(f_), globalN(0),
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


// -------------------- Build global index mapping --------------------

template<class GV, class F>
void P1ElementsPetsc<GV,F>::buildGlobalIndexMap(int mpiRank, int mpiSize)
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
    displs[p] = displs[p-1] + allCounts[p-1];
  int totalDoubles = displs[mpiSize-1] + allCounts[mpiSize-1];

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
    std::cout << "Global index map built: " << globalN << " unique global vertices" << std::endl;
}


// -------------------- Adjacency pattern (local indices) --------------------
template<class GV, class F>
void P1ElementsPetsc<GV,F>::determineAdjacencyPattern()
{
  const int N = gv.size(dim);
  adjacencyPattern.clear();
  adjacencyPattern.resize(N);

  const LeafIndexSet& set = gv.indexSet();
  const LeafIterator itend = gv.template end<0>();

  for (LeafIterator it = gv.template begin<0>(); it != itend; ++it)
  {
    auto geo = it->geometry();
    auto ref = Dune::referenceElement(geo);

    const IntersectionIterator isend = gv.iend(*it);
    for (IntersectionIterator is = gv.ibegin(*it); is != isend; ++is)
    {
      int vertexsize = ref.size(is->indexInInside(), 1, dim);
      for (int i = 0; i < vertexsize; i++)
      {
        int indexi = set.subIndex(*it, ref.subEntity(is->indexInInside(),1,i,dim), dim);
        for (int j = 0; j < vertexsize; j++)
        {
          int indexj = set.subIndex(*it, ref.subEntity(is->indexInInside(),1,j,dim), dim);
          adjacencyPattern[indexi].insert(indexj);
        }
      }
    }
  }
}


// -------------------- Assemble --------------------
template<class GV, class F>
void P1ElementsPetsc<GV,F>::assemble(int& argc, char**& argv)
{
  const int N = gv.size(dim);  // local vertex count (including ghosts)
  const LeafIndexSet& set = gv.indexSet();
  const LeafIterator itend = gv.template end<0>();

  static char help[] = "poisson Solver";

  PetscFunctionBeginUser;
  ierr = PetscInitialize(&argc, &argv, (char*)0, help);

  PetscMPIInt rank, size;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  MPI_Comm_size(PETSC_COMM_WORLD, &size);

  std::cout << "Rank " << rank << ": N=" << N << " local vertices, "
            << gv.size(0) << " local elements" << std::endl;

  // ---- Build the global index mapping ----
  buildGlobalIndexMap(rank, size);

  // ---- Determine adjacency pattern (local) ----
  determineAdjacencyPattern();

  // ---- P1 basis ----
  P1ShapeFunctionSet<ctype,ctype,dim> basis = P1ShapeFunctionSet<ctype,ctype,dim>::instance();

 

  // ---- Create ISLocalToGlobalMapping from our l2g array ----
  ISLocalToGlobalMapping lgmap;
  ISLocalToGlobalMappingCreate(PETSC_COMM_WORLD, 1, N, l2g.data(),
                               PETSC_COPY_VALUES, &lgmap);

  // ---- Create PETSc Vec b and u using the local-to-global mapping ----
  VecCreate(PETSC_COMM_WORLD, &b);
  VecSetSizes(b, PETSC_DECIDE, globalN);
  VecSetFromOptions(b);
  VecSetLocalToGlobalMapping(b, lgmap);
  VecDuplicate(b, &u);
  VecSet(b, 0.0);

  if (rank == 0)
    std::cout << "PETSc vectors created (globalN=" << globalN << ")" << std::endl;

  // ---- Create PETSc Mat A with local-to-global mapping ----
  MatCreate(PETSC_COMM_WORLD, &A);
  MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, globalN, globalN);
  MatSetFromOptions(A);

  // Set the local-to-global mapping for both rows and columns
  MatSetLocalToGlobalMapping(A, lgmap, lgmap);

  // Preallocation: use a generous per-row estimate
  {
    int maxNnz = 0;
    for (int i = 0; i < N; ++i) {
      int sz = (int)adjacencyPattern[i].size();
      if (sz > maxNnz) maxNnz = sz;
    }
    // Use maxNnz for both diagonal and off-diagonal blocks
    MatMPIAIJSetPreallocation(A, maxNnz, NULL, maxNnz, NULL);
    MatSeqAIJSetPreallocation(A, maxNnz, NULL);
  }

  // Safety net: allow extra nonzeros without error
  MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);

  if (rank == 0)
    std::cout << "PETSc matrix created with local-to-global mapping." << std::endl;

  // ---- Element loop: assemble stiffness + RHS directly into PETSc ----
  for (LeafIterator it = gv.template begin<0>(); it != itend; ++it)
  {
    auto geo = it->geometry();
    auto ref = Dune::referenceElement(geo);
    int vertexsize = ref.size(dim);

    // Collect local vertex indices for this element
    std::vector<PetscInt> localRows(vertexsize);
    for (int i = 0; i < vertexsize; ++i)
      localRows[i] = (PetscInt)set.subIndex(*it, i, dim);

    // Compute element stiffness matrix
    std::vector<PetscScalar> Ke(vertexsize * vertexsize, 0.0);

    const auto& rule = Dune::QuadratureRules<ctype,dim>::rule(it->type(), 1);
    for (auto r = rule.begin(); r != rule.end(); ++r)
    {
      JacobianInverseTransposed JinvT =
        geo.jacobianInverseTransposed(r->position());

      ctype w    = r->weight();
      ctype detJ = geo.integrationElement(r->position());

      for (int i = 0; i < vertexsize; ++i)
      {
        Dune::FieldVector<ctype,dim> grad_i;
        JinvT.mv(basis[i].evaluateGradient(r->position()), grad_i);

        for (int j = 0; j < vertexsize; ++j)
        {
          Dune::FieldVector<ctype,dim> grad_j;
          JinvT.mv(basis[j].evaluateGradient(r->position()), grad_j);

          Ke[i * vertexsize + j] += (PetscScalar)((grad_i * grad_j) * w * detJ);
        }
      }
    }

   
    MatSetValuesLocal(A, vertexsize, localRows.data(),
                         vertexsize, localRows.data(),
                         Ke.data(), ADD_VALUES);

    // ---- RHS vector assembly for this element ----
    std::vector<PetscScalar> be(vertexsize, 0.0);
    const auto& rule2 = Dune::QuadratureRules<ctype,dim>::rule(it->type(), 2);
    for (auto r = rule2.begin(); r != rule2.end(); ++r)
    {
      ctype w    = r->weight();
      ctype detJ = geo.integrationElement(r->position());
      auto x     = geo.global(r->position());

      for (int i = 0; i < vertexsize; ++i)
      {
        ctype phi_i = basis[i].evaluateFunction(r->position());
        be[i] += (PetscScalar)(phi_i * f(x) * w * detJ);
      }
    }

    VecSetValuesLocal(b, vertexsize, localRows.data(), be.data(), ADD_VALUES);
  }

  if (rank == 0)
    std::cout << "Element loop complete. Assembling..." << std::endl;

  // ---- Assemble PETSc matrix and vector ----
  MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

  VecAssemblyBegin(b);
  VecAssemblyEnd(b);

  if (rank == 0)
    std::cout << "Matrix and vector assembled." << std::endl;

  // ---- Dirichlet BC (u = 0 on boundary) ----
  
  std::vector<PetscInt> bcRows;
  for (LeafIterator it = gv.template begin<0>(); it != itend; ++it)
  {
    auto ref = Dune::referenceElement(it->geometry());

    for (IntersectionIterator is = gv.ibegin(*it); is != gv.iend(*it); ++is)
    {
      if (!is->boundary()) continue;

      int nv = ref.size(is->indexInInside(), 1, dim);
      for (int k = 0; k < nv; ++k)
      {
        int lv = ref.subEntity(is->indexInInside(), 1, k, dim);
        int li = set.subIndex(*it, lv, dim);
        PetscInt gi = l2g[li];  // global index for MatZeroRows
        bcRows.push_back(gi);
      }
    }
  }

  // Remove duplicates
  std::sort(bcRows.begin(), bcRows.end());
  bcRows.erase(std::unique(bcRows.begin(), bcRows.end()), bcRows.end());

  // Zero the BC rows in A, putting 1.0 on the diagonal
  if (!bcRows.empty())
  {
    MatZeroRows(A, (PetscInt)bcRows.size(), bcRows.data(), 1.0, NULL, NULL);
  }

  // Set b = 0 at BC DOFs (use global indices with VecSetValue)
  for (auto gi : bcRows)
  {
    VecSetValue(b, gi, 0.0, INSERT_VALUES);
  }

  VecAssemblyBegin(b);
  VecAssemblyEnd(b);

  // Clean up the mapping
  ISLocalToGlobalMappingDestroy(&lgmap);

  if (rank == 0)
    std::cout << "Dirichlet BCs applied. Assembly complete." << std::endl;
}


// -------------------- Solve --------------------
template<class GV, class F>
void P1ElementsPetsc<GV,F>::solve(int& argc, char**& argv)
{
  PetscMPIInt rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  KSP ksp;
  KSPCreate(PETSC_COMM_WORLD, &ksp);
  KSPSetOperators(ksp, A, A);
  KSPSetType(ksp, KSPCG);
  KSPSetFromOptions(ksp);
  KSPSetUp(ksp);

  if (rank == 0)
    std::cout << "KSP solver set up. Solving..." << std::endl;

  KSPSolve(ksp, b, u);

  // Print convergence info
  PetscInt its;
  KSPConvergedReason reason;
  KSPGetIterationNumber(ksp, &its);
  KSPGetConvergedReason(ksp, &reason);
  if (rank == 0)
    std::cout << "KSP converged in " << its << " iterations (reason=" << reason << ")" << std::endl;

  KSPDestroy(&ksp);
}


// --------------------------- RHS function ---------------------------
template<class ctype, int dim>
class Bump {
public:
  ctype operator()(Dune::FieldVector<ctype,dim> x) const
  {
    ctype result = 0;
    for (int i = 0; i < dim; i++)
      result += 2.0 * x[i] * (1.0 - x[i]);
    return result;
  }
};


// ========================= MAIN =========================
using namespace Dune;

int main(int argc, char** argv)
{
  try {
    // DUNE MPI helper (also initializes MPI)
    Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);
    SlepcInitialize(&argc, &argv, NULL, NULL);

    // ── Performance Timer ──
    QuDSimTimer timer(PETSC_COMM_WORLD);

    if (helper.rank() == 0)
      std::cout << "Poisson FEM (P1) with DUNE + PETSc + MPI" << std::endl;

    if (Dune::MPIHelper::isFake)
      std::cout << "This is a sequential program." << std::endl;
    else
      std::cout << "I am rank " << helper.rank() << " of " << helper.size()
                << " processes!" << std::endl;

    // ---- Grid setup ----
    timer.start("Mesh Loading");
    constexpr int dim = 2;
    using Grid = Dune::ALUGrid<dim, 2, Dune::simplex, Dune::conforming>;

    auto grid = Dune::GmshReader<Grid>::read("/home/athira/dune-projects/dune-p1/src/rectangle.msh");
    grid->loadBalance();

    auto gv = grid->leafGridView();

    typedef Dune::ALUGrid<dim, 2, Dune::simplex, Dune::conforming> GridType;
    typedef GridType::LeafGridView GV;
    typedef GridType::ctype ctype;
    typedef Bump<ctype,dim> Func;

    if (helper.rank() == 0)
      std::cout << "Grid loaded and balanced." << std::endl;
    timer.stop("Mesh Loading");

    std::cout << "Rank " << helper.rank() << ": "
              << gv.size(0) << " elements, "
              << gv.size(dim) << " vertices (local)" << std::endl;

    // ---- FEM assembly + solve ----
    Func f;
    P1ElementsPetsc<GV,Func> p1(gv, f);

    timer.start("FEM Assembly");
    p1.assemble(argc, argv);
    timer.stop("FEM Assembly");

    timer.start("KSP Solve");
    p1.solve(argc, argv);
    timer.stop("KSP Solve");

    if (helper.rank() == 0) std::cout << "Solve done! Writing VTK..." << std::endl;

    // ---- VTK output ----
    timer.start("VTK Output");
    
    Vec u_all = nullptr;
    VecScatter scatterAll = nullptr;

    VecScatterCreateToAll(p1.u, &scatterAll, &u_all);
    VecScatterBegin(scatterAll, p1.u, u_all, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(scatterAll, p1.u, u_all, INSERT_VALUES, SCATTER_FORWARD);

    // Each rank reads the full global solution and maps to its local vertices
    PetscInt totalN;
    VecGetSize(p1.u, &totalN);

    const PetscScalar* arr = nullptr;
    VecGetArrayRead(u_all, &arr);

    std::vector<double> u_vertex(gv.size(dim), 0.0);
    for (int i = 0; i < (int)gv.size(dim); ++i)
      u_vertex[i] = PetscRealPart(arr[p1.l2g[i]]);

    VecRestoreArrayRead(u_all, &arr);

    // VTKWriter::write is collective — all ranks must call it
    Dune::VTKWriter<GV> vtkwriter(gv);
    vtkwriter.addVertexData(u_vertex, "u");
    vtkwriter.write("poisson_petsc", Dune::VTK::appendedraw);

    if (helper.rank() == 0)
      std::cout << "VTK file written: poisson_petsc" << std::endl;

    VecScatterDestroy(&scatterAll);
    VecDestroy(&u_all);
    timer.stop("VTK Output");

    // ── Print performance report ──
    timer.report();
    timer.reportCSV("timings.csv");
    timer.appendScalingLine("scaling_results.csv");

    // Cleanup
    p1.destroyPetsc();
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
