

#ifndef POISSONFEM_PARALLEL_HH
#define POISSONFEM_PARALLEL_HH

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/geometry/referenceelements.hh>
//#include <dune/grid/common/quadraturerules.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/ilu.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/io.hh>
#include <mpi.h>
#include "shapefunctions.hh"

template<class GV, class Vector>
class PoissonFEM_Parallel
{
  private:
    static const int dim = GV::dimension;
    typedef typename GV::ctype ctype;
    typedef typename GV::template Codim<0>::Iterator LeafIterator;
    typedef typename GV::template Codim<dim>::Iterator VertexIterator;
    typedef typename GV::IntersectionIterator IntersectionIterator;
    typedef typename GV::IndexSet LeafIndexSet;
    typedef Dune::BCRSMatrix<Dune::FieldMatrix<ctype,1,1> > Matrix;
    typedef Dune::BlockVector<Dune::FieldVector<ctype,1> > ScalarField;

    const GV& gv;
    Vector& charge;
    Vector& potential;
    double& Vss;
    const LeafIndexSet& set;
    Vector& deltaRho;
    Vector& deltapotential;
    int it0, N;
    int rank, nprocs;

    std::vector<std::set<int> > adjacencyPattern;
    void determineAdjacencyPattern();
    void assemble();
    void apply_bc();
    void solve();
    void solve_NR();

    Matrix A, NR;
    ScalarField b, R, deltaRho_b, AV, vertex_share, cnt, vertex_xcoord;
    ScalarField u;

  public:
    PoissonFEM_Parallel(const GV& gv_, Vector& charge_, Vector& potential_,
                        double& Vss, Vector& deltaRho_, Vector& deltapotential_,
                        int& it0_);
    void apply();
};


// =============================================================================
// Constructor
// =============================================================================
template<class GV, class Vector>
PoissonFEM_Parallel<GV,Vector>::PoissonFEM_Parallel(
    const GV& gv_, Vector& charge_, Vector& potential_, double& Vss_,
    Vector& deltaRho_, Vector& deltapotential_, int& it0_)
  : gv(gv_), charge(charge_), potential(potential_), set(gv.indexSet()),
    Vss(Vss_), deltaRho(deltaRho_), deltapotential(deltapotential_), it0(it0_)
{
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  N = gv.size(dim);
  if (rank == 0)
    std::cout << "  [PoissonFEM_Parallel] dim = " << dim << ", N_local = " << N << std::endl;

  determineAdjacencyPattern();

  // Setup matrices
  auto setupMatrix = [&](Matrix& M) {
    M.setSize(N, N, N + 2 * gv.size(dim - 1));
    M.setBuildMode(Matrix::random);
    for (int i = 0; i < N; i++)
      M.setrowsize(i, adjacencyPattern[i].size());
    M.endrowsizes();
    for (int i = 0; i < N; i++) {
      for (auto j : adjacencyPattern[i])
        M.addindex(i, j);
    }
    M.endindices();
  };

  setupMatrix(A);
  setupMatrix(NR);

  b.resize(N);
  R.resize(N);
  deltaRho_b.resize(N);
  AV.resize(N);
  u.resize(N);
  vertex_share.resize(N);
  cnt.resize(N);
  vertex_xcoord.resize(N);

  A = 0.0;   NR = 0.0;
  b = 0.0;   deltaRho_b = 0.0;
  R = 0.0;   AV = 0.0;
  u = 2.0;
  vertex_share = 0.0;
  cnt = 0.0;  vertex_xcoord = 0.0;

  // vertex share information
  for (LeafIterator it = gv.template begin<0>(); it != gv.template end<0>(); ++it) {
    Dune::GeometryType gt = it->type();
    const auto& ref =
        Dune::ReferenceElements<ctype,dim>::general(gt);
    int vertexsize = ref.size(dim);
    for (int k = 0; k < vertexsize; k++) {
      vertex_share[set.subIndex(*it, k, dim)] += 1;
      vertex_xcoord[set.subIndex(*it, k, dim)] = it->geometry().corner(k)[0];
    }
  }

  if (rank == 0) {
    std::ofstream vertexshareinfo("vertexshareinfo.dat");
    for (int i = 0; i < N; i++)
      vertexshareinfo << vertex_share[i] << "           " << i << std::endl;
    vertexshareinfo.close();
  }

  // Assemble stiffness matrix A
  const P1ShapeFunctionSet<ctype,ctype,dim>& basis = P1ShapeFunctionSet<ctype,ctype,dim>::instance();
  for (LeafIterator it = gv.template begin<0>(); it != gv.template end<0>(); ++it) {
    Dune::GeometryType gt = it->type();
    const auto& ref =
        Dune::ReferenceElements<ctype,dim>::general(gt);
    int vertexsize = ref.size(dim);

    const Dune::QuadratureRule<ctype,dim>& rule =
        Dune::QuadratureRules<ctype,dim>::rule(gt, 5);
    for (typename Dune::QuadratureRule<ctype,dim>::const_iterator r = rule.begin();
         r != rule.end(); ++r) {
      Dune::FieldMatrix<ctype,dim,dim> jacInvTra =
          it->geometry().jacobianInverseTransposed(r->position());
      ctype weight = r->weight();
      ctype detjac = it->geometry().integrationElement(r->position());

      for (int i = 0; i < vertexsize; i++) {
        Dune::FieldVector<ctype,dim> grad1;
        jacInvTra.mv(basis[i].evaluateGradient(r->position()), grad1);
        for (int j = 0; j < vertexsize; j++) {
          Dune::FieldVector<ctype,dim> grad2;
          jacInvTra.mv(basis[j].evaluateGradient(r->position()), grad2);
          A[set.subIndex(*it, i, dim)][set.subIndex(*it, j, dim)]
              += (grad1 * grad2) * weight * detjac;
        }
      }
    }
  }
}


// =============================================================================
// Determine adjacency pattern (local partition)
// =============================================================================
template<class GV, class Vector>
void PoissonFEM_Parallel<GV,Vector>::determineAdjacencyPattern()
{
  adjacencyPattern.resize(N);
  if (rank == 0)
    std::cout << "  [PoissonFEM] Number of local unknowns: " << N << std::endl;

  const LeafIterator itend = gv.template end<0>();
  for (LeafIterator it = gv.template begin<0>(); it != itend; ++it) {
    Dune::GeometryType gt = it->type();
    const auto& ref =
        Dune::ReferenceElements<ctype,dim>::general(gt);

    const IntersectionIterator isend = gv.iend(*it);
    for (IntersectionIterator is = gv.ibegin(*it); is != isend; ++is) {
      int vertexsize = ref.size(is->indexInInside(), 1, dim);
      for (int i = 0; i < vertexsize; i++) {
        int indexi = set.subIndex(*it, ref.subEntity(is->indexInInside(), 1, i, dim), dim);
        for (int j = 0; j < vertexsize; j++) {
          int indexj = set.subIndex(*it, ref.subEntity(is->indexInInside(), 1, j, dim), dim);
          adjacencyPattern[indexi].insert(indexj);
        }
      }
    }
  }
}


// =============================================================================
// Assemble RHS and Newton-Raphson matrix
// =============================================================================
template<class GV, class Vector>
void PoissonFEM_Parallel<GV,Vector>::assemble()
{
  b = 0.0;
  NR = 0.0;
  deltaRho_b = 0.0;
  cnt = 0.0;

  const LeafIterator itend = gv.template end<0>();
  const P1ShapeFunctionSet<ctype,ctype,dim>& basis = P1ShapeFunctionSet<ctype,ctype,dim>::instance();

  for (LeafIterator it = gv.template begin<0>(); it != itend; ++it) {
    Dune::GeometryType gt = it->type();
    const auto& ref =
        Dune::ReferenceElements<ctype,dim>::general(gt);
    int vertexsize = ref.size(dim);

    const Dune::QuadratureRule<ctype,dim>& rule2 =
        Dune::QuadratureRules<ctype,dim>::rule(gt, 5);
    for (typename Dune::QuadratureRule<ctype,dim>::const_iterator r = rule2.begin();
         r != rule2.end(); ++r) {
      ctype weight = r->weight();
      ctype detjac = it->geometry().integrationElement(r->position());

      for (int i = 0; i < vertexsize; i++) {
        ctype fval = basis[i].evaluateFunction(r->position())
                   * charge[set.subIndex(*it, i, dim)];
        ctype fval1 = basis[i].evaluateFunction(r->position())
                    * deltaRho[set.subIndex(*it, i, dim)];

        b[set.subIndex(*it, i, dim)]        += fval  * weight * detjac;
        deltaRho_b[set.subIndex(*it, i, dim)] += fval1 * weight * detjac;
      }
    }
  }

  NR = A;

  // Modify NR for Newton-Raphson
  for (LeafIterator it = gv.template begin<0>(); it != itend; ++it) {
    Dune::GeometryType gt = it->type();
    const auto& ref =
        Dune::ReferenceElements<ctype,dim>::general(gt);
    int vertexsize = ref.size(dim);

    for (int i = 0; i < vertexsize; i++) {
      for (int j = 0; j < vertexsize; j++) {
        int si = set.subIndex(*it, i, dim);
        int sj = set.subIndex(*it, j, dim);
        if ((si == sj) && (NR[si][sj] == A[si][sj])) {
          NR[si][sj] -= deltaRho_b[si];
        }
      }
    }
  }
}


// =============================================================================
// Apply boundary conditions
// =============================================================================
template<class GV, class Vector>
void PoissonFEM_Parallel<GV,Vector>::apply_bc()
{
  LeafIterator itend = gv.template end<0>();

  for (LeafIterator it = gv.template begin<0>(); it != itend; ++it) {
    const IntersectionIterator isend = gv.iend(*it);
    for (IntersectionIterator is = gv.ibegin(*it); is != isend; ++is) {
      Dune::GeometryType gt = it->type();
      const auto& ref =
          Dune::ReferenceElements<ctype,dim>::general(gt);

      if (is->boundary()) {
        for (int i = 0; i < ref.size(is->indexInInside(), 1, dim); i++) {
          int indexi = set.subIndex(*it, ref.subEntity(is->indexInInside(), 1, i, dim), dim);
          Dune::FieldVector<ctype, dim> point = is->geometry().center();

          if (point[0] == 0.0 || point[0] == 1.0) {
            A[indexi]         = 0.0;
            A[indexi][indexi]  = 1.0;
            b[indexi]         = 0;
            R[indexi]         = 0.0;
            NR[indexi][indexi] = 1.0;
          }
        }
      }
    }
  }
}


// =============================================================================
// Solve (first iteration - direct solve)
// =============================================================================
template<class GV, class Vector>
void PoissonFEM_Parallel<GV,Vector>::solve()
{
  if (rank == 0)
    std::cout << "  [PoissonFEM] Solving (initial iteration)..." << std::endl;

  Dune::MatrixAdapter<Matrix, ScalarField, ScalarField> op(A);
  Dune::SeqILU<Matrix, ScalarField, ScalarField> ilu1(A, 1, 0.92);
  Dune::BiCGSTABSolver<ScalarField> bcgs(op, ilu1, 1e-35, 5000, (rank == 0) ? 1 : 0);
  Dune::InverseOperatorResult r;

  potential.resize(b.N());
  potential = 1.0;
  bcgs.apply(potential, b, r);

  if (rank == 0)
    std::cout << "  [PoissonFEM] Initial solve converged: " << r.converged
              << " in " << r.iterations << " iterations." << std::endl;
  it0 = 1;
}


// =============================================================================
// Solve NR (Newton-Raphson iterations)
// =============================================================================
template<class GV, class Vector>
void PoissonFEM_Parallel<GV,Vector>::solve_NR()
{
  if (rank == 0)
    std::cout << "  [PoissonFEM] Solving Newton-Raphson..." << std::endl;

  Dune::MatrixAdapter<Matrix, ScalarField, ScalarField> op_NR(NR);
  Dune::SeqILU<Matrix, ScalarField, ScalarField> ilu1_NR(NR, 1, 0.1);
  Dune::BiCGSTABSolver<ScalarField> bcgs_NR(op_NR, ilu1_NR, 1e-100, 10000, (rank == 0) ? 1 : 0);
  Dune::InverseOperatorResult r;

  deltapotential.resize(R.N());
  deltapotential = 1.0;
  bcgs_NR.apply(deltapotential, R, r);

  if (rank == 0)
    std::cout << "  [PoissonFEM] NR solve converged: " << r.converged
              << " in " << r.iterations << " iterations." << std::endl;
}


// =============================================================================
// Apply: Run the full Poisson solve
// =============================================================================
template<class GV, class Vector>
void PoissonFEM_Parallel<GV,Vector>::apply()
{
  assemble();
  apply_bc();
  A.mv(potential, AV);
  R = b;
  R -= AV;
  apply_bc();

  if (it0 == 0)
    solve();
  else
    solve_NR();
}

#endif // POISSONFEM_PARALLEL_HH
