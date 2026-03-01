// =============================================================================
// PostProcessing utilities — DUNE 2.10 parallel-compatible
// Changes from serial:
//   1. SingleCodimSingleGeomTypeMapper -> MultipleCodimMultipleGeomTypeMapper
//   2. GenericReferenceElements -> ReferenceElements
//   3. .map() -> .index()
//   4. Constructor args: mcmgElementLayout(), mcmgVertexLayout()
// =============================================================================

#ifndef POSTPROCESSING_HH
#define POSTPROCESSING_HH

#include <fstream>
#include <cmath>
#include <set>
#include <vector>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/ilu.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/io.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/quadraturerules.hh>

#include "shapefunctions.hh"

using namespace std;


// =============================================================================
// GradCalculator: computes element-wise gradient and optionally smoothens
// =============================================================================
template<class GV, class Vector>
class GradCalculator
{
  static const int dim = GV::dimension;
  typedef typename GV::ctype ctype;
  typedef typename GV::template Codim<dim>::Iterator VertexIterator;
  typedef typename GV::template Codim<0>::Iterator LeafIterator;
  typedef typename GV::IntersectionIterator IntersectionIterator;
  typedef typename GV::IndexSet LeafIndexSet;
  typedef Dune::BCRSMatrix<Dune::FieldMatrix<ctype,1,1> > Matrix;

  // DUNE 2.10: MultipleCodimMultipleGeomTypeMapper
  typedef Dune::MultipleCodimMultipleGeomTypeMapper<GV> VertexMap;
  typedef Dune::MultipleCodimMultipleGeomTypeMapper<GV> ElementMap;
  typedef Dune::BlockVector<Dune::FieldVector<ctype,1> > ScalarField;

private:
  std::vector<std::set<int> > adjacencyPattern;

  Vector& Q1;
  Vector& gradele;
  Vector grad;
  Matrix M;
  Vector c;
  Vector grad_ele;

  const GV& gv;
  const LeafIndexSet& set;

  void determineAdjacencyPattern();
  void eval_grad_ele();
  void assemble();
  void smoothen();

public:
  GradCalculator(const GV& gv_, Vector& Q1_, Vector& gradele_);
  void calculate();
};


template<class GV, class Vector>
GradCalculator<GV,Vector>::GradCalculator(const GV& gv_, Vector& Q1_, Vector& gradele_)
  : gv(gv_), Q1(Q1_), gradele(gradele_), set(gv.indexSet())
{
  determineAdjacencyPattern();

  int N = gv.size(dim);
  M.setSize(N, N, N + 2 * gv.size(dim - 1));
  M.setBuildMode(Matrix::random);
  for (int i = 0; i < N; i++)
    M.setrowsize(i, adjacencyPattern[i].size());
  M.endrowsizes();

  for (int i = 0; i < N; i++) {
    for (auto jj : adjacencyPattern[i])
      M.addindex(i, jj);
  }
  M.endindices();

  // DUNE 2.10: layout argument required
  ElementMap elementmap(gv, Dune::mcmgElementLayout());
  VertexMap vertexmap(gv, Dune::mcmgVertexLayout());

  c.resize(N, false);
  grad_ele.resize(elementmap.size());
  grad.resize(N, false);
}


template<class GV, class Vector>
void GradCalculator<GV,Vector>::determineAdjacencyPattern()
{
  const int N = gv.size(dim);
  adjacencyPattern.resize(N);
  std::cout << "Number of unknowns: " << N << std::endl;

  const LeafIterator itend = gv.template end<0>();
  for (LeafIterator it = gv.template begin<0>(); it != itend; ++it)
  {
    Dune::GeometryType gt = it->type();
    // DUNE 2.10: ReferenceElements replaces GenericReferenceElements
    const auto& ref = Dune::ReferenceElements<ctype,dim>::general(gt);

    const IntersectionIterator isend = gv.iend(*it);
    for (IntersectionIterator is = gv.ibegin(*it); is != isend; ++is)
    {
      int vertexsize = ref.size(is->indexInInside(), 1, dim);
      for (int i = 0; i < vertexsize; i++)
      {
        int indexi = set.subIndex(*it, ref.subEntity(is->indexInInside(), 1, i, dim), dim);
        for (int j = 0; j < vertexsize; j++)
        {
          int indexj = set.subIndex(*it, ref.subEntity(is->indexInInside(), 1, j, dim), dim);
          adjacencyPattern[indexi].insert(indexj);
        }
      }
    }
  }
}


template<class GV, class Vector>
void GradCalculator<GV,Vector>::calculate()
{
  VertexMap vertexmap(gv, Dune::mcmgVertexLayout());
  VertexIterator vtend = gv.template end<dim>();
  eval_grad_ele();
}


template<class GV, class Vector>
void GradCalculator<GV,Vector>::eval_grad_ele()
{
  int N = gv.size(dim);
  // DUNE 2.10: layout argument required
  ElementMap elementmap(gv, Dune::mcmgElementLayout());
  VertexMap vertexmap(gv, Dune::mcmgVertexLayout());

  grad_ele = 0.0;
  P1ShapeFunctionSet<ctype,ctype,dim> basis = P1ShapeFunctionSet<ctype,ctype,dim>::instance();

  VertexIterator vtend = gv.template end<dim>();
  LeafIterator itend = gv.template end<0>();

  for (LeafIterator it = gv.template begin<0>(); it != itend; ++it)
  {
    Dune::GeometryType gt = it->type();
    const auto& ref = Dune::ReferenceElements<ctype,dim>::general(gt);
    int vertexsize = ref.size(dim);

    Dune::FieldVector<ctype,dim> ctr_pt;
    ctr_pt = it->geometry().center();
    Dune::FieldMatrix<ctype,dim,dim> jacInvTra =
        it->geometry().jacobianInverseTransposed(ref.position(0, 0));

    for (int i = 0; i < vertexsize; i++)
    {
      Dune::FieldVector<ctype,dim> grad1;
      jacInvTra.mv(basis[i].evaluateGradient(ref.position(0, 0)), grad1);
      // DUNE 2.10: .index() replaces .map()
      gradele[elementmap.index(*it)] += Q1[set.subIndex(*it, i, dim)] * grad1[0];
    }
  }

  std::ofstream grad_plot("gradele_plot.dat");
  for (LeafIterator it = gv.template begin<0>(); it != itend; ++it)
  {
    grad_plot << gradele[elementmap.index(*it)] << "        "
              << it->geometry().center()[0] << "        "
              << it->geometry().center()[1] << endl;
  }
  grad_plot.close();
}


template<class GV, class Vector>
void GradCalculator<GV,Vector>::assemble()
{
  P1ShapeFunctionSet<ctype,ctype,dim> basis = P1ShapeFunctionSet<ctype,ctype,dim>::instance();
  ElementMap elementmap(gv, Dune::mcmgElementLayout());

  // Assemble mass matrix M
  for (LeafIterator it = gv.template begin<0>(); it != gv.template end<0>(); ++it)
  {
    Dune::GeometryType gt = it->type();
    const auto& ref = Dune::ReferenceElements<ctype,dim>::general(gt);
    int vertexsize = ref.size(dim);

    const Dune::QuadratureRule<ctype,dim>& rule =
        Dune::QuadratureRules<ctype,dim>::rule(gt, 5);
    for (typename Dune::QuadratureRule<ctype,dim>::const_iterator r = rule.begin();
         r != rule.end(); ++r)
    {
      Dune::FieldMatrix<ctype,dim,dim> jacInvTra =
          it->geometry().jacobianInverseTransposed(r->position());
      ctype weight = r->weight();
      ctype detjac = it->geometry().integrationElement(r->position());

      for (int i = 0; i < vertexsize; i++)
      {
        ctype phii = basis[i].evaluateFunction(r->position());
        for (int j = 0; j < vertexsize; j++)
        {
          ctype phij = basis[j].evaluateFunction(r->position());
          M[set.subIndex(*it, i, dim)][set.subIndex(*it, j, dim)]
              += (phii * phij) * weight * detjac;
        }
      }
    }
  }

  // Assemble RHS c
  for (LeafIterator it = gv.template begin<0>(); it != gv.template end<0>(); ++it)
  {
    Dune::GeometryType gt = it->type();
    const auto& ref = Dune::ReferenceElements<ctype,dim>::general(gt);
    int vertexsize = ref.size(dim);

    const Dune::QuadratureRule<ctype,dim>& rule =
        Dune::QuadratureRules<ctype,dim>::rule(gt, 5);
    for (typename Dune::QuadratureRule<ctype,dim>::const_iterator r = rule.begin();
         r != rule.end(); ++r)
    {
      Dune::FieldMatrix<ctype,dim,dim> jacInvTra =
          it->geometry().jacobianInverseTransposed(r->position());
      ctype weight = r->weight();
      ctype detjac = it->geometry().integrationElement(r->position());

      for (int i = 0; i < vertexsize; i++)
      {
        ctype phii = basis[i].evaluateFunction(r->position());
        c[set.subIndex(*it, i, dim)] +=
            phii * grad_ele[elementmap.index(*it)] * weight * detjac;
      }
    }
  }
}


template<class GV, class Vector>
void GradCalculator<GV,Vector>::smoothen()
{
  assemble();

  Dune::MatrixAdapter<Matrix, ScalarField, ScalarField> op(M);
  Dune::SeqJac<Matrix, ScalarField, ScalarField, 1> ilu1(M, 1, 0.92);
  Dune::BiCGSTABSolver<ScalarField> bcgs(op, ilu1, 1e-35, 5000, 0);
  Dune::InverseOperatorResult r;

  grad.resize(c.N(), false);
  grad = 1.0;
  bcgs.apply(grad, c, r);
}


// =============================================================================
// VectorIntegrator: integrates a vertex-based field over the mesh
// =============================================================================
template<class GV, class Vector>
class VectorIntegrator
{
  static const int dim = GV::dimension;
  typedef typename GV::ctype ctype;
  typedef typename GV::template Codim<dim>::Iterator VertexIterator;
  typedef typename GV::template Codim<0>::Iterator LeafIterator;
  typedef typename GV::IntersectionIterator IntersectionIterator;
  typedef typename GV::IndexSet LeafIndexSet;
  typedef Dune::BCRSMatrix<Dune::FieldMatrix<ctype,1,1> > Matrix;

  // DUNE 2.10: MultipleCodimMultipleGeomTypeMapper
  typedef Dune::MultipleCodimMultipleGeomTypeMapper<GV> VertexMap;
  typedef Dune::MultipleCodimMultipleGeomTypeMapper<GV> ElementMap;
  typedef Dune::BlockVector<Dune::FieldVector<ctype,1> > ScalarField;

private:
  const GV& gv;
  const LeafIndexSet& set;

public:
  VectorIntegrator(const GV& gv_);
  void integrate(Vector& Q1, double& integral);
};


template<class GV, class Vector>
VectorIntegrator<GV,Vector>::VectorIntegrator(const GV& gv_)
  : gv(gv_), set(gv.indexSet())
{
}


template<class GV, class Vector>
void VectorIntegrator<GV,Vector>::integrate(Vector& Q1, double& integral)
{
  double Q1_avg = 0.0;
  integral = 0.0;

  const LeafIterator itend = gv.template end<0>();
  for (LeafIterator it = gv.template begin<0>(); it != itend; ++it)
  {
    Dune::GeometryType gt = it->type();
    // DUNE 2.10: ReferenceElements
    const auto& ref = Dune::ReferenceElements<ctype,dim>::general(gt);
    int vertexsize = ref.size(dim);

    const Dune::QuadratureRule<ctype,dim>& rule =
        Dune::QuadratureRules<ctype,dim>::rule(gt, 2);
    for (typename Dune::QuadratureRule<ctype,dim>::const_iterator r = rule.begin();
         r != rule.end(); ++r)
    {
      ctype weight = r->weight();
      ctype detjac = it->geometry().integrationElement(r->position());

      Q1_avg = 0.0;
      for (int i = 0; i < vertexsize; i++)
        Q1_avg += Q1[set.subIndex((*it), i, dim)] / vertexsize;

      integral += Q1_avg * weight * detjac;
    }
  }
}


#endif // POSTPROCESSING_HH
