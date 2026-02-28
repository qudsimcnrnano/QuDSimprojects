#ifndef POSTPROCESSING_HH
#define POSTPROCESSING_HH

#include <dune/grid/common/mcmgmapper.hh>
#include <dune/geometry/referenceelements.hh>

// template<class GV, class Vector>
// class PostProcessing
// {
//   static const int dim = GV::dimension;
//   typedef typename GV::ctype ctype;
//   typedef typename GV::template Codim<dim>::Iterator VertexIterator;
//   typedef typename GV::template Codim<0>::Iterator LeafIterator;
//   typedef typename GV::IntersectionIterator IntersectionIterator;
//   //typedef typename Dune::template SingleCodimSingleGeomTypeMapper<GV,dim> VertexMap;
//   typedef typename GV::IndexSet LeafIndexSet;
//   typedef Dune::BCRSMatrix<Dune::FieldMatrix<ctype,1,1> > Matrix;
//   typedef Dune::MultipleCodimMultipleGeomTypeMapper<GV> VertexMap;
//   typedef Dune::MultipleCodimMultipleGeomTypeMapper<GV> ElementMap;
//   typedef Dune::BlockVector<Dune::FieldVector<ctype,1> > ScalarField;

// private:

//   //Vector Q1;
//   //Vector grad;

//   Vector grad;
//   std::vector < std::set<int> > adjacencyPattern;

// public:
//   const GV& gv;
//   const LeafIndexSet& set;
//   PostProcessing(const GV& gv);
//   void eval_surface_xgrad(Vector& Q1, double& xgrad);

// };

// template<class GV,class Vector>
// PostProcessing<GV,Vector>::PostProcessing(const GV& gv_):gv(gv_),set(gv.indexSet())
// {
//   grad.resize(gv.size(dim));

// }

// template<class GV,class Vector>
// void PostProcessing<GV,Vector>::eval_surface_xgrad(Vector& Q1, double& xgrad)
// {

//   int N=gv.size(dim);
//   ElementMap elementmap(gv, Dune::mcmgElementLayout());
//   VertexMap vertexmap(gv, Dune::mcmgVertexLayout());
//   VertexIterator vtend=gv.template end<dim>();
//   LeafIterator itend = gv.template end<0>();
// //   GradCalculator<GV,Vector> gc(gv,Q1,grad);
// //   gc.calculate();


//   double cum_grad=0;
//   double pts=0;
//   for(VertexIterator vt=gv.template begin<dim>();vt!=vtend;++vt)
//     if(vt->geometry().corner(0)[0]==0)
//       {
// 	cum_grad+=grad[vertexmap.index(*vt)];
// 	pts++;
//       }
  
//   cum_grad=sqrt(cum_grad)/pts;
//   cout<<"cum_grad = "<<cum_grad<<std::endl;
  

// //   std::ofstream grad_plot("grad_plot.dat");
// //   for (LeafIterator it = gv.template begin<0>(); it != itend; ++it) /*@\label{fem:loop1}@*/
// //     {
// //       grad_plot<<grad[elementmap.index(*it)]<<"        " <<it->geometry().center()[0]<<"        "<<it->geometry().center()[1]<<std::endl;
// //     }
//   std::ofstream grad_vtplot("grad_vtplot.dat");
//   for(VertexIterator vt=gv.template begin<dim>();vt!=vtend;++vt)
//     grad_vtplot<<grad[vertexmap.index(*vt)]<<"        " <<vt->geometry().corner(0)[0]<<"        "<<vt->geometry().corner(0)[1]<<std::endl;
  
// }



template<class GV, class Vector>
class GradCalculator
{

  static const int dim = GV::dimension;
  typedef typename GV::ctype ctype;
  typedef typename GV::template Codim<dim>::Iterator VertexIterator;
  typedef typename GV::template Codim<0>::Iterator LeafIterator;
  typedef typename GV::IntersectionIterator IntersectionIterator;
  //typedef typename Dune::template SingleCodimSingleGeomTypeMapper<GV,dim> VertexMap;
  typedef typename GV::IndexSet LeafIndexSet;
  typedef Dune::BCRSMatrix<Dune::FieldMatrix<ctype,1,1> > Matrix;
  typedef Dune::MultipleCodimMultipleGeomTypeMapper<GV> VertexMap;
  typedef Dune::MultipleCodimMultipleGeomTypeMapper<GV> ElementMap;
  typedef Dune::BlockVector<Dune::FieldVector<ctype,1> > ScalarField;

private:


  std::vector < std::set<int> > adjacencyPattern;


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
GradCalculator<GV, Vector>::GradCalculator(const GV& gv_, Vector& Q1_, Vector& gradele_):gv(gv_),Q1(Q1_),gradele(gradele_),set(gv.indexSet())
{


  determineAdjacencyPattern();

  int N=gv.size(dim);
  M.setSize(N, N, N + 2*gv.size(dim-1));
  M.setBuildMode(Matrix::random);
  for (int i = 0; i < N; i++)
    M.setrowsize(i,adjacencyPattern[i].size());
  M.endrowsizes();
  
  // set sparsity pattern of A with the information gained in determineAdjacencyPattern 
  for (int i = 0; i < N; i++) /*@\label{fem:setpattern}@*/
    {
      std::template set<int>::iterator setend = adjacencyPattern[i].end();
	    for (std::template set<int>::iterator setit = adjacencyPattern[i].begin();
		 setit != setend; ++setit)
	      M.addindex(i,*setit);
    }
  
  M.endindices(); /*@\label{fem:endindices}@*/
  
  ElementMap elementmap(gv, Dune::mcmgElementLayout());
  VertexMap vertexmap(gv, Dune::mcmgVertexLayout());
  c.resize(N);
  grad_ele.resize(elementmap.size());
  grad.resize(N);



}

template<class GV, class Vector>
void GradCalculator<GV, Vector>::determineAdjacencyPattern()
{
	const int N = gv.size(dim);
	adjacencyPattern.resize(N);
	std::cout<<"Number of unknowns: "<<N<<std::endl;
	
	const LeafIndexSet& set = gv.indexSet();
    const LeafIterator itend = gv.template end<0>();

    for (LeafIterator it = gv.template begin<0>(); it != itend; ++it)
    {
        Dune::GeometryType gt = it->type();
        const auto& ref =
            Dune::ReferenceElements<ctype,dim>::general(gt);

        // traverse all codim-1-entities of the current element and store all
        // pairs of vertices in adjacencyPattern
        const IntersectionIterator isend = gv.iend(*it);
        for (IntersectionIterator is = gv.ibegin(*it) ; is != isend ; ++is)
        {
            int vertexsize = ref.size(is->indexInInside(),1,dim);
            for (int i=0; i < vertexsize; i++)
            {
                int indexi = set.subIndex(*it,ref.subEntity(is->indexInInside(),1,i,dim),dim);
                for (int j=0; j < vertexsize; j++)
                {
                    int indexj = set.subIndex(*it,ref.subEntity(is->indexInInside(),1,j,dim),dim);
                    adjacencyPattern[indexi].insert(indexj);
                }
            }
        }
    }
}


template<class GV, class Vector>
void GradCalculator<GV, Vector>::calculate()
{
  VertexMap vertexmap(gv, Dune::mcmgVertexLayout());
  VertexIterator vtend=gv.template end<dim>();
  eval_grad_ele();
  //smoothen();

}

template<class GV, class Vector>
void GradCalculator<GV, Vector>::eval_grad_ele()
{
  int N=gv.size(dim);
  ElementMap elementmap(gv, Dune::mcmgElementLayout());
  VertexMap vertexmap(gv, Dune::mcmgVertexLayout());
  grad_ele=0.0;
  const P1ShapeFunctionSet<ctype,ctype,dim>& basis = P1ShapeFunctionSet<ctype,ctype,dim>::instance();
  VertexIterator vtend=gv.template end<dim>();
  LeafIterator itend = gv.template end<0>();
    for (LeafIterator it = gv.template begin<0>(); it != itend; ++it) /*@\label{fem:loop1}@*/
      {
        Dune::GeometryType gt = it->type();
        const auto& ref = Dune::ReferenceElements<ctype,dim>::general(gt);
        int vertexsize = ref.size(dim);
//         const Dune::QuadratureRule<ctype,dim>& rule = Dune::QuadratureRules<ctype,dim>::rule(gt,1);
//         for (typename Dune::QuadratureRule<ctype,dim>::const_iterator r = rule.begin(); r != rule.end() ; ++r)
// 	  {
	    Dune::FieldVector<ctype,dim> ctr_pt;
	    ctr_pt=it->geometry().center();
	    Dune::FieldMatrix<ctype,dim,dim> jacInvTra =  it->geometry().jacobianInverseTransposed(ref.position(0,0));
// 	    ctype detjac = it->geometry().integrationElement(r->position());
	    for(int i=0;i<vertexsize;i++)
	      {
		Dune::FieldVector<ctype,dim> grad1;
		jacInvTra.mv(basis[i].evaluateGradient(ref.position(0,0)),grad1);
		gradele[elementmap.index(*it)]+=Q1[set.subIndex(*it,i,dim)]*grad1[0];
	      }
// 	  }

      }

  std::ofstream grad_plot("gradele_plot.dat");
  for (LeafIterator it = gv.template begin<0>(); it != itend; ++it) /*@\label{fem:loop1}@*/
    {
      grad_plot<<gradele[elementmap.index(*it)]<<"        " <<it->geometry().center()[0]<<"        "<<it->geometry().center()[1]<<std::endl;
    }
  grad_plot.close();

}

template<class GV, class Vector>
void GradCalculator<GV, Vector>::assemble()
{

  const P1ShapeFunctionSet<ctype,ctype,dim>& basis = P1ShapeFunctionSet<ctype,ctype,dim>::instance();
  ElementMap elementmap(gv, Dune::mcmgElementLayout());
  //assemble M & c

  for (LeafIterator it = gv.template begin<0>(); it != gv.template end<0>(); ++it) /*@\label{fem:loop1}@*/
    {
      // determine geometry type of the current element and get the matching reference element
      Dune::GeometryType gt = it->type();
      const auto& ref =
      Dune::ReferenceElements<ctype,dim>::general(gt);
      int vertexsize = ref.size(dim);
      
      // get a quadrature rule of order one for the given geometry type
      const Dune::QuadratureRule<ctype,dim>& rule = Dune::QuadratureRules<ctype,dim>::rule(gt,5);
      for (typename Dune::QuadratureRule<ctype,dim>::const_iterator r = rule.begin();
	   r != rule.end() ; ++r)
	{
	  // compute the jacobian inverse transposed to transform the gradients
	  Dune::FieldMatrix<ctype,dim,dim> jacInvTra =
	    it->geometry().jacobianInverseTransposed(r->position());
	  // cout<<"r position "<<r->position()<<std::endl;
	  // get the weight at the current quadrature point
	  ctype weight = r->weight();
	  
	  // compute Jacobian determinant for the transformation formula
	  ctype detjac = it->geometry().integrationElement(r->position());
	  for (int i = 0; i < vertexsize; i++)
	    {
	      ctype phii = basis[i].evaluateFunction(r->position());
	      for (int j = 0; j < vertexsize; j++)
		{
		  ctype phij = basis[j].evaluateFunction(r->position());
		  M[set.subIndex(*it,i,dim)][set.subIndex(*it,j,dim)] /*@\label{fem:calca}@*/
 		    += (phii*phij) * weight *detjac;
		  
		  //     cout<<"A value " << (grad1*grad2) * weight * detjac<<std::endl;
		}


	    }
	}
      
    }

  for (LeafIterator it = gv.template begin<0>(); it != gv.template end<0>(); ++it) /*@\label{fem:loop1}@*/
    {
      // determine geometry type of the current element and get the matching reference element
      Dune::GeometryType gt = it->type();
      const auto& ref =
      Dune::ReferenceElements<ctype,dim>::general(gt);
      int vertexsize = ref.size(dim);
      
      // get a quadrature rule of order one for the given geometry type
      const Dune::QuadratureRule<ctype,dim>& rule = Dune::QuadratureRules<ctype,dim>::rule(gt,5);
      for (typename Dune::QuadratureRule<ctype,dim>::const_iterator r = rule.begin();
	   r != rule.end() ; ++r)
	{
	  // compute the jacobian inverse transposed to transform the gradients
	  Dune::FieldMatrix<ctype,dim,dim> jacInvTra =
	    it->geometry().jacobianInverseTransposed(r->position());
	  // cout<<"r position "<<r->position()<<std::endl;
	  // get the weight at the current quadrature point
	  ctype weight = r->weight();
	  
	  // compute Jacobian determinant for the transformation formula
	  ctype detjac = it->geometry().integrationElement(r->position());
	  for (int i = 0; i < vertexsize; i++)
	    {
	      ctype phii = basis[i].evaluateFunction(r->position());
	      c[set.subIndex(*it,i,dim)] += phii * grad_ele[elementmap.index(*it)] * weight *detjac;
	    }
	}
      
    }


}

template<class GV, class Vector>
void GradCalculator<GV, Vector>::smoothen()
{

  assemble();
  
  // make linear operator from M
  Dune::MatrixAdapter<Matrix,ScalarField,ScalarField> op(M);
  

  // initialize preconditioner
//   Dune::SeqILUn<Matrix,ScalarField,ScalarField> ilu1(M, 1, 0.92);
  Dune::SeqJac<Matrix,ScalarField,ScalarField> ilu1(M, 1, 0.92);

  // the inverse operator
  // Dune::BiCGSTABSolver<ScalarField> bcgs(op, ilu1, 1e-15, 5000, 0);
  
  Dune::BiCGSTABSolver<ScalarField> bcgs(op, ilu1, 1e-35, 5000, 0);
  Dune::InverseOperatorResult r;

  // initialue u to some arbitrary value to avoid u being the exact
  // solution
    
  grad.resize(c.N());
  grad = 1.0;

  // finally solve the system
  bcgs.apply(grad, c, r);

//   VertexMap vertexmap(gv, Dune::mcmgVertexLayout());
//   VertexIterator vtend=gv.template end<dim>();
// 		std::ofstream grad_plot("grad_vtplot.dat");
// 		for(VertexIterator vt=gv.template begin<dim>();vt!=vtend;++vt)
// 		  grad_plot<<vt->geometry().corner(0)[0]<<"       "<<vt->geometry().corner(0)[1]<<"        "<<grad[vertexmap.index(*vt)]<<std::endl;
// 		grad_plot.close();

}


template<class GV, class Vector>
class VectorIntegrator
{

  static const int dim = GV::dimension;
  typedef typename GV::ctype ctype;
  typedef typename GV::template Codim<dim>::Iterator VertexIterator;
  typedef typename GV::template Codim<0>::Iterator LeafIterator;
  typedef typename GV::IntersectionIterator IntersectionIterator;
  //typedef typename Dune::template SingleCodimSingleGeomTypeMapper<GV,dim> VertexMap;
  typedef typename GV::IndexSet LeafIndexSet;
  typedef Dune::BCRSMatrix<Dune::FieldMatrix<ctype,1,1> > Matrix;
  typedef Dune::MultipleCodimMultipleGeomTypeMapper<GV> VertexMap;
  typedef Dune::MultipleCodimMultipleGeomTypeMapper<GV> ElementMap;
  typedef Dune::BlockVector<Dune::FieldVector<ctype,1> > ScalarField;

private:

//   Vector& Q1;

  const GV& gv;
  const LeafIndexSet& set;

public:

  VectorIntegrator(const GV& gv_);
  void integrate(Vector& Q1, double& integral);


};

template<class GV, class Vector>
VectorIntegrator<GV, Vector>::VectorIntegrator(const GV& gv_):gv(gv_),set(gv.indexSet())
{
}

template<class GV, class Vector>
void VectorIntegrator<GV, Vector>::integrate(Vector& Q1, double& integral)
{

 double Q1_avg=0.0;
 integral=0.0;
 const LeafIterator itend = gv.template end<0>();
 for (LeafIterator it = gv.template begin<0>(); it != itend; ++it) /*@\label{fem:loop1}@*/
   {
     // determine geometry type of the current element and get the matching reference element
     Dune::GeometryType gt = it->type();
        const auto& ref = Dune::ReferenceElements<ctype,dim>::general(gt);
        int vertexsize = ref.size(dim);
   
        // get a quadrature rule of order one for the given geometry type
        const Dune::QuadratureRule<ctype,dim>& rule = Dune::QuadratureRules<ctype,dim>::rule(gt,2);
        for (typename Dune::QuadratureRule<ctype,dim>::const_iterator r = rule.begin(); r != rule.end() ; ++r)
        {
	  // get the weight at the current quadrature point
	  ctype weight = r->weight();
	  // compute Jacobian determinant for the transformation formula
	  ctype detjac = it->geometry().integrationElement(r->position());

	  Q1_avg=0.0;
	  for(int i=0;i<vertexsize;i++)
	    Q1_avg+=Q1[set.subIndex((*it),i,dim)]/vertexsize;
	  
	  integral+=Q1_avg*weight*detjac;

	  
	}
	
   }

}

// template<class GV, class Vector>
// VectorIntegrator<GV, Vector>::VectorIntegrator(GV& gv_:gv(gv_),set(gv.indexSet())
// 					       {}
		


#endif //POSTPROCESSING_HH
