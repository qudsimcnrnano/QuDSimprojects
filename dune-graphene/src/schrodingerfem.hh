

#ifndef SCHRODINGERFEM_HH
#define SCHRODINGERFEM_HH

#include <sstream>
#include <iostream>
#include "petscmat.h"
#include "slepceps.h"
#include <dune/common/fvector.hh>
#include <dune/geometry/quadraturerules.hh>
#include "shapefunctions.hh"
#include "pk2dlocalbasis.hh"
#include <dune/istl/bvector.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/istl/io.hh>
#include <dune/grid/common/mcmgmapper.hh>
//#include <dune/geometry/sgrid.hh>
#include"wavefunction.hh"



// #include <sstream>
// #include "petscmat.h"
// #include "slepceps.h"
// #include <dune/common/fvector.hh>
// #include <dune/grid/common/quadraturerules.hh>
// #include "shapefunctions.hh"
// #include <dune/istl/bvector.hh>
// #include <dune/istl/bcrsmatrix.hh>
// #include <dune/common/fvector.hh>
// #include <dune/common/fmatrix.hh>
// #include <dune/istl/io.hh>
// #include <dune/istl/bvector.hh>
// #include <dune/istl/bcrsmatrix.hh>
// #include <dune/istl/ilu.hh>
// #include <dune/istl/operators.hh>
// #include <dune/istl/solvers.hh>
// #include <dune/istl/preconditioners.hh>
// #include <dune/istl/io.hh>
// #include <dune/istl/io.hh>
// #include <petscksp.h>





template<class GV, class Vector>
class SchrodingerFEM
{
	private:
		static const int dim = GV::dimension;
		static const int dim1 = 2;   // for surface integral
		typedef typename GV::ctype ctype;
                typedef typename GV::template Codim<dim>::Iterator VertexIterator;
		typedef typename GV::template Codim<0>::Iterator LeafIterator;
              	typedef typename GV::template Codim<1>::Iterator triangleIterator;       // for surface integral
	   	typedef typename GV::IntersectionIterator IntersectionIterator;
  		//typedef typename Dune::template SingleCodimSingleGeomTypeMapper<GV,dim> VertexMap;
        	typedef typename GV::IndexSet LeafIndexSet;
                typedef typename GV::IndexSet TriangleIndexSet;
		typedef Dune::BCRSMatrix<Dune::FieldMatrix<ctype,1,1> > Matrix;
		
		typedef Dune::SingleCodimSingleGeomTypeMapper<GV,dim> VertexMap;

	

		
		const GV& gv;
		Vector& charge;
		Vector charge_tmp;
		Vector wavevector;
		Vector& potential;
		Vector& wavefunction_bc;
                std::vector<int> periodicmap;
  //		Vector& bc;
		const LeafIndexSet& set;
		//const TriangleIndexSet& tset;
		PetscInt nnz[];
		std::vector < std::set<int> > adjacencyPattern;
		void SchrodingerFEM<GV, Vector>::determineAdjacencyPattern();
		void assemble(double m);
		void apply_bc();
                void solve(double m, double g);
		 


		Matrix A;
		Matrix B;
		//Matrix Cb;
		Matrix C1;
                Matrix C2;
		//int twf[21011];
		Mat         	 A_p,B_p;		  
  		EPS         	 eps;		  
  		const EPSType    type;
                PetscReal   	 error, tol, re, im,mod_xr;
                PetscScalar 	 kr, ki,xr_val;
  		PetscErrorCode   ierr;
  		PetscInt    	 nev, maxit, i, its, lits, nconv;
                Vec               xr,xi;
                float             tmp;
                PetscInt          N_xr;
                int carrier_type;
                int car_type;
                double Ef;
		
	        public:
  SchrodingerFEM(const GV& gv_, Vector& charge_, Vector& potential_, Vector& wavefunction_bc, int carrier_type_, double Ef_, int& argc, char**& argv);
		void apply();
		
};




template<class GV, class Vector>
SchrodingerFEM<GV,Vector>::SchrodingerFEM(const GV& gv_, Vector& charge_, Vector& potential_, Vector& wavefunction_bc_, int carrier_type_, double Ef_, int& argc, char**& argv) : gv(gv_), charge(charge_), potential(potential_), wavefunction_bc(wavefunction_bc_), set(gv.indexSet())
{
  carrier_type=carrier_type_;
  Ef=Ef_;
	const int N = gv.size(dim);
	PetscInt nnz[N];


	
	
	determineAdjacencyPattern();
	
	for(int i=0; i<N ; ++i)
 	 nnz[i]=adjacencyPattern[i].size();
         
	
	char** null_argv;
	int null_argc= 0;
	static char help[] = "Schrodinger Solver";
	SlepcInitialize(&argc,&argv,(char*)0,help);

	MatCreate(MPI_COMM_WORLD, &A_p);
	MatCreate(MPI_COMM_WORLD, &B_p);
  	MatSetType(A_p, MATSEQAIJ);	
 	MatSetType(B_p, MATSEQAIJ);
	//MatSetType(A_p, MATMPIBAIJ);	
	//	MatSetType(B_p, MATMPIBAIJ);
	MatSetSizes(A_p, N, N, N, N);
	MatSetSizes(B_p, N, N, N, N);
	
	MatSeqAIJSetPreallocation(A_p,0,nnz);
	MatSeqAIJSetPreallocation(B_p,0,nnz);
	
	A.setSize(N, N, N + 2*gv.size(dim-1));
	A.setBuildMode(Matrix::random);

    	B.setSize(N, N, N + 2*gv.size(dim-1));
    	B.setBuildMode(Matrix::random);

	//Cb.setSize(N, N, N + 2*gv.size(dim-1));
    	//Cb.setBuildMode(Matrix::random);

	C1.setSize(N, N, N + 2*gv.size(dim-1));
    	C1.setBuildMode(Matrix::random);
        
	C2.setSize(N, N, N + 2*gv.size(dim-1));
    	C2.setBuildMode(Matrix::random);  
    for (int i = 0; i < N; i++)
	A.setrowsize(i,adjacencyPattern[i].size());
    	A.endrowsizes();

    for (int i = 0; i < N; i++)
    	B.setrowsize(i,adjacencyPattern[i].size());
   	B.endrowsizes();

   // for (int i = 0; i < N; i++)
    //	Cb.setrowsize(i,adjacencyPattern[i].size());
   //	Cb.endrowsizes();

    for (int i = 0; i < N; i++)
    	C1.setrowsize(i,adjacencyPattern[i].size());
   	C1.endrowsizes();

  for (int i = 0; i < N; i++)
    	C2.setrowsize(i,adjacencyPattern[i].size());
   	C2.endrowsizes();

    


    // set sparsity pattern of A with the information gained in determineAdjacencyPattern 
    for (int i = 0; i < N; i++) 
    {
        std::template set<int>::iterator setend = adjacencyPattern[i].end();
        for (std::template set<int>::iterator setit = adjacencyPattern[i].begin();
             setit != setend; ++setit)
            A.addindex(i,*setit);
    }
    A.endindices(); 




    // set sparsity pattern of B with the information gained in determineAdjacencyPattern 
    for (int i = 0; i < N; i++) 
    {
        std::template set<int>::iterator setend = adjacencyPattern[i].end();
        for (std::template set<int>::iterator setit = adjacencyPattern[i].begin();
             setit != setend; ++setit)
            B.addindex(i,*setit);
    }
    B.endindices(); 

    // set sparsity pattern of Cb with the information gained in determineAdjacencyPattern 
    //for (int i = 0; i < N; i++) 
   // {
     //   std::template set<int>::iterator setend = adjacencyPattern[i].end();
       // for (std::template set<int>::iterator setit = adjacencyPattern[i].begin();
   //          setit != setend; ++setit)
   //         Cb.addindex(i,*setit);
   // }
   // Cb.endindices(); 

    // set sparsity pattern of C with the information gained in determineAdjacencyPattern 
    for (int i = 0; i < N; i++) 
    {
        std::template set<int>::iterator setend = adjacencyPattern[i].end();
        for (std::template set<int>::iterator setit = adjacencyPattern[i].begin();
             setit != setend; ++setit)
            C1.addindex(i,*setit);
    }
    C1.endindices(); 


for (int i = 0; i < N; i++) 
    {
        std::template set<int>::iterator setend = adjacencyPattern[i].end();
        for (std::template set<int>::iterator setit = adjacencyPattern[i].begin();
             setit != setend; ++setit)
            C2.addindex(i,*setit);
    }
    C2.endindices(); 





    charge_tmp.resize(N, false);
    wavevector.resize(N, false);

    // initialize A and b
    A = 0.0;
    B = 0.0;
   // Cb = 0.0;
    C1 = 0.0;
    C2 = 0.0;
    charge_tmp=0.0;
    wavevector=0.0;



}


template<class GV, class Vector>
void SchrodingerFEM<GV, Vector>::determineAdjacencyPattern()
{
	const int N = gv.size(dim);
	adjacencyPattern.resize(N);
	std::cout<<"Number of unknowns: "<<N<<std::endl;
	
	const LeafIndexSet& set = gv.indexSet();
        //const TriangleIndexSet& set = gv.indexSet();
    const LeafIterator itend = gv.template end<0>();

    for (LeafIterator it = gv.template begin<0>(); it != itend; ++it)
    {
        Dune::GeometryType gt = it->type();
        const Dune::template ReferenceElements<ctype,dim> &ref =
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
void SchrodingerFEM<GV,Vector>::assemble(double m)
{


VertexMap vertexmap(gv);
Vector meff(vertexmap.size());



 		//VertexIterator vtend=gv.end<dim>();
		VertexIterator vtend=gv.template end<dim>();

		//for(VertexIterator vt=gv.begin<dim>();vt!=vtend;++vt)


             for(VertexIterator vt=gv.template begin<dim>();vt!=vtend;++vt)
		{
		double x_cord1 = vt->geometry().corner(0)[0];
		double y_cord1 = vt->geometry().corner(0)[1];
		//double z_cord1 = vt->geometry().corner(0)[2];
		double rho = sqrt(pow(x_cord1,2)+pow(y_cord1,2));
		double R1 = 60;
                //double R21 = 11.9929/0.052;
 		// Boundary belongs to core

			//if(rho<=R1)
			//{
                        //cout<<"I am in Dot  "<< rho<< endl;
			meff= 1;

			//}
			//else
		//	{
		//	meff[vertexmap.map(*vt)] = 1;
			// cout<<"I am in CdS  "<< rho<< endl;
		//	}
		}





    A = 0.0;
    B = 0.0;
   // Cb = 0.0;
    C1 = 0.0;
    C2 = 0.0;
    
  	
	



  // Volume integral //
double x_coordinate1=0.0;
double x_coordinate2=0.0;
double x_cord1=0.0;
double x_cord2=0.0;
double y_coordinate=0.0;
double y_cord=0.0;
double r_cord=0.0;
  const LeafIterator itend = gv.template end<0>();

	Dune::Pk2DLocalBasis<ctype,ctype,dim> basisL = Dune::Pk2DLocalBasis<ctype,ctype,dim>::instance();

    P1ShapeFunctionSet<ctype,ctype,dim> basis = P1ShapeFunctionSet<ctype,ctype,dim>::instance();

    for (LeafIterator it = gv.template begin<0>(); it != itend; ++it) //@\label{fem:loop1}@
      {

        // determine geometry type of the current element and get the matching reference element
        Dune::GeometryType gt = it->type();
        const Dune::template ReferenceElements<ctype,dim> &ref = Dune::ReferenceElements<ctype,dim>::general(gt);
        int vertexsize = ref.size(dim);
   
        // get a quadrature rule of order one for the given geometry type
        const Dune::QuadratureRule<ctype,dim>& rule = Dune::QuadratureRules<ctype,dim>::rule(gt,2);

        for (typename Dune::QuadratureRule<ctype,dim>::const_iterator r = rule.begin(); r != rule.end() ; ++r)
	  {
            // compute the jacobian inverse transposed to transform the gradients
            Dune::FieldMatrix<ctype,dim,dim> jacInvTra =  it->geometry().jacobianInverseTransposed(r->position());
	    
            // get the weight at the current quadrature point
            ctype weight = r->weight();
	    
            // compute Jacobian determinant for the transformation formula
            ctype detjac = it->geometry().integrationElement(r->position());

            for (int i = 0; i < vertexsize; i++)
	      {

			x_coordinate1 = it->geometry().corner(i)[0];
                        //y_coordinate = it->geometry().corner(i)[1];

                // compute transformed gradients
                Dune::FieldVector<ctype,dim> grad1;
                jacInvTra.mv(basisL[i].evaluateGradient(r->position()),grad1);	      // Try replacing the center method with a corner method
    //         if(x_coordinate<0)
      //        grad1[0]=grad1[0]*(-1);

                ctype phii = basisL[i].evaluateFunction(r->position()) ;
     //cout<<basis[i]<<endl;
		//cout<<basis[i].evaluateFunction(it->geometry().corner(i)); 

                for (int j = 0; j < vertexsize; j++)
		  {
                    Dune::FieldVector<ctype,dim> grad2;
                    jacInvTra.mv(basisL[j].evaluateGradient(r->position()),grad2);		// Try replacing the center method with a corner method
			
			
		    //grad2*=0.5;
		    ctype phij = basisL[j].evaluateFunction(r->position()) ;
			
			//x_coordinate = it->geometry().corner(i)[0];
                        y_coordinate = it->geometry().corner(i)[1];
			x_coordinate2 = it->geometry().corner(j)[0];


x_cord1=x_coordinate1;
x_cord2=x_coordinate2;
y_cord=y_coordinate;		//cout<<x_cord<<endl;



		




//WORKING CODE FOR

if(fabs(x_cord1)==0)
{
//cout<<x_cord<<endl;
A[set.subIndex(*it,i,dim)][set.subIndex(*it,j,dim)] += (grad1[1]*grad2[1])*weight*detjac +(2*phii* (potential[set.subIndex(*it,i,dim)])*phij) * weight * detjac;

else
A[set.subIndex(*it,i,dim)][set.subIndex(*it,j,dim)] += (grad1*grad2)*weight*detjac +(2*phii* (potential[set.subIndex(*it,i,dim)])*phij) * weight * detjac+((-grad1[0]*phij)/x_cord)*weight*detjac;
B[set.subIndex(*it,i,dim)][set.subIndex(*it,j,dim)] += 2*(phii*phij) * weight * detjac;

}
//END OF WORKING PART



		  }
	      }
	  }
      }	
    
}

template<class GV, class Vector>
void SchrodingerFEM<GV,Vector>::apply_bc()
{
 C1 = 0.0;
C2 = 0.0;
double pi=3.142;
                          double R2=120;
double Xr=125;
double Yr=76;
double Zr=76;
double XrL=-356;
double x_cord,y_cord;
int out=0;
int out1=0;
double rr_check=0.0;
double rr_check1=0.0;
  const int N = gv.size(dim);
  


  LeafIterator itend  = gv.template end<0>();



for ( LeafIterator it = gv.template begin<0>() ; it != itend ; ++it) 
    {
      const IntersectionIterator isend = gv.iend(*it);
      for (IntersectionIterator is = gv.ibegin(*it) ; is != isend ; ++is)
        {
	  // determine geometry type of the current element and get the matching reference element
	  Dune::GeometryType gt = it->type();
	  const Dune::template ReferenceElements<ctype,dim> &ref =
	  Dune::ReferenceElements<ctype,dim>::general(gt);
	  int vertexsize = ref.size(dim);
	  // check whether current intersection is on the boundary
	  if ( is->boundary() )
            {
	      Dune::FieldVector< ctype, dim > point = is -> geometry().center();
	      for (int i=0; i < ref.size(is->indexInInside(),1,dim); i++)

                {
		  // and replace the associated line of A and b with a trivial one
		  int indexi = set.subIndex(*it,ref.subEntity(is->indexInInside(),1,i,dim),dim);
		  Dune::FieldVector< ctype, dim > point = is -> geometry().center();
		

		   if(is->boundary()) 
		   
		{
	     // Dune::FieldVector< ctype, dim > point = is -> geometry().center();
	        for (int i=0; i < ref.size(is->indexInInside(),1,dim); i++)
                {
		  // and replace the associated line of A and b with a trivial one
		  int indexi = set.subIndex(*it,ref.subEntity(is->indexInInside(),1,i,dim),dim);
		  Dune::FieldVector< ctype, dim > point = is -> geometry().center();
                  rr_check=sqrt(pow(point[0],2)+pow(point[1],2));
			//Dune::FieldVector< ctype, dim > point1 = is -> geometry().corner();
                  //rr_check1=sqrt(pow(point1[0],2)+pow(point1[1],2));
		//  if((rr_check>(R2-1))||(rr_check<1))
 		if((rr_check>(R2-1)))
		 {     
			if (out==0)
			{
			
 			std::cout<<"Applying boundary conditions at radius   "<<rr_check<<std::endl; 
			out++;
			}
                     // C1[indexi][indexi] = 2*pi*R1;
		      //  C2[indexi][indexi] = 2*pi*R1*R1;
			

                       A[indexi] = 0.0; 
		      A[indexi][indexi] = 1.0;
		      B[indexi] = 0.0;

		 }
/*
out1=0;
if((rr_check<4))
		 {     
			if (out1==0)
			{
			
 			cout<<"Applying boundary conditions at radius   "<<rr_check<<endl; 
			out1++;
			}
                     
			

                       A[indexi] = 0.0; 
		      A[indexi][indexi] = 1.0;
		      B[indexi] = 0.0;

		 }

*/


		}
	      
               }



	      
            } 
        }
    }
}





  for(int i=0;i<N;i++)
	  {
	    if(A[i][i]==1)
	      {	   
		//	cout<<" i is "<<i<<endl;
		//	while(1){}
		for(int j=0;j<N;j++)
		  {
		     if(i!=j)
		       {
			if(A.exists(j,i))
			  A[j][i]=0;
			
			if(B.exists(j,i))
			  B[j][i]=0;
			 }
		  }
	      }
	  }


	double A_norm,B_norm;
	A_norm=(A.frobenius_norm2());
	  B_norm=(B.frobenius_norm2());
// 	   A/=(A_norm);
// 	   B/=(B_norm);
//	  cout<<" A and B norm "<<A_norm<<"   "<<B_norm<<endl;
// 	  A/=sqrt(A_norm);
// 	  B/=sqrt(B_norm);
	 
	
	Dune::writeMatrixToMatlab(A,"MatrixA");
	Dune::writeMatrixToMatlab(B,"MatrixB");
	//Dune::writeMatrixToMatlab(Cb,"MatrixCb");
        //Dune::writeMatrixToMatlab(C1,"MatrixC1");
	//Dune::writeMatrixToMatlab(C2,"MatrixC2");


	for (LeafIterator it = gv.template begin<0>(); it != itend; ++it)
	  {
	    Dune::GeometryType gt = it->type();
	    const Dune::template ReferenceElements<ctype,dim> &ref = Dune::ReferenceElements<ctype,dim>::general(gt);
	    int vertexsize = ref.size(dim);
	    for (int i = 0; i < vertexsize; i++)
	      {
		for (int j = 0; j < vertexsize; j++)
		  {
		    int indexi = set.subIndex(*it,i,dim); 
		    int indexj = set.subIndex(*it,j,dim);

 
		    MatSetValue(A_p,indexi,indexj,double(A[indexi][indexj]),INSERT_VALUES);
		    MatSetValue(B_p,indexi,indexj,double(B[indexi][indexj]),INSERT_VALUES);
		  }
	      }   
	  }


    MatAssemblyBegin(A_p,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A_p,MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(B_p,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(B_p,MAT_FINAL_ASSEMBLY);	
}


template<class GV, class Vector>
void SchrodingerFEM<GV,Vector>::solve(double m, double g)
{
 
  VertexMap vertexmap(gv);
  static const int N = gv.size(dim);
  
  // PetscViewer    viewer;
  PetscViewer view_out;
  PetscScalar    *avec;
 
  double k=8.61735E-5;
  double T=300;
  double Hr=27.21; 
  double beta=Hr/(k*T);
  double pai=3.142;
  double fermi;

        std::ofstream wf_plot;
	std::ofstream target;


        EPSCreate(MPI_COMM_WORLD,&eps);

        ierr = MatGetVecs(A_p,PETSC_NULL,&xr);
        ierr = MatGetVecs(A_p,PETSC_NULL,&xi);
	// PetscTruth Tr=PETSC_TRUE;
// 	EPSIsHermitian(eps,&Tr );
	//	EPSCreate(PETSC_COMM_WORLD,&eps); 

	EPSSetOperators(eps,A_p,B_p);
	EPSSetFromOptions(eps);
	EPSSetDimensions(eps, 1, PETSC_DECIDE, PETSC_DECIDE);	
	EPSSolve(eps);
	EPSGetIterationNumber(eps, &its);
 	//PetscPrintf(MPI_COMM_WORLD," Number of iterations of the method: %d\n",its);
  	//EPSGetOperationCounters(eps,PETSC_NULL,PETSC_NULL,&lits);
  	//PetscPrintf(MPI_COMM_WORLD," Number of linear iterations of the method: %d\n",lits);
  	EPSGetType(eps,&type);
  	//PetscPrintf(MPI_COMM_WORLD," Solution method: %s\n\n",type);
  	EPSGetDimensions(eps,&nev,PETSC_NULL,PETSC_NULL);
  	//PetscPrintf(MPI_COMM_WORLD," Number of requested eigenvalues: %d\n",nev);
  	EPSGetTolerances(eps,&tol,&maxit);
  	//PetscPrintf(MPI_COMM_WORLD," Stopping condition: tol=%.4g, maxit=%d\n",tol,maxit);
	EPSGetConverged(eps,&nconv);
 	//PetscPrintf(PETSC_COMM_WORLD," Number of converged approximate eigenpairs: %d\n\n",nconv);
	//PetscPrintf(PETSC_COMM_WORLD,
	//"           k             ||Ax-kBx||/||kx||\n"
	//"  --------------------- ------------------\n" );

double eigen[nconv];
	//Dune::VTKWriter<GV> vtkwriter_wf(gv,Dune::VTKOptions::conforming);
	Vector wf[nconv];
	for(int i=0; i<nconv; i++ )
	  {
	    ierr= EPSGetEigenpair(eps,i,&kr,&ki,xr,xi);
	    EPSComputeRelativeError(eps,i,&error);
            #if defined(PETSC_USE_COMPLEX)
	    re = PetscRealPart(kr);
	    im = PetscImaginaryPart(kr);
            #else
	    re = kr;
	    im = ki;
            #endif
	    //    while(1){}
	    eigen[i]=re;
	    
 	    //if( im != 0.0 ) 
 	      {
 		ierr = PetscPrintf(PETSC_COMM_WORLD," % 6f %+12g i",re*27.2,im*27.2);
 	      } 
 	    //else 
 	      {
 		//PetscPrintf(PETSC_COMM_WORLD,"    Real Part of Eigen   % 6f  eV    ",re*27.21);   // output in eV
 	      }
 	    PetscPrintf(PETSC_COMM_WORLD," % 12g   %d\n",error,i);
	    
	    
	    
	    
	    VecGetArray(xr,&avec);
	    for(int r=0; r<N; r++)
	      {
		wavevector[r]= PetscRealPart(*avec);//*(*avec);
		avec++;     
	      }
	    
	    
	    VecRestoreArray(xr,&avec);

 //wavefunction(gv,wavevector)

	    /*
	 //   wavefunction_bc<GV,Vector>;
	    wf[i]=wavevector;
	    if(i<40)
	    //  wavefunction_bc.wf_VTKwrite(vtkwriter, wf[i], i);
// 	    wavefunction.write_wf_gnuplot("wf_plot.dat");
// 	    while(1){}
	    if(wavefunction_bc.is_y_fundamental())
	      {
		//cout<<"this is a true mode"<<endl;
		wavefunction_bc.only_x_dependant();
		wavefunction_bc.pnormalize();
		wavefunction_bc.write_wf_gnuplot("wf_plot.dat");
		
     
	    
	    //    }
		double Ey=pai*pai/(2*m*50*50);
		if(carrier_type==0)
		  fermi=(g*m/(2*pai*beta))*log(1+exp((Ef-(re-Ey))*beta));
		if(carrier_type==1)
		  {
		    if(m==0.29)
		      {
			std::cout<<"Split off hole"<<std::endl;
			re=re+0.044/Hr;
		      }
		    fermi=(g*m/(2*pai*beta))*log(1+exp((Ef-(re-Ey))*beta));
		  }
			 //fermi=(g*m/beta)*(2*pai)*log(1+exp((Ef-(re-Ey))*beta));
		// 	     if(fabs(re)<0.6/27.211)
		// 	       {
// 		fermi=g/(1+(exp((re-Ef)*beta)));
		for(int p=0;p<N;p++)
		  {
		    charge_tmp[p]+=fermi*wavevector[p]*wavevector[p]*50;
		    //cout<<charge_tmp[p]<<endl;
		    // 		   }
		  } 
	      }
	  }
//	vtkwriter_wf.write("wavefunction", Dune::VTK::binaryappended);
	//	std::cout << "visualizing... wavefunction" << "\n";
//	 Dune::VTKWriter<GridType::LeafGridView> vtkwriter_wf(grid->leafGridView());
	//vtkwriter_wf.write("wavefunction", Dune::VTK::appendedraw);
	
	

  
  
  
  
  
target.open("target.dat");
        for (int kk=0;kk<nconv;kk++)
            target<<eigen[kk]<< std::endl; // write in hartree
	target.close();


*/

}
}
	

template<class GV, class Vector>
void SchrodingerFEM<GV,Vector>::apply()
{
//   determineAdjacencyPattern();

  double me[2]={0.916,0.19};
  double me_d[2]={0.19,0.417};
  double mh[3]={0.291,0.2,0.29};
  double mh_d[3]={0.645,0.251,0.29};
  double ge[2]={2,4};
  double gh[3]={1,1,1};
  charge_tmp=0.0;
  std::cout<<"carrier type "<<carrier_type<<"Ef = "<<Ef<<std::endl;
 // if(carrier_type==0)
    {

      //     cout<<"here"<<endl;
      for(int i=0;i<1;i++)
	{
    
	  assemble(me[i]);
	  std::cout<<"Assembled   "<<std::endl; 
	  apply_bc();
	//  solve(me_d[i],ge[i]);

// 	  assemble(1);
// 	  apply_bc();
std::cout<<"Applying boundary conditions    "<<std::endl; 
 	  solve(1,1);
std::cout<<"solving "<<std::endl; 
	}
    }

  if(carrier_type==1)
    {
      for(int i=0;i<1;i++)
	{
	  assemble(mh[i]);
	  apply_bc();
	  std::cout<<"here is"<<i; 
	  //solve(mh_d[i],gh[i]);
	}
    }

  charge=charge_tmp;

	//	cout<<gv.size(dim)<<endl;

       
}
   
   
   //....changes 1end
   
   #endif	// SCHRODINGERFEM_HH
   
   



