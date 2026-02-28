#ifndef WAVEFUNCTION_HH
#define WAVEFUNCTION_HH

#include <dune/grid/common/mcmgmapper.hh>
#include <dune/geometry/referenceelements.hh>



template<class GV, class Vector>
class WaveFunction
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
  

private: 
 const GV& gv;
 Vector& wavevector;
 const LeafIndexSet& set;
  double max; 
  Dune::FieldVector<ctype,dim> max_coords;
  void find_maximum();
 
public:
  WaveFunction(const GV& gv_, Vector& wavevector_);
  void pnormalize();
  bool is_y_fundamental();
  void only_x_dependant();
  void write_wf_gnuplot(const char* filename);
  void wf_VTKwrite(Dune::VTKWriter<GV>& vtkwriter_wf, Vector& wf, int i);

  
};



template<class GV, class Vector>
WaveFunction<GV,Vector>::WaveFunction(const GV& gv_, Vector& wavevector_):gv(gv_),wavevector(wavevector_),set(gv.indexSet())
{
  max=0;
  find_maximum();

}


template<class GV, class Vector>
void WaveFunction<GV,Vector>::find_maximum()
  {
    VertexMap vertexmap(gv, Dune::mcmgVertexLayout());
    VertexIterator vtend=gv.template end<dim>();
    for(VertexIterator vt=gv.template begin<dim>();vt!=vtend;++vt)
      {
	if(fabs(wavevector[vertexmap.index(*vt)])>max)
	  {
	    max=fabs(wavevector[vertexmap.index(*vt)]);
	    max_coords=vt->geometry().corner(0);
	    //	    max_y=vt->geometry().corner(0)[1];
	    // cout<<"wavefunction maximum "<<fabs(wavevector[vertexmap.index(*vt)])<<" @  "<<max_coords<<std::endl;
	    
	  }
      }
  }


template<class GV, class Vector>
bool WaveFunction<GV,Vector>::is_y_fundamental()
{

/*  This function determines wether wavefunction has fundamental mode in y-direction
*/
  bool true_mode=1;
  int sign_tmp=0;
  int sign=0;
  VertexMap vertexmap(gv, Dune::mcmgVertexLayout());
  VertexIterator vtend=gv.template end<dim>();
  for(VertexIterator vt=gv.template begin<dim>();vt!=vtend;++vt)
    {
      if(fabs(vt->geometry().corner(0)[0]-max_coords[0])<1e-8)
	{
	  if((fabs(vt->geometry().corner(0)[1])!=0)&&(fabs(vt->geometry().corner(0)[1])!=100))
	    {
	      if(wavevector[vertexmap.index(*vt)]<0)
		sign=-1;
	      if(wavevector[vertexmap.index(*vt)]>0)
		sign=1;
	      sign_tmp=sign;
	      break;
	    }
	}
      
    }

  for(VertexIterator vt=gv.template begin<dim>();vt!=vtend;++vt)
    {
      if(fabs(vt->geometry().corner(0)[0]-max_coords[0])<1e-8)
	{
	  if((fabs(vt->geometry().corner(0)[1])!=0)&&(fabs(vt->geometry().corner(0)[1])!=100))
	    {
	      if(wavevector[vertexmap.index(*vt)]<0)
		sign=-1;
	      if(wavevector[vertexmap.index(*vt)]>0)
		sign=1;
	      if(sign_tmp!=sign)
		true_mode=0;
	    }
	}
      
    }

  return true_mode;

//   if(true_mode=1)
//     return true;
//   if(true_mode=0)
//     return false;

  
}

template<class GV, class Vector>
void WaveFunction<GV,Vector>::only_x_dependant()
{
  VertexMap vertexmap(gv, Dune::mcmgVertexLayout());
  VertexIterator vtend=gv.template end<dim>();
  for(VertexIterator vt=gv.template begin<dim>();vt!=vtend;++vt)
    {
      for(VertexIterator vjt=gv.template begin<dim>();vjt!=vtend;++vjt)
	{
	  if(fabs(vjt->geometry().corner(0)[1]-max_coords[1])<1e-8)
	    {
	      if(fabs(vt->geometry().corner(0)[0]-vjt->geometry().corner(0)[0])<1e-8)
		{
		  wavevector[vertexmap.index(*vt)]=wavevector[vertexmap.index(*vjt)];
		}
	    }
	}
    }
  
}

template<class GV, class Vector>
void WaveFunction<GV,Vector>::pnormalize()
{

 double wf_avg=0.0;
 double norm=0.0;
 const LeafIterator itend = gv.template end<0>();
 for (LeafIterator it = gv.template begin<0>(); it != itend; ++it) /*@\label{fem:loop1}@*/
   {
     // determine geometry type of the current element and get the matching reference element
     Dune::GeometryType gt = it->type();
        const auto& ref = Dune::ReferenceElements<ctype,dim>::general(gt);
        int vertexsize = ref.size(dim);
   
        // get a quadrature rule of order one for the given geometry type
        const Dune::QuadratureRule<ctype,dim>& rule = Dune::QuadratureRules<ctype,dim>::rule(gt,1);
        for (typename Dune::QuadratureRule<ctype,dim>::const_iterator r = rule.begin(); r != rule.end() ; ++r)
        {
	  // get the weight at the current quadrature point
	  ctype weight = r->weight();
	  // compute Jacobian determinant for the transformation formula
	  ctype detjac = it->geometry().integrationElement(r->position());

 	  wf_avg=0.0;
	  for(int i=0;i<vertexsize;i++)
	    wf_avg+=wavevector[set.subIndex((*it),i,dim)]*wavevector[set.subIndex((*it),i,dim)]/vertexsize;
	  
	  norm+=wf_avg*weight*detjac;

	  
	}
	
   }
 norm=sqrt(norm);
 wavevector/=norm;
}


template<class GV, class Vector>
void WaveFunction<GV,Vector>::write_wf_gnuplot(const char* filename)
{
  VertexMap vertexmap(gv, Dune::mcmgVertexLayout());
  VertexIterator vtend=gv.template end<dim>();
  std::ofstream wf_plot;
  wf_plot.open(filename);
  // VertexIterator vtend=gv. template end<dim>();
  for(VertexIterator vt=gv.template begin<dim>();vt!=vtend;++vt)
    {
      wf_plot<<vt->geometry().corner(0)[0]<<"    "<<vt->geometry().corner(0)[1]<<"    "<<0.0<<"    "<<wavevector[vertexmap.index(*vt)]<<std::endl;
    }
  wf_plot.close();
}


template<class GV, class Vector>
void WaveFunction<GV,Vector>::wf_VTKwrite(Dune::VTKWriter<GV>& vtkwriter_wf, Vector& wf, int i)
{
				 std::string wf_str("wavefunction");
				 std::stringstream ss;
				 ss<<i;
				 wf_str+=ss.str();
				 vtkwriter_wf.addVertexData(wf, wf_str);
}

#endif	// WAVEFUNCTION_HH
