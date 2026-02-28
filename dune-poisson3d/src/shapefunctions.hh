/**
 * P1 (linear) shape functions for simplices in 2D (triangles) and 3D (tetrahedra).
 *
 * For a simplex with (dim+1) vertices, the P1 basis functions are:
 *   phi_0 = 1 - xi_1 - xi_2 - ... - xi_dim
 *   phi_i = xi_i   for i = 1, ..., dim
 *
 * This file provides:
 *   - P1ShapeFunction<CT,RT,dim>      : a single basis function
 *   - P1ShapeFunctionSet<CT,RT,dim>   : the set of (dim+1) basis functions
 *
 * Compatible with DUNE's FieldVector for positions and gradients.
 */

#ifndef SHAPEFUNCTIONS_HH
#define SHAPEFUNCTIONS_HH

#include <dune/common/fvector.hh>

// ─── Single P1 shape function ───────────────────────────────────────────────
template<class CT, class RT, int dim>
class P1ShapeFunction
{
public:
  /**
   * @param index  Vertex index in the reference simplex (0 .. dim)
   */
  P1ShapeFunction() : index_(0) {}
  explicit P1ShapeFunction(int index) : index_(index) {}

  /** Evaluate the shape function at local coordinate xi */
  RT evaluateFunction(const Dune::FieldVector<CT, dim>& xi) const
  {
    if (index_ == 0) {
      RT val = 1.0;
      for (int d = 0; d < dim; ++d)
        val -= xi[d];
      return val;
    } else {
      return xi[index_ - 1];
    }
  }

  /** Evaluate the gradient (in reference coordinates) at local coordinate xi */
  Dune::FieldVector<RT, dim> evaluateGradient(const Dune::FieldVector<CT, dim>& /*xi*/) const
  {
    Dune::FieldVector<RT, dim> grad(0.0);
    if (index_ == 0) {
      for (int d = 0; d < dim; ++d)
        grad[d] = -1.0;
    } else {
      grad[index_ - 1] = 1.0;
    }
    return grad;
  }

private:
  int index_;
};


// ─── Set of all (dim+1) P1 shape functions on a simplex ─────────────────────
template<class CT, class RT, int dim>
class P1ShapeFunctionSet
{
public:
  static const int numBasis = dim + 1;

  static P1ShapeFunctionSet& instance()
  {
    static P1ShapeFunctionSet singleton;
    return singleton;
  }

  const P1ShapeFunction<CT, RT, dim>& operator[](int i) const
  {
    return basis_[i];
  }

  int size() const { return numBasis; }

private:
  P1ShapeFunctionSet()
  {
    for (int i = 0; i < numBasis; ++i)
      basis_[i] = P1ShapeFunction<CT, RT, dim>(i);
  }

  P1ShapeFunction<CT, RT, dim> basis_[numBasis];
};

#endif // SHAPEFUNCTIONS_HH
