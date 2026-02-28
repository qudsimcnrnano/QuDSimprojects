
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
