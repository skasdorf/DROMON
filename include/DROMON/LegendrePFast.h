//
// Created by Jake J. Harmon (jake.harmon@ieee.org) on 5/1/22.
//

#ifndef DROMON_LEGENDREPFAST_H
#define DROMON_LEGENDREPFAST_H
#include "config.h"

DROMON_NAMESPACE_OPEN

namespace utility
{
// Legendre polynomials computed in their Horner form
template <class Real>
Real legendre_p_0(const Real&)
{
  return Real(1);
}
template <class Real>
Real legendre_p_1(const Real& x)
{
  return x;
}
template <class Real>
Real legendre_p_2(const Real& x)
{
  return Real(1.5)*x*x - Real(0.5);
}
template <class Real>
Real legendre_p_3(const Real& x)
{
  return x*(Real(2.5)*x*x - Real(1.5));
}
template <class Real>
Real legendre_p_4(const Real& x)
{
  return x*x*(Real(4.375)*x*x - Real(3.75)) + Real(0.375);
}
template <class Real>
Real legendre_p_5(const Real& x)
{
  return x*(x*x*(Real(7.875)*x*x - Real(8.75)) + Real(1.875));
}

// Similarly, we now have functions for computing the Q functions from the MaxOrthoBasis

template <class Real>
Real legendre_q_0(const Real& x)
{
  return (Real)(1) - x;
}
template <class Real>
Real legendre_q_1(const Real& x)
{
  return (Real)(1) + x;
}
template <class Real>
Real legendre_q_2(const Real& x)
{
  return Real(1.5)*(x*x - Real(1.0));
}
template <class Real>
Real legendre_q_3(const Real& x)
{
  return Real(2.5)*x*(x*x - Real(1.0));
}
template <class Real>
Real legendre_q_4(const Real& x)
{
  return x*x*(Real(4.375)*x*x - Real(5.0)) + Real(0.625);
}


template <class Real>
class LegendrePContainer
{
public:
  using _legendre_p_function = Real(*) (const Real&);
  Real operator()(const unsigned int& order, const Real& x) const
  {
    #ifdef DEBUG
    assert(order < n_orders && "Order exceeds the currently provided number of functions in this container!");
    #endif
    return legendre_polys[order](x);
  }
private:
  _legendre_p_function legendre_polys [6]
  = {legendre_p_0<Real>, legendre_p_1<Real>,legendre_p_2<Real>,legendre_p_3<Real>,legendre_p_4<Real>, legendre_p_5<Real>};
  const unsigned int n_orders = 6;
};


}
DROMON_NAMESPACE_CLOSE
#endif // DROMON_LEGENDREPFAST_H
