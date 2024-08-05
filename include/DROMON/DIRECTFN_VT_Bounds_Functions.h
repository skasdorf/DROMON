//
// Created by Jake J. Harmon (jake.harmon@ieee.org) on 5/1/22.
//

#ifndef DROMON_DIRECTFN_VT_BOUNDS_FUNCTIONS_H
#define DROMON_DIRECTFN_VT_BOUNDS_FUNCTIONS_H
#include "config.h"
#include "Point.h"

DROMON_NAMESPACE_OPEN
template <class Real>
void VT_psi_bounds_I(const Real& atan_Lq_over_lp, Real& psi_lower, Real& psi_upper)
{
  psi_lower = Real(0.0);
  psi_upper = atan_Lq_over_lp;
}

template <class Real>
void VT_psi_bounds_II(const Real& atan_Lq_over_lp, Real& psi_lower, Real& psi_upper)
{
  psi_lower = atan_Lq_over_lp;
  psi_upper = constants<Real>::PI / Real(2.0);
}

template <class Real>
void VT_Ltest_0_or_1(const Real& cos_theta_test, const Real&, Real& lp)
{
  lp = Real(2.0)/cos_theta_test;
}
template <class Real>
void VT_Ltest_2_or_3(const Real&, const Real& sin_theta_test, Real& lp)
{
  lp = Real(2.0)/sin_theta_test;
}

template <class Real>
void VT_Ltrial_0_or_2(const Real& cos_theta_trial, const Real&, Real& lp)
{
  lp = Real(2.0)/cos_theta_trial;
}
template <class Real>
void VT_Ltrial_1_or_3(const Real&, const Real& sin_theta_trial, Real& lp)
{
  lp = Real(2.0)/sin_theta_trial;
}

template <class Real>
void VT_theta_first(Real& theta_lower, Real&  theta_upper)
{
  theta_lower = Real(0.0);
  theta_upper = constants<Real>::PI / Real(4.0);
}
template <class Real>
void VT_theta_second(Real& theta_lower, Real&  theta_upper)
{
  theta_lower = constants<Real>::PI / Real(4.0);
  theta_upper = constants<Real>::PI / Real(2.0);
}

// For mapping the u and v coordinates based on the orientation of the adjacency
template <class Real>
void VT_map_u_and_v_coords_0(Point<2, Real>&)
{
  return;
}
template <class Real>
void VT_map_u_and_v_coords_1(Point<2, Real>& uv)
{
  const Real temp = uv(0);
  uv[0] = -uv(1);
  uv[1] = temp;
}
template <class Real>
void VT_map_u_and_v_coords_2(Point<2, Real>& uv)
{
  const Real temp = uv(0);
  uv[0] = uv(1);
  uv[1] = -temp;
}
template <class Real>
void VT_map_u_and_v_coords_3(Point<2, Real>& uv)
{
  uv[0] *= Real(-1.0);
  uv[1] *= Real(-1.0);
}

DROMON_NAMESPACE_CLOSE
#endif // DROMON_DIRECTFN_VT_BOUNDS_FUNCTIONS_H
