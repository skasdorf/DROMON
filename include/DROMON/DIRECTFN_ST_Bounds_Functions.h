//
// Created by Jake J. Harmon (jake.harmon@ieee.org) on 4/29/22.
//

#ifndef DROMON_DIRECTFN_BOUNDS_FUNCTIONS_H
#define DROMON_DIRECTFN_BOUNDS_FUNCTIONS_H
#include "config.h"
#include "Point.h"
#include <boost/math/special_functions/sin_pi.hpp>
#include <boost/math/special_functions/cos_pi.hpp>

#include <gsl/gsl_sf_trig.h>

DROMON_NAMESPACE_OPEN
template <class Real>
void ST_psi_bounds_0(Real& psi_lower, Real& psi_upper)
{
  psi_lower = Real(0);
  psi_upper = constants<Real>::PI/Real(4.);
}

template <class Real>
void ST_psi_bounds_1_or_2(Real& psi_lower, Real& psi_upper)
{
  psi_lower = constants<Real>::PI/Real(4.);
  psi_upper = constants<Real>::PI/Real(2.);
}

template <class Real>
void ST_psi_bounds_3_or_4(Real& psi_lower, Real& psi_upper)
{
  psi_lower = constants<Real>::PI/Real(2.);
  psi_upper = Real(3.) * constants<Real>::PI / Real(4.);
}

template <class Real>
void ST_psi_bounds_5(Real& psi_lower, Real& psi_upper)
{
  psi_lower = Real(3.) * constants<Real>::PI / Real(4.);
  psi_upper = constants<Real>::PI;
}

template <class Real>
void ST_u_bounds_0_or_5(const Real&, Real& u_lower, Real& u_upper)
{
  u_lower = Real(-1.);
  u_upper = Real(1.);
}
template <class Real>
void ST_u_bounds_1(const Real& tan_psi_less_pi_half, Real& u_lower, Real& u_upper)
{
  u_lower = Real(2.) * tan_psi_less_pi_half + Real(1);
  u_upper = Real(1.);
}
template <class Real>
void ST_u_bounds_2(const Real& tan_psi_less_pi_half, Real& u_lower, Real& u_upper)
{
  u_lower = Real(-1.);
  u_upper = Real(2.) * tan_psi_less_pi_half + Real(1);
}
template <class Real>
void ST_u_bounds_3(const Real& tan_psi_less_pi_half, Real& u_lower, Real& u_upper)
{
  u_lower = Real(2.) * tan_psi_less_pi_half - Real(1);
  u_upper = Real(1.);
}
template <class Real>
void ST_u_bounds_4(const Real& tan_psi_less_pi_half, Real& u_lower, Real& u_upper)
{
  u_lower = Real(-1.);
  u_upper = Real(2.) * tan_psi_less_pi_half - Real(1);
}

template <class Real>
void ST_lambda_bounds_0_or_1(const Real& cos_psi, const Real&, const Real& u, Real& lambda_lower, Real& lambda_upper, const Real& psi)
{
  lambda_lower = Real(0.0);
  lambda_upper = (Real(1.) - u) / cos_psi;
}

template <class Real>
void ST_lambda_bounds_2_or_3(const Real&, const Real& sin_psi, const Real& u, Real& lambda_lower, Real& lambda_upper, const Real& psi)
{
  lambda_lower = Real(0.0);
  lambda_upper = (Real(2.)) / sin_psi;
}

template <class Real>
void ST_lambda_bounds_4_or_5(const Real& cos_psi, const Real&, const Real& u, Real& lambda_lower, Real& lambda_upper, const Real& psi)
{
//  const auto test = -1.0/cos_psi; // Second most accurate
//  const auto test2 = 1.0/cos(constants<Real>::PI - psi);
//  const auto test3 = 1.0/boost::math::cos_pi((1.0 - psi/constants<Real>::PI));
//  const auto test4 = -1.0/gsl_sf_cos(psi); // Most accurate
//  const auto test5 = 1.0/gsl_sf_cos(constants<Real>::PI - psi);

  lambda_lower = Real(0.0);
  //lambda_upper = -(u + Real(1.0)) / cos_psi;
  lambda_upper = (u+Real(1.0))/cos(constants<Real>::PI-psi);
}

template <class Real>
void ST_rho_bounds(const Real& lambda, Real& rho_lower, Real& rho_upper)
{
  rho_lower = Real(0.0);
  rho_upper = lambda;
}

template <class Real>
Point<2, Real> ST_map_u_and_v_coords_0(const Point<2, Real>& uv)
{
  return uv;
}

template <class Real>
Point<2, Real> ST_map_u_and_v_coords_1(const Point<2, Real>& uv)
{
  Point<2, Real> out;
  out[0] = Real(-1.) * uv(1);
  out[1] = uv(0);
  return out;
}

template <class Real>
Point<2, Real> ST_map_u_and_v_coords_2(const Point<2, Real>& uv)
{
  Point<2, Real> out;
  out[0] = Real(-1.) * uv(0);
  out[1] = Real(-1.) * uv(1);
  return out;
}

template <class Real>
Point<2, Real> ST_map_u_and_v_coords_3(const Point<2, Real>& uv)
{
  Point<2, Real> out;
  out[0] = uv(1);
  out[1] = Real(-1.) * uv(0);
  return out;
}
DROMON_NAMESPACE_CLOSE
#endif // DROMON_DIRECTFN_BOUNDS_FUNCTIONS_H
