//
// Created by Jake J. Harmon (jake.harmon@ieee.org) on 4/30/22.
//

#ifndef DROMON_DIRECTFN_ET_BOUNDS_FUNCTIONS_H
#define DROMON_DIRECTFN_ET_BOUNDS_FUNCTIONS_H
#include "config.h"
#include "Point.h"
DROMON_NAMESPACE_OPEN
// Bounds functions for theta
template <class Real>
void ET_theta_bounds_0_or_1(Real &theta_lower, Real &theta_upper) {
  theta_lower = Real(0.);
  theta_upper = Real(constants<Real>::PI) / Real(4.);
}
template <class Real>
void ET_theta_bounds_2_or_3(Real &theta_lower, Real &theta_upper) {
  theta_lower = Real(constants<Real>::PI) / Real(4.);
  theta_upper = Real(constants<Real>::PI) / Real(2.);
}
template <class Real>
void ET_theta_bounds_4_or_5(Real &theta_lower, Real &theta_upper) {
  theta_lower = Real(constants<Real>::PI) / Real(2.);
  theta_upper = Real(3.) * constants<Real>::PI / Real(4.);
}
template <class Real>
void ET_theta_bounds_6_or_7(Real &theta_lower, Real &theta_upper) {
  theta_lower = Real(3.) * constants<Real>::PI / Real(4.);
  theta_upper = constants<Real>::PI;
}

// Bounds functions for psi
template <class Real>
void ET_psi_bounds_0(const Real &atan_cos_theta, const Real &, Real &psi_lower,
                     Real &psi_upper) {
  psi_lower = Real(0.);
  psi_upper = atan_cos_theta;
}
template <class Real>
void ET_psi_bounds_1(const Real &atan_cos_theta, const Real &, Real &psi_lower,
                     Real &psi_upper) {
  psi_lower = atan_cos_theta;
  psi_upper = constants<Real>::PI / Real(2.);
}
template <class Real>
void ET_psi_bounds_2(const Real &, const Real &atan_sin_theta, Real &psi_lower,
                     Real &psi_upper) {
  psi_lower = Real(0.);
  psi_upper = atan_sin_theta;
}
template <class Real>
void ET_psi_bounds_3(const Real &, const Real &atan_sin_theta, Real &psi_lower,
                     Real &psi_upper) {
  psi_lower = atan_sin_theta;
  psi_upper = constants<Real>::PI / Real(2.);
}
template <class Real>
void ET_psi_bounds_4(const Real &, const Real &atan_sin_theta, Real &psi_lower,
                     Real &psi_upper) {
  psi_lower = Real(0.);
  psi_upper = atan_sin_theta;
}
template <class Real>
void ET_psi_bounds_5(const Real &, const Real &atan_sin_theta, Real &psi_lower,
                     Real &psi_upper) {
  psi_lower = atan_sin_theta;
  psi_upper = constants<Real>::PI / Real(2.);
}
template <class Real>
void ET_psi_bounds_6(const Real &atan_cos_theta, const Real &, Real &psi_lower,
                     Real &psi_upper) {
  psi_lower = Real(0.);
  psi_upper = -atan_cos_theta;
}
template <class Real>
void ET_psi_bounds_7(const Real &atan_cos_theta, const Real &, Real &psi_lower,
                     Real &psi_upper) {
  psi_lower = -atan_cos_theta;
  psi_upper = constants<Real>::PI / Real(2.);
}

// For the u bounds
template <class Real>
void ET_u_bounds_0_or_6(const Real &, const Real &, const Real&,
                        Real &u_lower, Real &u_upper) {
  u_lower = Real(-1.);
  u_upper = Real(1.);
}
template <class Real>
void ET_u_bounds_1_lower(const Real &cos_theta, const Real &,
                         const Real &tan_psi, Real &u_lower, Real &u_upper) {
  u_lower = Real(-1.);
  u_upper = Real(2.) * cos_theta / tan_psi - Real(1.);
}
template <class Real>
void ET_u_bounds_2_lower(const Real &, const Real &tan_theta,
                         const Real &, Real &u_lower, Real &u_upper) {
  u_lower = Real(-1.);
  u_upper = Real(2.0) / tan_theta - Real(1.0);
}
template <class Real>
void ET_u_bounds_3_lower(const Real &cos_theta, const Real &,
                         const Real &tan_psi, Real &u_lower, Real &u_upper) {
  u_lower = Real(-1.);
  u_upper = Real(2.) * cos_theta / tan_psi - Real(1.); // u_1_theta
}
template <class Real>
void ET_u_bounds_4_lower(const Real &, const Real &tan_theta,
                         const Real &, Real &u_lower, Real &u_upper) {
  u_lower = Real(-1.);
  u_upper = Real(2.) / tan_theta + Real(1.); // u_2_theta;
}
template <class Real>
void ET_u_bounds_5_or_7_lower(const Real &cos_theta, const Real &,
                              const Real &tan_psi, Real &u_lower,
                              Real &u_upper) {
  u_lower = Real(-1.);
  u_upper = Real(2.) * cos_theta / tan_psi + Real(1.); // u_2_psi;
}

template <class Real>
void ET_u_bounds_1_upper(const Real &cos_theta, const Real &,
                         const Real &tan_psi, Real &u_lower, Real &u_upper) {
  u_lower = Real(2.) * cos_theta / tan_psi - Real(1.);
  u_upper = Real(1.);
}
template <class Real>
void ET_u_bounds_2_upper(const Real &, const Real &tan_theta,
                         const Real &, Real &u_lower, Real &u_upper) {
  u_lower = Real(2.) / tan_theta - Real(1.); // u_1_theta;
  u_upper = Real(1.);
}
template <class Real>
void ET_u_bounds_3_upper(const Real &cos_theta, const Real &,
                         const Real &tan_psi, Real &u_lower, Real &u_upper) {
  u_lower = Real(2.) * cos_theta / tan_psi - Real(1.); // u_1_psi;
  u_upper = Real(1.);
}
template <class Real>
void ET_u_bounds_4_upper(const Real &, const Real &tan_theta,
                         const Real &, Real &u_lower, Real &u_upper) {
  u_lower = Real(2.) / tan_theta + Real(1.);
  u_upper = Real(1.);
}
template <class Real>
void ET_u_bounds_5_or_7_upper(const Real &cos_theta, const Real &,
                              const Real &tan_psi, Real &u_lower,
                              Real &u_upper) {
  u_lower = Real(2.) * cos_theta / tan_psi + Real(1.);
  u_upper = Real(1.);
}

// For the lambda bounds
template <class Real>
void ET_lambda_bounds_0(const Real &cos_theta, const Real &sin_theta,
                                 const Real &cos_psi, const Real &sin_psi,
                                 const Real &u, Real &lambda_lower,
                                 Real &lambda_upper) {
  lambda_lower = Real(0.);
  lambda_upper = (u + Real(1.)) / (cos_theta * cos_psi);
}

template <class Real>
void ET_lambda_bounds_6(const Real &cos_theta, const Real &sin_theta,
                        const Real &cos_psi, const Real &sin_psi,
                        const Real &u, Real &lambda_lower,
                        Real &lambda_upper) {
  lambda_lower = Real(0.);
  lambda_upper =  (u - Real(1.)) / (cos_theta * cos_psi);
}

template <class Real>
void ET_lambda_bounds_1_2_or_3_region_0(const Real &cos_theta, const Real &sin_theta,
                        const Real &cos_psi, const Real &sin_psi,
                        const Real &u, Real &lambda_lower,
                        Real &lambda_upper) {
  lambda_lower = Real(0.);
  lambda_upper =  (u + Real(1.0)) / (cos_theta * cos_psi);
}

template <class Real>
void ET_lambda_bounds_4_region_0(const Real &cos_theta, const Real &sin_theta,
                                        const Real &cos_psi, const Real &sin_psi,
                                        const Real &u, Real &lambda_lower,
                                        Real &lambda_upper) {
  lambda_lower = Real(0.);
  lambda_upper =  Real(2.) / (sin_theta * cos_psi);
}

template <class Real>
void ET_lambda_bounds_5_region_0(const Real &cos_theta, const Real &sin_theta,
                                 const Real &cos_psi, const Real &sin_psi,
                                 const Real &u, Real &lambda_lower,
                                 Real &lambda_upper) {
  lambda_lower = Real(0.);
  lambda_upper =  Real(2.) / sin_psi;
}

template <class Real>
void ET_lambda_bounds_7_region_0(const Real &cos_theta, const Real &sin_theta,
                                 const Real &cos_psi, const Real &sin_psi,
                                 const Real &u, Real &lambda_lower,
                                 Real &lambda_upper) {
  lambda_lower = Real(0.);
  lambda_upper =  (u - Real(1.)) / (cos_theta * cos_psi);
}

template <class Real>
void ET_lambda_bounds_1_3_or_7_region_1(const Real &cos_theta, const Real &sin_theta,
                                        const Real &cos_psi, const Real &sin_psi,
                                        const Real &u, Real &lambda_lower,
                                        Real &lambda_upper) {
  lambda_lower = Real(0.);
  lambda_upper =  Real(2.) / sin_psi;
}

template <class Real>
void ET_lambda_bounds_2_region_1(const Real &cos_theta, const Real &sin_theta,
                                 const Real &cos_psi, const Real &sin_psi,
                                 const Real &u, Real &lambda_lower,
                                 Real &lambda_upper) {
  lambda_lower = Real(0.);
  lambda_upper =  Real(2.) / (sin_theta * cos_psi);
}

template <class Real>
void ET_lambda_bounds_4_or_5_region_1(const Real &cos_theta, const Real &sin_theta,
                                 const Real &cos_psi, const Real &sin_psi,
                                 const Real &u, Real &lambda_lower,
                                 Real &lambda_upper) {
  lambda_lower = Real(0.);
  lambda_upper =  (u - Real(1.)) / (cos_theta * cos_psi);
}

// For the u_and_v_coords mapping

// Since the algorithm expects a particular orientation of the two cells with
// adjacent edges, we must update these intermediate uv_test and uv_trial
// values to the correct ones

/* The expected orientation is

----3----
|       | uv_test
2 |       | 1
|___0___|
|   0   |
| 1   2 |
|       | uv_trial
----3----
 */
template <class Real>
void ET_map_u_and_v_coords_0(Point<2, Real>&)
{
  return;
}

template <class Real>
void ET_map_u_and_v_coords_1(Point<2, Real>& uv)
{
  Real temp = uv(0);
  uv[0] = -uv(1);
  uv[1] = temp;
}
template <class Real>
void ET_map_u_and_v_coords_2(Point<2, Real>& uv)
{
  Real temp = uv(0);
  uv[0] = uv(1);
  uv[1] = -temp;
}
template <class Real>
void ET_map_u_and_v_coords_3(Point<2, Real>& uv)
{
  uv[0] *= -1;
  uv[1] *= -1;
}
DROMON_NAMESPACE_CLOSE
#endif // DROMON_DIRECTFN_ET_BOUNDS_FUNCTIONS_H
