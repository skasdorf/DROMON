//
// Created by Jake J. Harmon (jake.harmon@ieee.org) on 4/30/22.
// Edited by Christopher A. Erickson (Christopher.Erickson@ieee.org) on 7/2/25
//
/**
 * @file
 * @defgroup ETBounds
 * @defgroup ETMapping
 * @brief Defines helper functions for computing angular and parametric bounds
 *        in element integration regions for direct function evaluations.
 *
 * These functions partition theta, psi, u, and lambda integration ranges into subregions
 * 
 * 
 * Also includes Functions for reorienting u and v parametric coordinates.
 *
 * These transformations adjust local element coordinates to enforce consistent orientation between adjacent cells.
 */
#ifndef DROMON_DIRECTFN_ET_BOUNDS_FUNCTIONS_H
#define DROMON_DIRECTFN_ET_BOUNDS_FUNCTIONS_H
#include "config.h"
#include "Point.h"
DROMON_NAMESPACE_OPEN
// Bounds functions for theta
/**
 * @brief Computes the lower and upper bounds for theta in region 0 or 1.
 * @ingroup ETBounds
 * @tparam Real Floating point type.
 * @param[out] theta_lower Lower bound of theta.
 * @param[out] theta_upper Upper bound of theta.
 */
template <class Real>
void ET_theta_bounds_0_or_1(Real &theta_lower, Real &theta_upper) {
  theta_lower = Real(0.);
  theta_upper = Real(constants<Real>::PI) / Real(4.);
}
/**
 * @brief Computes the lower and upper bounds for theta in region 2 or 3.
 * @ingroup ETBounds
 * @tparam Real Floating point type.
 * @param[out] theta_lower Lower bound of theta.
 * @param[out] theta_upper Upper bound of theta.
 */
template <class Real>
void ET_theta_bounds_2_or_3(Real &theta_lower, Real &theta_upper) {
  theta_lower = Real(constants<Real>::PI) / Real(4.);
  theta_upper = Real(constants<Real>::PI) / Real(2.);
}
/**
 * @brief Computes the lower and upper bounds for theta in region 4 or 5.
 * @ingroup ETBounds
 * @tparam Real Floating point type.
 * @param[out] theta_lower Lower bound of theta.
 * @param[out] theta_upper Upper bound of theta.
 */
template <class Real>
void ET_theta_bounds_4_or_5(Real &theta_lower, Real &theta_upper) {
  theta_lower = Real(constants<Real>::PI) / Real(2.);
  theta_upper = Real(3.) * constants<Real>::PI / Real(4.);
}
/**
 * @brief Computes the lower and upper bounds for theta in region 6 or 7.
 * @ingroup ETBounds
 * @tparam Real Floating point type.
 * @param[out] theta_lower Lower bound of theta.
 * @param[out] theta_upper Upper bound of theta.
 */
template <class Real>
void ET_theta_bounds_6_or_7(Real &theta_lower, Real &theta_upper) {
  theta_lower = Real(3.) * constants<Real>::PI / Real(4.);
  theta_upper = constants<Real>::PI;
}

// Bounds functions for psi
/**
 * @brief Computes the psi bounds for region 0.
 * @ingroup ETBounds
 * Sets the integration interval [0, arctan(cos(theta))].
 *
 * @tparam Real Floating point type.
 * @param atan_cos_theta Precomputed arctangent of cos(theta).
 * @param unused Unused parameter.
 * @param[out] psi_lower Lower bound of psi.
 * @param[out] psi_upper Upper bound of psi.
 */
template <class Real>
void ET_psi_bounds_0(const Real &atan_cos_theta, const Real &, Real &psi_lower,
                     Real &psi_upper) {
  psi_lower = Real(0.);
  psi_upper = atan_cos_theta;
}
/**
 * @brief Computes the psi bounds for region 1.
 * @ingroup ETBounds
 * Sets the integration interval [arctan(cos(theta)), pi/2].
 *
 * @tparam Real Floating point type.
 * @param atan_cos_theta Precomputed arctangent of cos(theta).
 * @param unused Unused parameter.
 * @param[out] psi_lower Lower bound of psi.
 * @param[out] psi_upper Upper bound of psi.
 */
template <class Real>
void ET_psi_bounds_1(const Real &atan_cos_theta, const Real &, Real &psi_lower,
                     Real &psi_upper) {
  psi_lower = atan_cos_theta;
  psi_upper = constants<Real>::PI / Real(2.);
}
/**
 * @brief Computes the psi bounds for region 2.
 * @ingroup ETBounds
 * Sets the integration interval [0, arctan(sin(theta))].
 *
 * @tparam Real Floating point type.
 * @param unused Unused parameter.
 * @param atan_sin_theta Precomputed arctangent of sin(theta).
 * @param[out] psi_lower Lower bound of psi.
 * @param[out] psi_upper Upper bound of psi.
 */
template <class Real>
void ET_psi_bounds_2(const Real &, const Real &atan_sin_theta, Real &psi_lower,
                     Real &psi_upper) {
  psi_lower = Real(0.);
  psi_upper = atan_sin_theta;
}
/**
 * @brief Computes the psi bounds for region 3.
 * @ingroup ETBounds
 * Sets the integration interval [arctan(sin(theta)), pi/2].
 *
 * @tparam Real Floating point type.
 * @param unused Unused parameter.
 * @param atan_sin_theta Precomputed arctangent of sin(theta).
 * @param[out] psi_lower Lower bound of psi.
 * @param[out] psi_upper Upper bound of psi.
 */
template <class Real>
void ET_psi_bounds_3(const Real &, const Real &atan_sin_theta, Real &psi_lower,
                     Real &psi_upper) {
  psi_lower = atan_sin_theta;
  psi_upper = constants<Real>::PI / Real(2.);
}
/**
 * @brief Computes the psi bounds for region 4.
 * @ingroup ETBounds
 * Sets the integration interval [0, arctan(sin(theta))].
 *
 * @tparam Real Floating point type.
 * @param unused Unused parameter.
 * @param atan_sin_theta Precomputed arctangent of sin(theta).
 * @param[out] psi_lower Lower bound of psi.
 * @param[out] psi_upper Upper bound of psi.
 */
template <class Real>
void ET_psi_bounds_4(const Real &, const Real &atan_sin_theta, Real &psi_lower,
                     Real &psi_upper) {
  psi_lower = Real(0.);
  psi_upper = atan_sin_theta;
}
/**
 * @brief Computes the psi bounds for region 5.
 * @ingroup ETBounds
 * Sets the integration interval [arctan(sin(theta)), pi/2].
 *
 * @tparam Real Floating point type.
 * @param unused Unused parameter.
 * @param atan_sin_theta Precomputed arctangent of sin(theta).
 * @param[out] psi_lower Lower bound of psi.
 * @param[out] psi_upper Upper bound of psi.
 */
template <class Real>
void ET_psi_bounds_5(const Real &, const Real &atan_sin_theta, Real &psi_lower,
                     Real &psi_upper) {
  psi_lower = atan_sin_theta;
  psi_upper = constants<Real>::PI / Real(2.);
}
/**
 * @brief Computes the psi bounds for region 6.
 * @ingroup ETBounds
 * Sets the integration interval [0, -arctan(cos(theta))].
 *
 * @tparam Real Floating point type.
 * @param atan_cos_theta Precomputed arctangent of cos(theta).
 * @param unused Unused parameter.
 * @param[out] psi_lower Lower bound of psi.
 * @param[out] psi_upper Upper bound of psi.
 */
template <class Real>
void ET_psi_bounds_6(const Real &atan_cos_theta, const Real &, Real &psi_lower,
                     Real &psi_upper) {
  psi_lower = Real(0.);
  psi_upper = -atan_cos_theta;
}
/**
 * @brief Computes the psi bounds for region 7.
 * @ingroup ETBounds
 * Sets the integration interval [-arctan(cos(theta)), pi/2].
 *
 * @tparam Real Floating point type.
 * @param atan_cos_theta Precomputed arctangent of cos(theta).
 * @param unused Unused parameter.
 * @param[out] psi_lower Lower bound of psi.
 * @param[out] psi_upper Upper bound of psi.
 */
template <class Real>
void ET_psi_bounds_7(const Real &atan_cos_theta, const Real &, Real &psi_lower,
                     Real &psi_upper) {
  psi_lower = -atan_cos_theta;
  psi_upper = constants<Real>::PI / Real(2.);
}

// For the u bounds
/**
 * @brief Computes the full interval u bounds for regions 0 or 6.
 * @ingroup ETBounds
 * Sets u in [-1, 1].
 *
 * @tparam Real Floating point type.
 * @param unused1 Unused parameter.
 * @param unused2 Unused parameter.
 * @param unused3 Unused parameter.
 * @param[out] u_lower Lower bound of u.
 * @param[out] u_upper Upper bound of u.
 */
template <class Real>
void ET_u_bounds_0_or_6(const Real &, const Real &, const Real&,
                        Real &u_lower, Real &u_upper) {
  u_lower = Real(-1.);
  u_upper = Real(1.);
}
/**
 * @brief Computes the lower interval u bounds for region 1.
 * @ingroup ETBounds
 * @tparam Real Floating point type.
 * @param cos_theta Cosine of theta.
 * @param unused Unused parameter.
 * @param tan_psi Tangent of psi.
 * @param[out] u_lower Lower bound of u.
 * @param[out] u_upper Upper bound of u.
 */

template <class Real>
void ET_u_bounds_1_lower(const Real &cos_theta, const Real &,
                         const Real &tan_psi, Real &u_lower, Real &u_upper) {
  u_lower = Real(-1.);
  u_upper = Real(2.) * cos_theta / tan_psi - Real(1.);
}
/**
 * @brief Computes the lower interval u bounds for region 2.
 * @ingroup ETBounds
 * @tparam Real Floating point type.
 * @param unused1 Unused parameter.
 * @param tan_theta Tangent of theta.
 * @param unused2 Unused parameter.
 * @param[out] u_lower Lower bound of u.
 * @param[out] u_upper Upper bound of u.
 */

template <class Real>
void ET_u_bounds_2_lower(const Real &, const Real &tan_theta,
                         const Real &, Real &u_lower, Real &u_upper) {
  u_lower = Real(-1.);
  u_upper = Real(2.0) / tan_theta - Real(1.0);
}
/**
 * @brief Computes the lower interval u bounds for region 3.
 * @ingroup ETBounds
 * @tparam Real Floating point type.
 * @param cos_theta Cosine of theta.
 * @param unused Unused parameter.
 * @param tan_psi Tangent of psi.
 * @param[out] u_lower Lower bound of u.
 * @param[out] u_upper Upper bound of u.
 */

template <class Real>
void ET_u_bounds_3_lower(const Real &cos_theta, const Real &,
                         const Real &tan_psi, Real &u_lower, Real &u_upper) {
  u_lower = Real(-1.);
  u_upper = Real(2.) * cos_theta / tan_psi - Real(1.); // u_1_theta
}
/**
 * @brief Computes the lower interval u bounds for region 4.
 * @ingroup ETBounds
 * @tparam Real Floating point type.
 * @param unused1 Unused parameter.
 * @param tan_theta Tangent of theta.
 * @param unused2 Unused parameter.
 * @param[out] u_lower Lower bound of u.
 * @param[out] u_upper Upper bound of u.
 */

template <class Real>
void ET_u_bounds_4_lower(const Real &, const Real &tan_theta,
                         const Real &, Real &u_lower, Real &u_upper) {
  u_lower = Real(-1.);
  u_upper = Real(2.) / tan_theta + Real(1.); // u_2_theta;
}
/**
 * @brief Computes the lower interval u bounds for regions 5 or 7.
 * @ingroup ETBounds
 * @tparam Real Floating point type.
 * @param cos_theta Cosine of theta.
 * @param unused Unused parameter.
 * @param tan_psi Tangent of psi.
 * @param[out] u_lower Lower bound of u.
 * @param[out] u_upper Upper bound of u.
 */

template <class Real>
void ET_u_bounds_5_or_7_lower(const Real &cos_theta, const Real &,
                              const Real &tan_psi, Real &u_lower,
                              Real &u_upper) {
  u_lower = Real(-1.);
  u_upper = Real(2.) * cos_theta / tan_psi + Real(1.); // u_2_psi;
}
/**
 * @brief Computes the upper interval u bounds for region 1.
 * @ingroup ETBounds
 * @tparam Real Floating point type.
 * @param cos_theta Cosine of theta.
 * @param unused Unused parameter.
 * @param tan_psi Tangent of psi.
 * @param[out] u_lower Lower bound of u.
 * @param[out] u_upper Upper bound of u.
 */

template <class Real>
void ET_u_bounds_1_upper(const Real &cos_theta, const Real &,
                         const Real &tan_psi, Real &u_lower, Real &u_upper) {
  u_lower = Real(2.) * cos_theta / tan_psi - Real(1.);
  u_upper = Real(1.);
}
/**
 * @brief Computes the upper interval u bounds for region 2.
 * @ingroup ETBounds
 * @tparam Real Floating point type.
 * @param unused1 Unused parameter.
 * @param tan_theta Tangent of theta.
 * @param unused2 Unused parameter.
 * @param[out] u_lower Lower bound of u.
 * @param[out] u_upper Upper bound of u.
 */

template <class Real>
void ET_u_bounds_2_upper(const Real &, const Real &tan_theta,
                         const Real &, Real &u_lower, Real &u_upper) {
  u_lower = Real(2.) / tan_theta - Real(1.); // u_1_theta;
  u_upper = Real(1.);
}
/**
 * @brief Computes the upper interval u bounds for region 3.
 * @ingroup ETBounds
 * @tparam Real Floating point type.
 * @param cos_theta Cosine of theta.
 * @param unused Unused parameter.
 * @param tan_psi Tangent of psi.
 * @param[out] u_lower Lower bound of u.
 * @param[out] u_upper Upper bound of u.
 */

template <class Real>
void ET_u_bounds_3_upper(const Real &cos_theta, const Real &,
                         const Real &tan_psi, Real &u_lower, Real &u_upper) {
  u_lower = Real(2.) * cos_theta / tan_psi - Real(1.); // u_1_psi;
  u_upper = Real(1.);
}
/**
 * @brief Computes the upper interval u bounds for region 4.
 * @ingroup ETBounds
 * @tparam Real Floating point type.
 * @param unused1 Unused parameter.
 * @param tan_theta Tangent of theta.
 * @param unused2 Unused parameter.
 * @param[out] u_lower Lower bound of u.
 * @param[out] u_upper Upper bound of u.
 */

template <class Real>
void ET_u_bounds_4_upper(const Real &, const Real &tan_theta,
                         const Real &, Real &u_lower, Real &u_upper) {
  u_lower = Real(2.) / tan_theta + Real(1.);
  u_upper = Real(1.);
}
/**
 * @brief Computes the upper interval u bounds for regions 5 or 7.
 * @ingroup ETBounds
 * @tparam Real Floating point type.
 * @param cos_theta Cosine of theta.
 * @param unused Unused parameter.
 * @param tan_psi Tangent of psi.
 * @param[out] u_lower Lower bound of u.
 * @param[out] u_upper Upper bound of u.
 */

template <class Real>
void ET_u_bounds_5_or_7_upper(const Real &cos_theta, const Real &,
                              const Real &tan_psi, Real &u_lower,
                              Real &u_upper) {
  u_lower = Real(2.) * cos_theta / tan_psi + Real(1.);
  u_upper = Real(1.);
}

// For the lambda bounds
/**
 * @brief Computes the lambda bounds for region 0.
 * @ingroup ETBounds
 * Sets lambda in [0, (u + 1)/(cos(theta) * cos(psi))].
 *
 * @tparam Real Floating point type.
 * @param cos_theta Cosine of theta.
 * @param sin_theta Sine of theta.
 * @param cos_psi Cosine of psi.
 * @param sin_psi Sine of psi.
 * @param u Current u parameter.
 * @param[out] lambda_lower Lower bound of lambda.
 * @param[out] lambda_upper Upper bound of lambda.
 */

template <class Real>
void ET_lambda_bounds_0(const Real &cos_theta, const Real &sin_theta,
                                 const Real &cos_psi, const Real &sin_psi,
                                 const Real &u, Real &lambda_lower,
                                 Real &lambda_upper) {
  lambda_lower = Real(0.);
  lambda_upper = (u + Real(1.)) / (cos_theta * cos_psi);
}
/**
 * @brief Computes the lambda bounds for region 6.
 * @ingroup ETBounds
 * Sets lambda in [0, (u - 1)/(cos(theta) * cos(psi))].
 *
 * @tparam Real Floating point type.
 * @param cos_theta Cosine of theta.
 * @param sin_theta Sine of theta.
 * @param cos_psi Cosine of psi.
 * @param sin_psi Sine of psi.
 * @param u Current u parameter.
 * @param[out] lambda_lower Lower bound of lambda.
 * @param[out] lambda_upper Upper bound of lambda.
 */

template <class Real>
void ET_lambda_bounds_6(const Real &cos_theta, const Real &sin_theta,
                        const Real &cos_psi, const Real &sin_psi,
                        const Real &u, Real &lambda_lower,
                        Real &lambda_upper) {
  lambda_lower = Real(0.);
  lambda_upper =  (u - Real(1.)) / (cos_theta * cos_psi);
}
/**
 * @brief Computes the lambda bounds for regions 1, 2, or 3 in subregion 0.
 * @ingroup ETBounds
 * Sets lambda in [0, (u + 1)/(cos(theta) * cos(psi))].
 *
 * @tparam Real Floating point type.
 * @param cos_theta Cosine of theta.
 * @param sin_theta Sine of theta.
 * @param cos_psi Cosine of psi.
 * @param sin_psi Sine of psi.
 * @param u Current u parameter.
 * @param[out] lambda_lower Lower bound of lambda.
 * @param[out] lambda_upper Upper bound of lambda.
 */

template <class Real>
void ET_lambda_bounds_1_2_or_3_region_0(const Real &cos_theta, const Real &sin_theta,
                        const Real &cos_psi, const Real &sin_psi,
                        const Real &u, Real &lambda_lower,
                        Real &lambda_upper) {
  lambda_lower = Real(0.);
  lambda_upper =  (u + Real(1.0)) / (cos_theta * cos_psi);
}
/**
 * @brief Computes the lambda bounds for region 4 in subregion 0.
 * @ingroup ETBounds
 * Sets lambda in [0, 2/(sin(theta) * cos(psi))].
 *
 * @tparam Real Floating point type.
 * @param cos_theta Cosine of theta.
 * @param sin_theta Sine of theta.
 * @param cos_psi Cosine of psi.
 * @param sin_psi Sine of psi.
 * @param u Current u parameter.
 * @param[out] lambda_lower Lower bound of lambda.
 * @param[out] lambda_upper Upper bound of lambda.
 */

template <class Real>
void ET_lambda_bounds_4_region_0(const Real &cos_theta, const Real &sin_theta,
                                        const Real &cos_psi, const Real &sin_psi,
                                        const Real &u, Real &lambda_lower,
                                        Real &lambda_upper) {
  lambda_lower = Real(0.);
  lambda_upper =  Real(2.) / (sin_theta * cos_psi);
}
/**
 * @brief Computes the lambda bounds for region 5 in subregion 0.
 * @ingroup ETBounds
 * Sets lambda in [0, 2/sin(psi)].
 *
 * @tparam Real Floating point type.
 * @param cos_theta Cosine of theta.
 * @param sin_theta Sine of theta.
 * @param cos_psi Cosine of psi.
 * @param sin_psi Sine of psi.
 * @param u Current u parameter.
 * @param[out] lambda_lower Lower bound of lambda.
 * @param[out] lambda_upper Upper bound of lambda.
 */

template <class Real>
void ET_lambda_bounds_5_region_0(const Real &cos_theta, const Real &sin_theta,
                                 const Real &cos_psi, const Real &sin_psi,
                                 const Real &u, Real &lambda_lower,
                                 Real &lambda_upper) {
  lambda_lower = Real(0.);
  lambda_upper =  Real(2.) / sin_psi;
}
/**
 * @brief Computes the lambda bounds for region 7 in subregion 0.
 * @ingroup ETBounds
 * Sets lambda in [0, (u - 1)/(cos(theta) * cos(psi))].
 *
 * @tparam Real Floating point type.
 * @param cos_theta Cosine of theta.
 * @param sin_theta Sine of theta.
 * @param cos_psi Cosine of psi.
 * @param sin_psi Sine of psi.
 * @param u Current u parameter.
 * @param[out] lambda_lower Lower bound of lambda.
 * @param[out] lambda_upper Upper bound of lambda.
 */

template <class Real>
void ET_lambda_bounds_7_region_0(const Real &cos_theta, const Real &sin_theta,
                                 const Real &cos_psi, const Real &sin_psi,
                                 const Real &u, Real &lambda_lower,
                                 Real &lambda_upper) {
  lambda_lower = Real(0.);
  lambda_upper =  (u - Real(1.)) / (cos_theta * cos_psi);
}
/**
 * @brief Computes the lambda bounds for regions 1, 3, or 7 in subregion 1.
 * @ingroup ETBounds
 * Sets lambda in [0, 2/sin(psi)].
 *
 * @tparam Real Floating point type.
 * @param cos_theta Cosine of theta.
 * @param sin_theta Sine of theta.
 * @param cos_psi Cosine of psi.
 * @param sin_psi Sine of psi.
 * @param u Current u parameter.
 * @param[out] lambda_lower Lower bound of lambda.
 * @param[out] lambda_upper Upper bound of lambda.
 */

template <class Real>
void ET_lambda_bounds_1_3_or_7_region_1(const Real &cos_theta, const Real &sin_theta,
                                        const Real &cos_psi, const Real &sin_psi,
                                        const Real &u, Real &lambda_lower,
                                        Real &lambda_upper) {
  lambda_lower = Real(0.);
  lambda_upper =  Real(2.) / sin_psi;
}
/**
 * @brief Computes the lambda bounds for region 2 in subregion 1.
 * @ingroup ETBounds
 * Sets lambda in [0, 2/(sin(theta) * cos(psi))].
 *
 * @tparam Real Floating point type.
 * @param cos_theta Cosine of theta.
 * @param sin_theta Sine of theta.
 * @param cos_psi Cosine of psi.
 * @param sin_psi Sine of psi.
 * @param u Current u parameter.
 * @param[out] lambda_lower Lower bound of lambda.
 * @param[out] lambda_upper Upper bound of lambda.
 */

template <class Real>
void ET_lambda_bounds_2_region_1(const Real &cos_theta, const Real &sin_theta,
                                 const Real &cos_psi, const Real &sin_psi,
                                 const Real &u, Real &lambda_lower,
                                 Real &lambda_upper) {
  lambda_lower = Real(0.);
  lambda_upper =  Real(2.) / (sin_theta * cos_psi);
}
/**
 * @brief Computes the lambda bounds for regions 4 or 5 in subregion 1.
 * @ingroup ETBounds
 * Sets lambda in [0, (u - 1)/(cos(theta) * cos(psi))].
 *
 * @tparam Real Floating point type.
 * @param cos_theta Cosine of theta.
 * @param sin_theta Sine of theta.
 * @param cos_psi Cosine of psi.
 * @param sin_psi Sine of psi.
 * @param u Current u parameter.
 * @param[out] lambda_lower Lower bound of lambda.
 * @param[out] lambda_upper Upper bound of lambda.
 */

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

/**
 * @brief Leaves the u and v coordinates unchanged.
 *
 * This is used when no reorientation of the parametric coordinates is needed.
 *
 * @ingroup ETMapping
 *
 * @tparam Real Floating point type.
 * @param[in,out] uv The 2D coordinate vector.
 */

template <class Real>
void ET_map_u_and_v_coords_0(Point<2, Real>&)
{
  return;
}
/**
 * @brief Maps u and v coordinates for region 1.
 *
 * Rotates the coordinates by 90 degrees counterclockwise and negates u to align element orientation.
 *
 * @ingroup ETMapping
 *
 * @tparam Real Floating point type.
 * @param[in,out] uv The 2D coordinate vector to transform.
 */

template <class Real>
void ET_map_u_and_v_coords_1(Point<2, Real>& uv)
{
  Real temp = uv(0);
  uv[0] = -uv(1);
  uv[1] = temp;
}
/**
 * @brief Maps u and v coordinates for region 2.
 *
 * Rotates the coordinates by 90 degrees clockwise and negates v to align element orientation.
 *
 * @ingroup ETMapping
 *
 * @tparam Real Floating point type.
 * @param[in,out] uv The 2D coordinate vector to transform.
 */

template <class Real>
void ET_map_u_and_v_coords_2(Point<2, Real>& uv)
{
  Real temp = uv(0);
  uv[0] = uv(1);
  uv[1] = -temp;
}
/**
 * @brief Maps u and v coordinates for region 3.
 *
 * Rotates the coordinates by 180 degrees (negates both components) to align element orientation.
 *
 * @ingroup ETMapping
 *
 * @tparam Real Floating point type.
 * @param[in,out] uv The 2D coordinate vector to transform.
 */

template <class Real>
void ET_map_u_and_v_coords_3(Point<2, Real>& uv)
{
  uv[0] *= -1;
  uv[1] *= -1;
}

DROMON_NAMESPACE_CLOSE
#endif // DROMON_DIRECTFN_ET_BOUNDS_FUNCTIONS_H
