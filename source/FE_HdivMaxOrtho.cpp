//
// Created by Jake J. Harmon (jake.harmon@ieee.org) on 8/6/21.
//
#include "FE_HdivMaxOrtho.h"

template <>
double dromon::FE_HdivMaxOrtho<2, 3, double>::evaluate_shape_function(
    const DoFDirection &direction, const MultiIndex<2> &indices,
    const Point<2, double> &uv) {
  Point<2, double> uv_temp = uv;
  MultiIndex<2> indices_temp = indices;
  double out;

  // f_udir_ij = Q_i(u)*P_j(v),
  // where Q_i is the max-ortho basis function as in $[max_ortho_paper]$ and
  // P_j denotes the Legendre polynomial of degree j

  // Assumes that the direction is u-directed, then adjusts later
  if (direction != DoFDirection::u_dir) {
    uv_temp[1] = uv(0);
    uv_temp[0] = uv(1);

    indices_temp[1] = indices(0);
    indices_temp[0] = indices(1);
  }
  // First, collect the P_j(v) portion...
  out = boost::math::legendre_p(indices_temp(1), uv_temp(1));

  // Now compute the Q function portion
  switch (indices_temp(0)) {
  case 0:
    out *= (double)(1) - uv_temp(0);
    break;
  case 1:
    out *= (double)(1) + uv_temp(0);
    break;
  case 2:
//    out *= boost::math::legendre_p(2, uv_temp(0)) -
//           boost::math::legendre_p(0, uv_temp(0));
      out *= utility::legendre_q_2(uv_temp(0));
    break;
  case 3:
    out *= boost::math::legendre_p(3, uv_temp(0)) -
           boost::math::legendre_p(1, uv_temp(0));
    break;
  case 4:
    out *= ((double)(6) * boost::math::legendre_p(4, uv_temp(0)) -
            (double)(5) * boost::math::legendre_p(2, uv_temp(0)) -
            boost::math::legendre_p(0, uv_temp(0))) /
           ((double)(6));
    break;
  case 5:
    out *= ((double)(10) * boost::math::legendre_p(5, uv_temp(0)) -
            (double)(7) * boost::math::legendre_p(3, uv_temp(0)) -
            (double)(3) * boost::math::legendre_p(1, uv_temp(0))) /
           ((double)(10));
    break;
  case 6:
    out *= ((double)(15) * boost::math::legendre_p(6, uv_temp(0)) -
            (double)(9) * boost::math::legendre_p(4, uv_temp(0)) -
            (double)(5) * boost::math::legendre_p(2, uv_temp(0)) -
            boost::math::legendre_p(0, uv(0))) /
           ((double)(15));
    break;
  case 7:
    out *= ((double)(21) * boost::math::legendre_p(7, uv_temp(0)) -
            (double)(11) * boost::math::legendre_p(5, uv_temp(0)) -
            (double)(7) * boost::math::legendre_p(3, uv_temp(0)) -
            (double)(3) * boost::math::legendre_p(1, uv_temp(0))) /
           ((double)(21));
    break;
  case 8:
    out *= ((double)(28) * boost::math::legendre_p(8, uv_temp(0)) -
            (double)(13) * boost::math::legendre_p(6, uv_temp(0)) -
            (double)(9) * boost::math::legendre_p(4, uv_temp(0)) -
            (double)(5) * boost::math::legendre_p(2, uv_temp(0)) -
            boost::math::legendre_p(0, uv_temp(0))) /
           ((double)(28));
    break;
    // Found using the obvious "pattern" in the construction of these basis
    // functions
  case 9:
    out *= ((double)(36) * boost::math::legendre_p(9, uv_temp(0)) -
            (double)(15) * boost::math::legendre_p(7, uv_temp(0)) -
            (double)(11) * boost::math::legendre_p(5, uv_temp(0)) -
            (double)(7) * boost::math::legendre_p(3, uv_temp(0)) -
            (double)(3) * boost::math::legendre_p(1, uv_temp(0))) /
           ((double)(36));
    break;
  case 10:
    out *= ((double)(45) * boost::math::legendre_p(10, uv_temp(0)) -
            (double)(17) * boost::math::legendre_p(8, uv_temp(0)) -
            (double)(13) * boost::math::legendre_p(6, uv_temp(0)) -
            (double)(9) * boost::math::legendre_p(4, uv_temp(0)) -
            (double)(5) * boost::math::legendre_p(2, uv_temp(0)) -
            boost::math::legendre_p(0, uv_temp(0))) /
           ((double)(45));
    break;
  default:
    assert(false && "Case not supported (yet) for FE_HdivMaxOrtho!");
    break;
  }

  return euclidean_scaling[indices_temp(0)]*out;
}
#include "FE_HdivMaxOrtho.h"
