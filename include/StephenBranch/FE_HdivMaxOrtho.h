//
// Created by Jake J. Harmon (jake.harmon@ieee.org) on 8/6/21.
//

#ifndef DROMON_FE_HDIVMAXORTHO_H
#define DROMON_FE_HDIVMAXORTHO_H

#include "DoFBase.h"
#include "FEBase.h"
#include "MultiIndex.h"
#include "Point.h"
#include "config.h"
#include <complex>
#include "LegendrePFast.h"

// For evaluating the legendre polynomials and their deriviates
#include <boost/math/special_functions/legendre.hpp>

DROMON_NAMESPACE_OPEN

// Create a local index of all the potential degrees of freedom?
template <unsigned int dim, unsigned int spacedim, class Real = double>
class FE_HdivMaxOrtho : public FiniteElement<dim, Real> {
public:
  FE_HdivMaxOrtho(const unsigned int &degree, const CurrentType &type);

  // virtual Real evaluate_shape_function(const DoFDirection& direction, const
  // MultiIndex<dim>& indices, const Point<dim, Real>& uv) override;
  virtual Real evaluate_shape_function(const DoFDirection &direction,
                                       const MultiIndex<dim> &indices,
                                       const Point<dim, Real> &uv) override;
  virtual Real evaluate_shape_function(const DoFBase<dim> &dof,
                                       const Point<dim, Real> &uv) override;
  virtual Real
  evaluate_shape_function_divergence(const DoFDirection &direction,
                                     const MultiIndex<dim> &indices,
                                     const Point<dim, Real> &uv) override;
  virtual Real
  evaluate_shape_function_divergence(const DoFBase<dim> &dof,
                                     const Point<dim, Real> &uv) override;

  std::string get_FE_name() const override;

  virtual bool
  is_subset_of_max(const MultiIndex<dim> &indices, DoFDirection dir,
                   const unsigned int max_expasion_order) const override;

  virtual bool contains_this_dof(const DoFBase<dim> &dof) const override;

  // virtual std::pair<MultiIndex<dim, double>, bool> determine_pairing(const
  // unsigned int& local_face_index, const unsigned int & neighbor_face_index,
  // const bool& matched_face_orientation) const override;
  virtual void
  determine_pairing(const unsigned int &local_face_index,
                    const unsigned int &neighbor_face_index,
                    const bool &matched_face_orientation,
                    MultiIndex<dim, double> *coordinate_orientation,
                    bool *invert_value) const override;
  // Make part of the DoFHandler...
  //    template <class CoefficientType>
  //    CoefficientType evaluate_shape_function(const unsigned int direction,
  //    const MultiIndex<dim>& indices, const Point<dim, Real>& uv, const
  //    CoefficientType& coeff);
private:
  // Add memoization functions here...
  utility::LegendrePContainer<Real> legendrePContainer;

  //  const Real euclidean_scaling[11] = {0.6123724356957945245493210186,
  //                                      0.6123724356957945245493210186,
  //                                      0.6454972243679028141965442332,
  //                                      1.0246950765959598383221038680,
  //                                      1.3416407864998738178455042012,
  //                                      1.0,1.0,1.0,1.0,1.0,1.0};
  const Real euclidean_scaling[11] = {1.0,1.0,1.0,1.0,1.0,
                                      1.0,1.0,1.0,1.0,1.0,1.0};
};

template <unsigned int dim, unsigned int spacedim, class Real>
FE_HdivMaxOrtho<dim, spacedim, Real>::FE_HdivMaxOrtho(
    const unsigned int &degree, const CurrentType &type)
    : FiniteElement<dim, Real>(degree, dim, Hdiv, type) {}


// While memoization is not really possible in the general case (for the
// numerical integration via the DIRECTFN method), the regular integration
// can easily be memoized for a given quadrature...
template <unsigned int dim, unsigned int spacedim, class Real>
Real FE_HdivMaxOrtho<dim, spacedim, Real>::evaluate_shape_function(
    const DoFDirection &direction, const MultiIndex<dim> &indices,
    const Point<dim, Real> &uv) {
  Point<dim, Real> uv_temp = uv;
  MultiIndex<dim> indices_temp = indices;
  Real out;
  switch (dim) {
  case 2:
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
    // out = legendrePContainer(indices_temp(1), uv_temp(1));
    // Now compute the Q function portion
    switch (indices_temp(0)) {
    case 0:
      out *= (Real)(1) - uv_temp(0);
      break;
    case 1:
      out *= (Real)(1) + uv_temp(0);
      break;
    case 2:
      out *= boost::math::legendre_p(2, uv_temp(0)) -
             boost::math::legendre_p(0, uv_temp(0));
      // out *= utility::legendre_q_2(uv_temp(0));
      break;
    case 3:
      out *= boost::math::legendre_p(3, uv_temp(0)) -
             boost::math::legendre_p(1, uv_temp(0));
      break;
    case 4:
      out *= ((Real)(6) * boost::math::legendre_p(4, uv_temp(0)) -
              (Real)(5) * boost::math::legendre_p(2, uv_temp(0)) -
              boost::math::legendre_p(0, uv_temp(0))) /
             ((Real)(6));
      break;
    case 5:
      out *= ((Real)(10) * boost::math::legendre_p(5, uv_temp(0)) -
              (Real)(7) * boost::math::legendre_p(3, uv_temp(0)) -
              (Real)(3) * boost::math::legendre_p(1, uv_temp(0))) /
             ((Real)(10));
      break;
    case 6:
      out *= ((Real)(15) * boost::math::legendre_p(6, uv_temp(0)) -
              (Real)(9) * boost::math::legendre_p(4, uv_temp(0)) -
              (Real)(5) * boost::math::legendre_p(2, uv_temp(0)) -
              boost::math::legendre_p(0, uv(0))) /
             ((Real)(15));
      break;
    case 7:
      out *= ((Real)(21) * boost::math::legendre_p(7, uv_temp(0)) -
              (Real)(11) * boost::math::legendre_p(5, uv_temp(0)) -
              (Real)(7) * boost::math::legendre_p(3, uv_temp(0)) -
              (Real)(3) * boost::math::legendre_p(1, uv_temp(0))) /
             ((Real)(21));
      break;
    case 8:
      out *= ((Real)(28) * boost::math::legendre_p(8, uv_temp(0)) -
              (Real)(13) * boost::math::legendre_p(6, uv_temp(0)) -
              (Real)(9) * boost::math::legendre_p(4, uv_temp(0)) -
              (Real)(5) * boost::math::legendre_p(2, uv_temp(0)) -
              boost::math::legendre_p(0, uv_temp(0))) /
             ((Real)(28));
      break;
      // Found using the obvious "pattern" in the construction of these basis
      // functions
    case 9:
      out *= ((Real)(36) * boost::math::legendre_p(9, uv_temp(0)) -
              (Real)(15) * boost::math::legendre_p(7, uv_temp(0)) -
              (Real)(11) * boost::math::legendre_p(5, uv_temp(0)) -
              (Real)(7) * boost::math::legendre_p(3, uv_temp(0)) -
              (Real)(3) * boost::math::legendre_p(1, uv_temp(0))) /
             ((Real)(36));
      break;
    case 10:
      out *= ((Real)(45) * boost::math::legendre_p(10, uv_temp(0)) -
              (Real)(17) * boost::math::legendre_p(8, uv_temp(0)) -
              (Real)(13) * boost::math::legendre_p(6, uv_temp(0)) -
              (Real)(9) * boost::math::legendre_p(4, uv_temp(0)) -
              (Real)(5) * boost::math::legendre_p(2, uv_temp(0)) -
              boost::math::legendre_p(0, uv_temp(0))) /
             ((Real)(45));
      break;
    default:
      assert(false && "Case not supported (yet) for FE_HdivMaxOrtho!");
      break;
    }

    break;
  default:
    assert(false && "Case not supported (yet) for FE_HdivMaxOrtho!");
    break;
  }
  out *= euclidean_scaling[indices_temp(0)];
  return out;
}



template <unsigned int dim, unsigned int spacedim, class Real>
std::string FE_HdivMaxOrtho<dim, spacedim, Real>::get_FE_name() const {
  return "FE_HdivMaxOrtho";
}

template <unsigned int dim, unsigned int spacedim, class Real>
bool FE_HdivMaxOrtho<dim, spacedim, Real>::is_subset_of_max(
    const MultiIndex<dim> &indices, DoFDirection dir,
    const unsigned int max_expasion_order) const {
  switch (dir) {
  case (u_dir):
    for (unsigned int i = 0; i < dim; ++i) {
      if (i == 0) {
        if (indices(i) > max_expasion_order)
          return false;
      } else if (indices(i) >= max_expasion_order)
        return false;
    }
    return true;
    break;
  case (v_dir):
    for (unsigned int i = 0; i < dim; ++i) {
      if (i == 1) {
        if (indices(i) > max_expasion_order)
          return false;
      } else if (indices(i) >= max_expasion_order)
        return false;
    }
    return true;
    break;
  case (w_dir):
    assert(false && "Not implemented yet!");
    return false;
    break;
  default:
    assert(false && "Incorrect direction designation!");
    return false;
    break;
  }
}
template <unsigned int dim, unsigned int spacedim, class Real>
bool FE_HdivMaxOrtho<dim, spacedim, Real>::contains_this_dof(
    const DoFBase<dim> &dof) const {
  if (this->get_type() != dof.cur_type)
    return false;
  if (dof.direction == u_dir)
    return (dof.orders(0) <= this->get_degree() &&
            dof.orders(1) < this->get_degree());
  else if (dof.direction == v_dir)
    return (dof.orders(0) < this->get_degree() &&
            dof.orders(1) <= this->get_degree());
  else
    assert(false && "Direction not implemented!");
  return false;
}

template <unsigned int dim, unsigned int spacedim, class Real>
void FE_HdivMaxOrtho<dim, spacedim, Real>::determine_pairing(
    const unsigned int &local_face_index,
    const unsigned int &neighbor_face_index,
    const bool &matched_face_orientation,
    MultiIndex<dim, double> *coordinate_orientation, bool *invert_value) const {
  bool invert_sign = (local_face_index % 2 == neighbor_face_index % 2);
  MultiIndex<dim, double> invert_coordinates(1.0);
  if (!matched_face_orientation) {
    switch (dim) {
    case 2:
      // Depending on the local face index, we must either invert
      // the u_dir or the v_dir
      if (local_face_index == 0 || local_face_index == 3)
        invert_coordinates[0] = -1.0;
      else
        invert_coordinates[1] = -1.0;
      break;
    default:
      assert(false && "Dimension not implemented!");
    }
  }
  *coordinate_orientation = invert_coordinates;
  *invert_value = invert_sign;
}

/**
 * Since the shape functions defined on the reference cell are unidirectional
 * (Note: after the Piola transformation, on the real cell, they are not
 * necessarily unidirectional), we can take just a partial directive with
 * respect to the direction of the shape function to arrive at the divergence
 * @tparam dim
 * @tparam spacedim
 * @tparam Real
 * @param direction
 * @param indices
 * @param uv
 * @return
 */
template <unsigned int dim, unsigned int spacedim, class Real>
Real FE_HdivMaxOrtho<dim, spacedim, Real>::evaluate_shape_function_divergence(
    const DoFDirection &direction, const MultiIndex<dim> &indices,
    const Point<dim, Real> &uv) {
  Point<dim, Real> uv_temp = uv;
  MultiIndex<dim> indices_temp = indices;
  Real out;
  switch (dim) {
  case 2:
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
      out *= -(Real)(1);
      break;
    case 1:
      out *= (Real)(1);
      break;
    case 2:
      out *= boost::math::legendre_p_prime(2, uv_temp(0)) -
             boost::math::legendre_p_prime(0, uv_temp(0));
      break;
    case 3:
      out *= boost::math::legendre_p_prime(3, uv_temp(0)) -
             boost::math::legendre_p_prime(1, uv_temp(0));
      break;
    case 4:
      out *= ((Real)(6) * boost::math::legendre_p_prime(4, uv_temp(0)) -
              (Real)(5) * boost::math::legendre_p_prime(2, uv_temp(0)) -
              boost::math::legendre_p_prime(0, uv_temp(0))) /
             ((Real)(6));
      break;
    case 5:
      out *= ((Real)(10) * boost::math::legendre_p_prime(5, uv_temp(0)) -
              (Real)(7) * boost::math::legendre_p_prime(3, uv_temp(0)) -
              (Real)(3) * boost::math::legendre_p_prime(1, uv_temp(0))) /
             ((Real)(10));
      break;
    case 6:
      out *= ((Real)(15) * boost::math::legendre_p_prime(6, uv_temp(0)) -
              (Real)(9) * boost::math::legendre_p_prime(4, uv_temp(0)) -
              (Real)(5) * boost::math::legendre_p_prime(2, uv_temp(0)) -
              boost::math::legendre_p_prime(0, uv(0))) /
             ((Real)(15));
      break;
    case 7:
      out *= ((Real)(21) * boost::math::legendre_p_prime(7, uv_temp(0)) -
              (Real)(11) * boost::math::legendre_p_prime(5, uv_temp(0)) -
              (Real)(7) * boost::math::legendre_p_prime(3, uv_temp(0)) -
              (Real)(3) * boost::math::legendre_p_prime(1, uv_temp(0))) /
             ((Real)(21));
      break;
    case 8:
      out *= ((Real)(28) * boost::math::legendre_p_prime(8, uv_temp(0)) -
              (Real)(13) * boost::math::legendre_p_prime(6, uv_temp(0)) -
              (Real)(9) * boost::math::legendre_p_prime(4, uv_temp(0)) -
              (Real)(5) * boost::math::legendre_p_prime(2, uv_temp(0)) -
              boost::math::legendre_p_prime(0, uv_temp(0))) /
             ((Real)(28));
      break;
      // Found using the obvious "pattern" in the construction of these basis
      // functions
    case 9:
      out *= ((Real)(36) * boost::math::legendre_p_prime(9, uv_temp(0)) -
              (Real)(15) * boost::math::legendre_p_prime(7, uv_temp(0)) -
              (Real)(11) * boost::math::legendre_p_prime(5, uv_temp(0)) -
              (Real)(7) * boost::math::legendre_p_prime(3, uv_temp(0)) -
              (Real)(3) * boost::math::legendre_p_prime(1, uv_temp(0))) /
             ((Real)(36));
      break;
    case 10:
      out *= ((Real)(45) * boost::math::legendre_p_prime(10, uv_temp(0)) -
              (Real)(17) * boost::math::legendre_p_prime(8, uv_temp(0)) -
              (Real)(13) * boost::math::legendre_p_prime(6, uv_temp(0)) -
              (Real)(9) * boost::math::legendre_p_prime(4, uv_temp(0)) -
              (Real)(5) * boost::math::legendre_p_prime(2, uv_temp(0)) -
              boost::math::legendre_p_prime(0, uv_temp(0))) /
             ((Real)(45));
      break;
    default:
      assert(false && "Case not supported (yet) for FE_HdivMaxOrtho!");
      break;
    }

    break;
  default:
    assert(false && "Case not supported (yet) for FE_HdivMaxOrtho!");
    break;
  }
  out *= euclidean_scaling[indices_temp(0)];
  return out;
}
template <unsigned int dim, unsigned int spacedim, class Real>
Real FE_HdivMaxOrtho<dim, spacedim, Real>::evaluate_shape_function(
    const DoFBase<dim> &dof, const Point<dim, Real> &uv) {
  // Use the information contained in the DoF to call the other
  // evaluate_shape_function call...
  Point<dim, Real> point = uv;
  for (unsigned int i = 0; i < dim; ++i) {
    point[i] *= dof.coordinate_orientation.at(i);
  }
  if (dof.invert_value == true)
    return Real(-1) *
           evaluate_shape_function(dof.direction, dof.orders, point);
  else
    return evaluate_shape_function(dof.direction, dof.orders, point);
}
template <unsigned int dim, unsigned int spacedim, class Real>
Real FE_HdivMaxOrtho<dim, spacedim, Real>::evaluate_shape_function_divergence(
    const DoFBase<dim> &dof, const Point<dim, Real> &uv) {
  // Use the information contained in the DoF to call the other
  // evaluate_shape_function_divergence call...
  Point<dim, Real> point = uv;
  for (unsigned int i = 0; i < dim; ++i) {
    point[i] *= dof.coordinate_orientation.at(i);
  }
  if (dof.invert_value == true)
    return Real(-1) *
           evaluate_shape_function_divergence(dof.direction, dof.orders, point);
  else
    return evaluate_shape_function_divergence(dof.direction, dof.orders, point);
}

DROMON_NAMESPACE_CLOSE
#endif // DROMON_FE_HDIVMAXORTHO_H