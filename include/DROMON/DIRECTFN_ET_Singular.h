//
// Created by Jake J. Harmon (jake.harmon@ieee.org) on 4/30/22.
//

#ifndef DROMON_DIRECTFN_ET_SINGULAR_H
#define DROMON_DIRECTFN_ET_SINGULAR_H
#include "DIRECTFN_ET_Bounds_Functions.h"
#include "DoFBase.h"
#include "DoFGeom.h"
#include "DoFMask.h"
#include "Kernels.h"
#include "QuadratureCollection.h"
#include "SubMatrix.h"
#include "config.h"
#include <iomanip>

DROMON_NAMESPACE_OPEN
namespace NumericalIntegration {
template <class DoFCellType, class CoefficientType = std::complex<double>,
          class Real = double>
class DIRECTFN_ET2 {
public:
  DIRECTFN_ET2() {
    assert(false && "The general implementation of the DIRECTFN_ST numerical "
                    "integrator is not yet implemented!");
  }

private:
};

template <class CoefficientType, class Real>
class DIRECTFN_ET2<DoFParent<2, 3, Real>, CoefficientType, Real> {
public:
  /**
   *
   * @param N0 - Order of numerical integration for outermost integral
   * @param N1 - See above, for next integral
   * @param N2 - See above, for next integral
   * @param N3 - See above, for innermost integral
   */
  DIRECTFN_ET2(const unsigned int &N0, const unsigned int &N1,
               const unsigned int &N2, const unsigned int &N3);

  /**
   *
   * @tparam Integrand
   * @param integrand - A lambda function (or a functor) that evaluates the
   * integrand
   * @param cell_test
   * @param cell_trial
   * @param mask_test
   * @param mask_trial
   * @param output
   */
  template <class TestComponent, class TrialComponent, class Integrand,
            class Kernel>
  void integrate_singular(
      Integrand integrand, Kernel kernel, TestComponent testComponent,
      TrialComponent trialComponent, const DoFParent<2, 3, Real> &cell_test,
      const DoFParent<2, 3, Real> &cell_trial, const DoFMask &mask_test,
      const DoFMask &mask_trial,
      const MultiIndex<2, unsigned int> &local_adjacent_edge_indices,
      DenseSubMatrix<CoefficientType> *output);

private:
  using _theta_bounds_function = void (*)(Real &, Real &);
  _theta_bounds_function get_current_theta_bounds;

  using _psi_bounds_function = void (*)(const Real &, const Real &, Real &,
                                        Real &);
  _psi_bounds_function get_current_psi_bounds;

  using _u_bounds_function = void (*)(const Real &, const Real &, const Real &,
                                      Real &, Real &);
  _u_bounds_function get_current_u_bounds;

  using _lambda_bounds_function = void (*)(const Real &, const Real &,
                                           const Real &, const Real &,
                                           const Real &, Real &, Real &);
  _lambda_bounds_function get_current_lambda_bounds;

  using _uv_coords_mapping_function =
      void (*)(Point<2, Real> &uv);
  _uv_coords_mapping_function map_u_and_v_coords_test;
  _uv_coords_mapping_function map_u_and_v_coords_trial;

  const bool has_double[8] = {false, true, true, true, true, true, false, true};

  unsigned int int_0_quadrature_index;
  unsigned int int_1_quadrature_index;
  unsigned int int_2_quadrature_index;
  unsigned int int_3_quadrature_index;

  Quadrature::QuadratureCollection<Quadrature::GaussQuadrature> qcollection;

  const unsigned int n_sub_integrals = 8;
  void set_outer_integral_bounds_function_pointers(
      const unsigned int &sub_integral_index);
  void set_inner_integral_bounds(const unsigned int &sub_integral_index,
                                 const unsigned int &part_index);
  void set_u_and_v_coords_mapping(
      const MultiIndex<2, unsigned int> &local_adjacent_edge_indices);
};
template <class CoefficientType, class Real>
DIRECTFN_ET2<DoFParent<2, 3, Real>, CoefficientType, Real>::DIRECTFN_ET2(
    const unsigned int &N0, const unsigned int &N1, const unsigned int &N2,
    const unsigned int &N3)
{
  int_0_quadrature_index = N0;
  int_1_quadrature_index = N1;
  int_2_quadrature_index = N2;
  int_3_quadrature_index = N3;

  qcollection.push_back(N0);
  qcollection.push_back(N1);
  qcollection.push_back(N2);
  qcollection.push_back(N3);
}
template <class CoefficientType, class Real>
void DIRECTFN_ET2<DoFParent<2, 3, Real>, CoefficientType, Real>::
    set_u_and_v_coords_mapping(
        const MultiIndex<2, unsigned int> &local_adjacent_edge_indices) {
  switch (local_adjacent_edge_indices(0)) {
  case 1:
    map_u_and_v_coords_test = ET_map_u_and_v_coords_1<Real>;
    break;
  case 2:
    map_u_and_v_coords_test = ET_map_u_and_v_coords_2<Real>;
    break;
  case 3:
    map_u_and_v_coords_test = ET_map_u_and_v_coords_3<Real>;
    break;
  default:
    map_u_and_v_coords_test = ET_map_u_and_v_coords_0<Real>;
  }

  switch (local_adjacent_edge_indices(1)) {
  case 1:
    map_u_and_v_coords_trial = ET_map_u_and_v_coords_1<Real>;
    break;
  case 2:
    map_u_and_v_coords_trial = ET_map_u_and_v_coords_2<Real>;
    break;
  case 3:
    map_u_and_v_coords_trial = ET_map_u_and_v_coords_3<Real>;
    break;
  default:
    map_u_and_v_coords_trial = ET_map_u_and_v_coords_0<Real>;
  }
}

template <class CoefficientType, class Real>
void DIRECTFN_ET2<DoFParent<2, 3, Real>, CoefficientType, Real>::
    set_inner_integral_bounds(const unsigned int &sub_integral_index,
                              const unsigned int &part_index) {
  bool is_lower_region = true;
  if (sub_integral_index == 7) {
    if (part_index == 0)
      is_lower_region = false;
  } else if (part_index > 0)
    is_lower_region = false;
  switch (sub_integral_index) {
  case 0:
    get_current_u_bounds = ET_u_bounds_0_or_6<Real>;
    get_current_lambda_bounds = ET_lambda_bounds_0<Real>;
    break;
  case 1:
    if (is_lower_region) {
      get_current_u_bounds = ET_u_bounds_1_lower<Real>;
    } else {
      get_current_u_bounds = ET_u_bounds_1_upper<Real>;
    }
    if (part_index == 0) {
      get_current_lambda_bounds = ET_lambda_bounds_1_2_or_3_region_0<Real>;
    } else {
      get_current_lambda_bounds = ET_lambda_bounds_1_3_or_7_region_1<Real>;
    }
    break;
  case 2:
    if (is_lower_region) {
      get_current_u_bounds = ET_u_bounds_2_lower<Real>;
    } else {
      get_current_u_bounds = ET_u_bounds_2_upper<Real>;
    }
    if (part_index == 0) {
      get_current_lambda_bounds = ET_lambda_bounds_1_2_or_3_region_0<Real>;
    } else {
      get_current_lambda_bounds = ET_lambda_bounds_2_region_1<Real>;
    }

    break;
  case 3:
    if (is_lower_region) {
      get_current_u_bounds = ET_u_bounds_3_lower<Real>;
    } else {
      get_current_u_bounds = ET_u_bounds_3_upper<Real>;
    }
    if (part_index == 0) {
      get_current_lambda_bounds = ET_lambda_bounds_1_2_or_3_region_0<Real>;
    } else {
      get_current_lambda_bounds = ET_lambda_bounds_1_3_or_7_region_1<Real>;
    }
    break;
  case 4:
    if (is_lower_region) {
      get_current_u_bounds = ET_u_bounds_4_lower<Real>;
    } else {
      get_current_u_bounds = ET_u_bounds_4_upper<Real>;
    }
    if (part_index == 0) {
      get_current_lambda_bounds = ET_lambda_bounds_4_region_0<Real>;
    } else {
      get_current_lambda_bounds = ET_lambda_bounds_4_or_5_region_1<Real>;
    }
    break;
  case 5:
    if (is_lower_region) {
      get_current_u_bounds = ET_u_bounds_5_or_7_lower<Real>;
    } else {
      get_current_u_bounds = ET_u_bounds_5_or_7_upper<Real>;
    }
    if (part_index == 0) {
      get_current_lambda_bounds = ET_lambda_bounds_5_region_0<Real>;
    } else {
      get_current_lambda_bounds = ET_lambda_bounds_4_or_5_region_1<Real>;
    }
    break;
  case 6:
    get_current_lambda_bounds = ET_lambda_bounds_6<Real>;
    get_current_u_bounds = ET_u_bounds_0_or_6<Real>;
    break;
  case 7:
    if (is_lower_region) {
      get_current_u_bounds = ET_u_bounds_5_or_7_lower<Real>;
    } else {
      get_current_u_bounds = ET_u_bounds_5_or_7_upper<Real>;
    }
    if (part_index == 0) {
      get_current_lambda_bounds = ET_lambda_bounds_7_region_0<Real>;
    } else {
      get_current_lambda_bounds = ET_lambda_bounds_1_3_or_7_region_1<Real>;
    }
    break;
  default:
    assert(false && "Invalid sub_integral_index");
  }
}
template <class CoefficientType, class Real>
void DIRECTFN_ET2<DoFParent<2, 3, Real>, CoefficientType, Real>::
    set_outer_integral_bounds_function_pointers(
        const unsigned int &sub_integral_index) {
  switch (sub_integral_index) {
  case 0:
    get_current_theta_bounds = ET_theta_bounds_0_or_1<Real>;
    get_current_psi_bounds = ET_psi_bounds_0<Real>;
    break;
  case 1:
    get_current_theta_bounds = ET_theta_bounds_0_or_1<Real>;
    get_current_psi_bounds = ET_psi_bounds_1<Real>;
    break;
  case 2:
    get_current_theta_bounds = ET_theta_bounds_2_or_3<Real>;
    get_current_psi_bounds = ET_psi_bounds_2<Real>;
    break;
  case 3:
    get_current_theta_bounds = ET_theta_bounds_2_or_3<Real>;
    get_current_psi_bounds = ET_psi_bounds_3<Real>;
    break;
  case 4:
    get_current_theta_bounds = ET_theta_bounds_4_or_5<Real>;
    get_current_psi_bounds = ET_psi_bounds_4<Real>;
    break;
  case 5:
    get_current_theta_bounds = ET_theta_bounds_4_or_5<Real>;
    get_current_psi_bounds = ET_psi_bounds_5<Real>;
    break;
  case 6:
    get_current_theta_bounds = ET_theta_bounds_6_or_7<Real>;
    get_current_psi_bounds = ET_psi_bounds_6<Real>;
    break;
  case 7:
    get_current_theta_bounds = ET_theta_bounds_6_or_7<Real>;
    get_current_psi_bounds = ET_psi_bounds_7<Real>;
    break;
  default:
    assert(false && "Invalid sub_integral_index!");
  }
}
template <class CoefficientType, class Real>
template <class TestComponent, class TrialComponent, class Integrand,
          class Kernel>
void DIRECTFN_ET2<DoFParent<2, 3, Real>, CoefficientType, Real>::
    integrate_singular(
        Integrand integrand, Kernel kernel, TestComponent testComponent,
        TrialComponent trialComponent, const DoFParent<2, 3, Real> &cell_test,
        const DoFParent<2, 3, Real> &cell_trial, const DoFMask &mask_test,
        const DoFMask &mask_trial,
        const MultiIndex<2, unsigned int> &local_adjacent_edge_indices,
        DenseSubMatrix<CoefficientType> *output) {
  using ComponentsReturnType = typename std::result_of<TestComponent(
      const DoFParent<2, 3, Real> &, const unsigned int &,
      const Point<2, Real> &)>::type;

  this->set_u_and_v_coords_mapping(local_adjacent_edge_indices);
  // First, we create a temporary Matrix, as this is more efficient than
  // accessing the DenseSubMatrix frequently While the mask might have disabled
  // entries, given that the number of DoFs on cells will be fairly low, the
  // difference is not significant
  ContiguousMatrix<CoefficientType> temp_matrix(mask_test.size(),
                                                mask_trial.size());

  // This function is for the instances where one element was previously filled,
  // and has been refined. Therefore, a subset of the interactions have already been
  // filled.
  ContiguousMatrix<char> is_filled(mask_test.size(), mask_trial.size());
  output->generate_filled_data(cell_test, cell_trial, mask_test, mask_trial, &is_filled);

  std::vector<char> test_has_integrals_to_perform(mask_test.size(), false);
  std::vector<char> trial_has_integrals_to_perform(mask_trial.size(), false);
  for (unsigned int i = 0; i < mask_test.size(); ++i)
    for (unsigned int j = 0; j < mask_trial.size(); ++j)
      if (is_filled(i,j) == false)
      {
        test_has_integrals_to_perform[i] = true;
        break;
      }
  for (unsigned int j = 0; j < mask_trial.size(); ++j)
    for (unsigned int i = 0; i < mask_test.size(); ++i)
      if (is_filled(i,j) == false)
      {
        trial_has_integrals_to_perform[j] = true;
        break;
      }
 // KahanContiguousMatrix<CoefficientType> temp_matrix(mask_test.size(), mask_trial.size());

  auto cell_test_geom = cell_test.get_cell_on_mesh();
  auto cell_trial_geom = cell_trial.get_cell_on_mesh();
  std::complex<Real> complexj = {0., 1.};
  std::vector<double> *xgl_theta, *wgl_theta;
  qcollection.get_weights_and_points(this->int_0_quadrature_index, wgl_theta,
                                     xgl_theta);
  std::vector<double> *xgl_psi, *wgl_psi;
  qcollection.get_weights_and_points(this->int_1_quadrature_index, wgl_psi,
                                     xgl_psi);
  std::vector<double> *xgl_u, *wgl_u;
  qcollection.get_weights_and_points(this->int_2_quadrature_index, wgl_u,
                                     xgl_u);
  std::vector<double> *xgl_lambda, *wgl_lambda;
  qcollection.get_weights_and_points(this->int_3_quadrature_index, wgl_lambda,
                                     xgl_lambda);
  // We have the 8 primary sub integrals (some of which split into 2 additional
  // integrals)
  for (unsigned sub_integral_index = 0; sub_integral_index < n_sub_integrals;
       ++sub_integral_index) {
    this->set_outer_integral_bounds_function_pointers(sub_integral_index);
    // We have the outer integral over \theta
    Real theta_lower, theta_upper;
    this->get_current_theta_bounds(theta_lower, theta_upper);
    // Get transformation to work with Gauss quadrature...
    const Real theta_average = (theta_lower + theta_upper) / (Real(2.));
    const Real theta_half_length = (theta_upper - theta_lower) / (Real(2.));
    // Get Gauss quadrature parameters

    for (unsigned int theta_index = 0; theta_index < wgl_theta->size();
         ++theta_index) {
      const Real jac0 = theta_half_length * wgl_theta->operator[](theta_index);
      // Get the value of psi
      const Real theta =
          theta_half_length * xgl_theta->operator[](theta_index) +
          theta_average;
      const Real cos_theta = cos(theta);
      const Real sin_theta = sin(theta);
      const Real tan_theta = tan(theta);
      const Real atan_cos_theta = atan(cos(theta));
      const Real atan_sin_theta = atan(sin(theta));

      // Now we have the integral over \Psi
      Real psi_lower, psi_upper;
      this->get_current_psi_bounds(atan_cos_theta, atan_sin_theta, psi_lower,
                                   psi_upper);
      const Real psi_average = (psi_lower + psi_upper) / (Real(2.));
      const Real psi_half_length = (psi_upper - psi_lower) / (Real(2.));
      // Get Gauss quadrature parameters

      for (unsigned int psi_index = 0; psi_index < wgl_psi->size();
           ++psi_index) {
        const Real jac1 =
            psi_half_length * wgl_psi->operator[](psi_index) * jac0;
        const Real psi =
            psi_half_length * xgl_psi->operator[](psi_index) + psi_average;
        const Real cos_psi = cos(psi);
        const Real sin_psi = sin(psi);
        const Real tan_psi = tan(psi);
        // Now we have the integral over u, which, for certain sub integrals,
        // actually splits into a pair of integrals as determined by is_lower
        // region...
        Real u_lower, u_upper;
        const unsigned int n_parts =
            this->has_double[sub_integral_index] ? 2 : 1;
        for (unsigned part_index = 0; part_index < n_parts; ++part_index) {
          set_inner_integral_bounds(sub_integral_index, part_index);

          get_current_u_bounds(cos_theta, tan_theta, tan_psi, u_lower, u_upper);

          const Real u_average = (u_lower + u_upper) / (Real(2.));
          const Real u_half_length = (u_upper - u_lower) / (Real(2.));
          // Get Gauss quadrature parameters

          for (unsigned int u_index = 0; u_index < wgl_u->size(); ++u_index) {
            const Real jac2 = u_half_length * wgl_u->operator[](u_index) * jac1;
            const Real u =
                u_half_length * xgl_u->operator[](u_index) + u_average;
            Real lambda_lower, lambda_upper;

            get_current_lambda_bounds(cos_theta, sin_theta, cos_psi, sin_psi, u,
                                      lambda_lower, lambda_upper);
            const Real lambda_average =
                (lambda_lower + lambda_upper) / (Real(2.));
            const Real lambda_half_length =
                (lambda_upper - lambda_lower) / (Real(2.));
            // Get Gauss quadrature parameters

            for (unsigned int lambda_index = 0;
                 lambda_index < wgl_lambda->size(); ++lambda_index) {
              const Real lambda =
                  lambda_half_length * xgl_lambda->operator[](lambda_index) +
                  lambda_average;
              // Now we need to get the actual uv coordinates
              auto tester = sin(psi);
              Point<2, Real> uv_test = {u, lambda * sin_psi - Real(1.0)};
              Point<2, Real> uv_trial = {lambda * cos_psi * cos_theta - u,
                                         lambda * cos_psi * sin_theta -
                                             Real(1.)};
              const Real transformation_jacobian = lambda * lambda * cos_psi;
              const Real full_jacobian = transformation_jacobian *
                                         wgl_lambda->operator[](lambda_index) *
                                         lambda_half_length * jac2;
              // However, in constructing these points, we have assumed a
              // particular orientation
              // of the edge adjacency
              // Since this adjacency is not satisfied in most cases, we must
              // introduce the correct transformation
              this->map_u_and_v_coords_test(uv_test);
              this->map_u_and_v_coords_trial(uv_trial);

              const auto unitary_test_u =
                  cell_test_geom->unitary_vector(uv_test, u_dir);
              const auto unitary_test_v =
                  cell_test_geom->unitary_vector(uv_test, v_dir);
              const auto unitary_trial_u =
                  cell_trial_geom->unitary_vector(uv_trial, u_dir);
              const auto unitary_trial_v =
                  cell_trial_geom->unitary_vector(uv_trial, v_dir);
              const auto r_test = cell_test_geom->r(uv_test);
              const auto r_trial = cell_trial_geom->r(uv_trial);
              const Real R = (r_test - r_trial).norm();
              const CoefficientType kernel_value = kernel(R);

              // Now that we've computed uv_trial, we can precompute the values
              // of the shape functions
              ComponentsReturnType dof_trial_memoize[mask_trial.size()];
              for (unsigned int dof_trial_index = 0;
                   dof_trial_index < mask_trial.size(); ++dof_trial_index) {
                if (!trial_has_integrals_to_perform[dof_trial_index])
                  continue;

                dof_trial_memoize[dof_trial_index] =
                    trialComponent(cell_trial, dof_trial_index, uv_trial);
              }

              for (unsigned int dof_test_index = 0;
                   dof_test_index < mask_test.size(); ++dof_test_index) {
                if (!test_has_integrals_to_perform[dof_test_index])
                  continue;
                const DoFBase<cell_test.cell_dim> &dof_test =
                    cell_test.get_dof(dof_test_index);
                const auto test_component_value =
                    testComponent(cell_test, dof_test_index, uv_test);

                for (unsigned int dof_trial_index = 0;
                     dof_trial_index < mask_trial.size(); ++dof_trial_index) {

                  if (is_filled(dof_test_index, dof_trial_index))
                    continue;
                  const DoFBase<cell_trial.cell_dim> &dof_trial =
                      cell_trial.get_dof(dof_trial_index);

                  const auto &trial_component_value =
                      dof_trial_memoize[dof_trial_index];

                  const auto value_to_add = full_jacobian *
                                            integrand(test_component_value, trial_component_value,
                                                      dof_test.direction == u_dir ? unitary_test_u
                                                                                  : unitary_test_v,
                                                      dof_trial.direction == u_dir ? unitary_trial_u
                                                                                   : unitary_trial_v,
                                                      kernel_value, R);

                  // temp_matrix(dof_test_index, dof_trial_index) += value_to_add;
                  temp_matrix.accumulate(dof_test_index, dof_trial_index,value_to_add);

                }
              }
            }
          }
        }
      }
    }
  }
  temp_matrix.apply_carry();
  for (unsigned int i = 0; i < mask_test.size(); ++i) {
    if (!test_has_integrals_to_perform[i])
      continue;
    const auto &dof_test = cell_test.get_dof(i);
    for (unsigned int j = 0; j < mask_trial.size(); ++j) {
      if (is_filled(i,j))
        continue;
      const auto &dof_trial = cell_trial.get_dof(j);

      output->at(dof_test.global_index, dof_trial.global_index) +=
          temp_matrix(i, j);
      output->get_is_filled_at(dof_test.global_index, dof_trial.global_index) = true;
    }
  }
}


} // namespace NumericalIntegration
DROMON_NAMESPACE_CLOSE
#endif // DROMON_DIRECTFN_ET_SINGULAR_H
