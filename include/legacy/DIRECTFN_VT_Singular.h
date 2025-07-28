//
// Created by Jake J. Harmon (jake.harmon@ieee.org) on 5/1/22.
//

#ifndef DROMON_DIRECTFN_VT_SINGULAR_H
#define DROMON_DIRECTFN_VT_SINGULAR_H
#include "DIRECTFN_VT_Bounds_Functions.h"
#include "DoFBase.h"
#include "DoFGeom.h"
#include "DoFMask.h"
#include "Kernels.h"
#include "QuadratureCollection.h"
#include "SubMatrix.h"
#include "config.h"
DROMON_NAMESPACE_OPEN
namespace NumericalIntegration {
template <class DoFCellType, class CoefficientType = std::complex<double>,
          class Real = double>
class DIRECTFN_VT2 {
public:
  DIRECTFN_VT2() {
    assert(false && "The general implementation of the DIRECTFN_VT numerical "
                    "integrator is not yet implemented!");
  }

private:
};

template <class CoefficientType, class Real>
class DIRECTFN_VT2<DoFParent<2, 3, Real>, CoefficientType, Real> {
public:
  /**
   *
   * @param N0 - Order of numerical integration for outermost integral
   * @param N1 - See above, for next integral
   * @param N2 - See above, for next integral
   * @param N3 - See above, for innermost integral
   */
  DIRECTFN_VT2(const unsigned int &N0, const unsigned int &N1,
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
      const MultiIndex<2, unsigned int> &local_adjacent_vertex_indices,
      DenseSubMatrix<CoefficientType> *output);

private:
  using _theta_bounds_function = void (*)(Real &, Real &);
  _theta_bounds_function get_current_theta_test_bounds;
  _theta_bounds_function get_current_theta_trial_bounds;

  using _L_bounds_function = void (*)(const Real &, const Real &,
                                      Real &);
  _L_bounds_function get_current_L_test;
  _L_bounds_function get_current_L_trial;

  using _uv_coords_mapping_function = void (*)(Point<2, Real> &uv);
  _uv_coords_mapping_function map_u_and_v_coords_test;
  _uv_coords_mapping_function map_u_and_v_coords_trial;

  static constexpr unsigned int n_sub_ranges = 4;

  unsigned int int_0_quadrature_index;
  unsigned int int_1_quadrature_index;
  unsigned int int_2_quadrature_index;
  unsigned int int_3_quadrature_index;

  Quadrature::QuadratureCollection<Quadrature::GaussQuadrature> qcollection;

  void set_primary_bounds(const unsigned int &sub_range_index);
  void set_u_and_v_coords_mapping(
      const MultiIndex<2, unsigned int> &local_adjacent_vertex_indices);
};
template <class CoefficientType, class Real>
DIRECTFN_VT2<DoFParent<2, 3, Real>, CoefficientType, Real>::DIRECTFN_VT2(
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
void DIRECTFN_VT2<DoFParent<2, 3, Real>, CoefficientType, Real>::
    set_u_and_v_coords_mapping(
        const MultiIndex<2, unsigned int> &local_adjacent_vertex_indices) {
  switch (local_adjacent_vertex_indices(0)) {
  case 1:
    map_u_and_v_coords_test = VT_map_u_and_v_coords_1<Real>;
    break;
  case 2:
    map_u_and_v_coords_test = VT_map_u_and_v_coords_2<Real>;
    break;
  case 3:
    map_u_and_v_coords_test = VT_map_u_and_v_coords_3<Real>;
    break;
  default:
    map_u_and_v_coords_test = VT_map_u_and_v_coords_0<Real>;
  }
  switch (local_adjacent_vertex_indices(1)) {
  case 1:
    map_u_and_v_coords_trial = VT_map_u_and_v_coords_1<Real>;
    break;
  case 2:
    map_u_and_v_coords_trial = VT_map_u_and_v_coords_2<Real>;
    break;
  case 3:
    map_u_and_v_coords_trial = VT_map_u_and_v_coords_3<Real>;
    break;
  default:
    map_u_and_v_coords_trial = VT_map_u_and_v_coords_0<Real>;
  }
}

template <class CoefficientType, class Real>
template <class TestComponent, class TrialComponent, class Integrand,
          class Kernel>
void DIRECTFN_VT2<DoFParent<2, 3, Real>, CoefficientType, Real>::
    integrate_singular(
        Integrand integrand, Kernel kernel, TestComponent testComponent,
        TrialComponent trialComponent, const DoFParent<2, 3, Real> &cell_test,
        const DoFParent<2, 3, Real> &cell_trial, const DoFMask &mask_test,
        const DoFMask &mask_trial,
        const MultiIndex<2, unsigned int> &local_adjacent_vertex_indices,
        DenseSubMatrix<CoefficientType> *output) {
  using ComponentsReturnType = typename std::result_of<TestComponent(
      const DoFParent<2, 3, Real> &, const unsigned int &,
      const Point<2, Real> &)>::type;

  this->set_u_and_v_coords_mapping(local_adjacent_vertex_indices);

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
  std::vector<double> *xgl_theta_test, *wgl_theta_test;
  qcollection.get_weights_and_points(this->int_0_quadrature_index,
                                     wgl_theta_test, xgl_theta_test);
  std::vector<double> *xgl_theta_trial, *wgl_theta_trial;
  qcollection.get_weights_and_points(this->int_1_quadrature_index,
                                     wgl_theta_trial, xgl_theta_trial);
  std::vector<double> *xgl_psi, *wgl_psi;
  qcollection.get_weights_and_points(this->int_2_quadrature_index, wgl_psi,
                                     xgl_psi);
  std::vector<double> *xgl_lambda, *wgl_lambda;
  qcollection.get_weights_and_points(this->int_3_quadrature_index, wgl_lambda,
                                     xgl_lambda);
  //std::ofstream uv_out("uv_points_out.txt");
  //uv_out << std::setprecision(16);
  for (unsigned int sub_range_index = 0; sub_range_index < this->n_sub_ranges;
       ++sub_range_index) {

    this->set_primary_bounds(sub_range_index);
    Real theta_test_lower, theta_test_upper;
    this->get_current_theta_test_bounds(theta_test_lower, theta_test_upper);

    const Real theta_test_average =
        (theta_test_lower + theta_test_upper) / (Real(2.));
    const Real theta_test_half_length =
        (theta_test_upper - theta_test_lower) / (Real(2.));

    for (unsigned int theta_test_index = 0;
         theta_test_index < wgl_theta_test->size(); ++theta_test_index) {
      const Real jac0 =
          theta_test_half_length * wgl_theta_test->operator[](theta_test_index);
      const Real theta_test = theta_test_half_length *
                                  xgl_theta_test->operator[](theta_test_index) +
                              theta_test_average;
      const Real cos_theta_test = cos(theta_test);
      const Real sin_theta_test = sin(theta_test);
      Real L_test;
      this->get_current_L_test(cos_theta_test, sin_theta_test, L_test);

      Real theta_trial_lower, theta_trial_upper;
      this->get_current_theta_trial_bounds(theta_trial_lower,
                                           theta_trial_upper);

      const Real theta_trial_average =
          (theta_trial_lower + theta_trial_upper) / (Real(2.));
      const Real theta_trial_half_length =
          (theta_trial_upper - theta_trial_lower) / (Real(2.));
      // Get Gauss quadrature parameters

      for (unsigned int theta_trial_index = 0;
           theta_trial_index < wgl_theta_trial->size(); ++theta_trial_index) {
        const Real jac1 = jac0 * theta_trial_half_length *
                          wgl_theta_trial->operator[](theta_trial_index);
        const Real theta_trial =
            theta_trial_half_length *
                xgl_theta_trial->operator[](theta_trial_index) +
            theta_trial_average;
        const Real cos_theta_trial = cos(theta_trial);
        const Real sin_theta_trial = sin(theta_trial);

        Real L_trial;
        this->get_current_L_trial(cos_theta_trial, sin_theta_trial, L_trial);
        const Real atan_Ltrial_over_Ltest = atan(L_trial / L_test);

        Real psi_lower, psi_upper, psi_average, psi_half_length;
        // For the first integral part
        // this->get_psi_bounds_I(atan_Ltrial_over_Ltest, &psi_lower,
        // &psi_upper);
        VT_psi_bounds_I(atan_Ltrial_over_Ltest, psi_lower, psi_upper);
        psi_average = (psi_upper + psi_lower) / Real(2.0);
        psi_half_length = (psi_upper - psi_lower) / Real(2.0);
        for (unsigned int psi_index = 0; psi_index < wgl_psi->size();
             ++psi_index) {
          const Real jac2 =
              jac1 * psi_half_length * wgl_psi->operator[](psi_index);
          const Real psi =
              psi_half_length * xgl_psi->operator[](psi_index) + psi_average;
          const Real cos_psi = cos(psi);
          const Real sin_psi = sin(psi);

          const Real lambda_lower = Real(0.0);
          const Real lambda_upper = L_test / cos_psi;
          const Real lambda_average = lambda_upper / Real(2.0);
          const Real lambda_half_length = lambda_upper / Real(2.0);

          for (unsigned int lambda_index = 0; lambda_index < wgl_lambda->size();
               ++lambda_index) {

            const Real lambda =
                lambda_half_length *
                (xgl_lambda->operator[](lambda_index) + Real(1.0));

            Point<2, Real> uv_test = {
                lambda * cos_psi * cos_theta_test - Real(1.0),
                lambda * cos_psi * sin_theta_test - Real(1.0)};
            Point<2, Real> uv_trial = {
                lambda * sin_psi * cos_theta_trial - Real(1.0),
                lambda * sin_psi * sin_theta_trial - Real(1.0)};
//            uv_out << uv_test(0) << " " << uv_test(1) << " " <<
//                uv_trial(0) << " " << uv_trial(1) << std::endl;
            const Real transformation_jacobian =
                lambda * lambda * lambda * cos_psi * sin_psi;
            const Real full_jacobian = jac2 * transformation_jacobian *
                                       lambda_half_length *
                                       wgl_lambda->operator[](lambda_index);
            // Map these uv_coordinates to the necessary ones for the given cell
            // orientation...
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
                // Now, we evaluate the integrand and then multiply by the
                // integration weights...
                const auto value_to_add = full_jacobian *
                                          integrand(test_component_value, trial_component_value,
                                                    dof_test.direction == u_dir ? unitary_test_u
                                                                                : unitary_test_v,
                                                    dof_trial.direction == u_dir ? unitary_trial_u
                                                                                 : unitary_trial_v,
                                                    kernel_value);
                // temp_matrix(dof_test_index, dof_trial_index) += value_to_add;
                temp_matrix.accumulate(dof_test_index, dof_trial_index, value_to_add);

              }
            }
          }
        }

        // For the second integral part
        // this->get_psi_bounds_II(atan_Ltrial_over_Ltest, &psi_lower,
        // &psi_upper);
        VT_psi_bounds_II(atan_Ltrial_over_Ltest, psi_lower, psi_upper);
        psi_average = (psi_upper + psi_lower) / Real(2.0);
        psi_half_length = (psi_upper - psi_lower) / Real(2.0);
        for (unsigned int psi_index = 0; psi_index < wgl_psi->size();
             ++psi_index) {
          const Real jac2 =
              jac1 * psi_half_length * wgl_psi->operator[](psi_index);
          const Real psi =
              psi_half_length * xgl_psi->operator[](psi_index) + psi_average;
          const Real cos_psi = cos(psi);
          const Real sin_psi = sin(psi);

          const Real lambda_lower = Real(0.0);
          const Real lambda_upper = L_trial / sin_psi;
          const Real lambda_average = lambda_upper / Real(2.0);
          const Real lambda_half_length = lambda_upper / Real(2.0);

          for (unsigned int lambda_index = 0; lambda_index < wgl_lambda->size();
               ++lambda_index) {

            const Real lambda =
                lambda_half_length *
                (xgl_lambda->operator[](lambda_index) + Real(1.0));

            Point<2, Real> uv_test = {
                lambda * cos_psi * cos_theta_test - Real(1.0),
                lambda * cos_psi * sin_theta_test - Real(1.0)};
            Point<2, Real> uv_trial = {
                lambda * sin_psi * cos_theta_trial - Real(1.0),
                lambda * sin_psi * sin_theta_trial - Real(1.0)};
//            uv_out << uv_test(0) << " " << uv_test(1) << " " <<
//                uv_trial(0) << " " << uv_trial(1) << std::endl;
            const Real transformation_jacobian =
                lambda * lambda * lambda * cos(psi) * sin(psi);
            const Real full_jacobian = jac2 * transformation_jacobian *
                                       lambda_half_length *
                                       wgl_lambda->operator[](lambda_index);
            // Map these uv_coordinates to the necessary ones for the given cell
            // orientation...
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
                // Now, we evaluate the integrand and then multiply by the
                // integration weights...
                const auto value_to_add = full_jacobian *
                                          integrand(test_component_value, trial_component_value,
                                                    dof_test.direction == u_dir ? unitary_test_u
                                                                                : unitary_test_v,
                                                    dof_trial.direction == u_dir ? unitary_trial_u
                                                                                 : unitary_trial_v,
                                                    kernel_value);
                // temp_matrix(dof_test_index, dof_trial_index) += value_to_add;
                temp_matrix.accumulate(dof_test_index, dof_trial_index, value_to_add);
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

template <class CoefficientType, class Real>
void DIRECTFN_VT2<DoFParent<2, 3, Real>, CoefficientType, Real>::
    set_primary_bounds(const unsigned int &sub_range_index) {
  switch (sub_range_index) {
  case 0:
    get_current_theta_test_bounds = VT_theta_first<Real>;
    get_current_theta_trial_bounds = VT_theta_first<Real>;
    get_current_L_test = VT_Ltest_0_or_1<Real>;
    get_current_L_trial = VT_Ltrial_0_or_2<Real>;
    break;
  case 1:
    get_current_theta_test_bounds = VT_theta_first<Real>;
    get_current_theta_trial_bounds = VT_theta_second<Real>;
    get_current_L_test = VT_Ltest_0_or_1<Real>;
    get_current_L_trial = VT_Ltrial_1_or_3<Real>;
    break;
  case 2:
    get_current_theta_test_bounds = VT_theta_second<Real>;
    get_current_theta_trial_bounds = VT_theta_first<Real>;
    get_current_L_test = VT_Ltest_2_or_3<Real>;
    get_current_L_trial = VT_Ltrial_0_or_2<Real>;
    break;
  case 3:
    get_current_theta_test_bounds = VT_theta_second<Real>;
    get_current_theta_trial_bounds = VT_theta_second<Real>;
    get_current_L_test = VT_Ltest_2_or_3<Real>;
    get_current_L_trial = VT_Ltrial_1_or_3<Real>;
    break;
  default:
    assert(false && "Invalid sub_range_index");
  }
}
} // namespace NumericalIntegration
DROMON_NAMESPACE_CLOSE
#endif // DROMON_DIRECTFN_VT_SINGULAR_H
