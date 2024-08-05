//
// Created by Jake J. Harmon (jake.harmon@ieee.org) on 4/29/22.
//

#ifndef DROMON_DIRECTFN_SINGULAR_MOD_H
#define DROMON_DIRECTFN_SINGULAR_MOD_H
#include "DIRECTFN_ST_Bounds_Functions.h"
#include "DoFBase.h"
#include "DoFGeom.h"
#include "DoFMask.h"
#include "Kernels.h"
#include "QuadratureCollection.h"
#include "SubMatrix.h"
#include "config.h"
#include <functional>
#include <iomanip>
#include <boost/math/quadrature/gauss_kronrod.hpp>

DROMON_NAMESPACE_OPEN
namespace NumericalIntegration {
template <class DoFCellType, class CoefficientType = std::complex<double>,
          class Real = double>
class DIRECTFN_ST2 {
public:
  DIRECTFN_ST2() {
    assert(false && "The general implementation of the DIRECTFN_ST numerical "
                    "integrator is not yet implemented!");
  }

private:
};

template <class CoefficientType, class Real>
class DIRECTFN_ST2<DoFParent<2, 3, Real>, CoefficientType, Real> {
public:
  /**
   *
   * @param N0 - Order of numerical integration for outermost integral
   * @param N1 - See above, for next integral
   * @param N2 - See above, for next integral
   * @param N3 - See above, for innermost integral
   */
  DIRECTFN_ST2(const unsigned int &N0, const unsigned int &N1,
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
  void integrate_singular(Integrand integrand, Kernel kernel,
                          TestComponent testComponent,
                          TrialComponent trialComponent,
                          const DoFParent<2, 3, Real> &cell_test,
                          const DoFParent<2, 3, Real> &cell_trial,
                          const DoFMask &mask_test, const DoFMask &mask_trial,
                          DenseSubMatrix<CoefficientType> *output);
  template <class TestComponent, class TrialComponent, class Integrand,
            class Kernel>
  void integrate_singular_adaptive(Integrand integrand, Kernel kernel,
                                   TestComponent testComponent,
                                   TrialComponent trialComponent,
                                   const DoFParent<2, 3, Real> &cell_test,
                                   const DoFParent<2, 3, Real> &cell_trial,
                                   const DoFMask &mask_test, const DoFMask &mask_trial,
                                   DenseSubMatrix<CoefficientType> *output);

private:
  // Function Pointer for the rho bounds (the inner-most integral)
  // The first entry is the value of lambda
  using _rho_bounds_function = void (*)(const Real &, Real &, Real &);
  _rho_bounds_function get_current_rho_bounds;
  // Function Pointer for the lambda bounds (the penultimate integral)
  // The first value is cos(psi), the second is sin(psi), the third is u
  using _lambda_bounds_function = void (*)(const Real &, const Real &,
                                           const Real &, Real &, Real &,
                                           const Real &psi);
  _lambda_bounds_function get_current_lambda_bounds;

  // Same for the u integral
  // The first entry is tan(psi - constants<Real>::PI/2)
  using _u_bounds_function = void (*)(const Real &, Real &, Real &);
  _u_bounds_function get_current_u_bounds;

  using _psi_bounds_function = void (*)(Real &, Real &);
  _psi_bounds_function get_current_psi_bounds;

  using _uv_coords_mapping_function =
      Point<2, Real> (*)(const Point<2, Real> &uv);
  _uv_coords_mapping_function map_u_and_v_coords;

  void
  set_integral_bounds_function_pointers(const unsigned int &sub_triangle_index,
                                        const unsigned int &sub_integral_index);

  //  void get_psi_bounds(const unsigned int &sub_integral_index, Real
  //  *psi_lower,
  //                      Real *psi_upper) const;
  //
  //  void get_u_bounds(const unsigned int &sub_integral_index, const Real &psi,
  //                    Real *u_lower, Real *u_upper) const;
  //
  //  void get_lambda_bounds(const unsigned int &sub_integral_index,
  //                         const Real &psi, const Real &u, Real *lambda_lower,
  //                         Real *lambda_upper) const;

  //  Point<2, Real> map_u_and_v_coords(const unsigned int &sub_triangle_index,
  //                                    const Point<2, Real> &uv) const;

  const unsigned int n_sub_triangles = 4;
  const unsigned int n_sub_integrals = 6;

  unsigned int int_0_quadrature_index;
  unsigned int int_1_quadrature_index;
  unsigned int int_2_quadrature_index;
  unsigned int int_3_quadrature_index;

  Quadrature::QuadratureCollection<Quadrature::GaussQuadrature> qcollection;
};
template <class CoefficientType, class Real>
void DIRECTFN_ST2<DoFParent<2, 3, Real>, CoefficientType, Real>::
    set_integral_bounds_function_pointers(
        const unsigned int &sub_triangle_index,
        const unsigned int &sub_integral_index) {
  get_current_rho_bounds = ST_rho_bounds<Real>;
  switch (sub_integral_index) {
  case 0:
    get_current_psi_bounds = ST_psi_bounds_0<Real>;
    get_current_u_bounds = ST_u_bounds_0_or_5<Real>;
    get_current_lambda_bounds = ST_lambda_bounds_0_or_1<Real>;
    break;
  case 1:
    get_current_psi_bounds = ST_psi_bounds_1_or_2<Real>;
    get_current_u_bounds = ST_u_bounds_1<Real>;
    get_current_lambda_bounds = ST_lambda_bounds_0_or_1<Real>;
    break;
  case 2:
    get_current_psi_bounds = ST_psi_bounds_1_or_2<Real>;
    get_current_u_bounds = ST_u_bounds_2<Real>;
    get_current_lambda_bounds = ST_lambda_bounds_2_or_3<Real>;
    break;
  case 3:
    get_current_psi_bounds = ST_psi_bounds_3_or_4<Real>;
    get_current_u_bounds = ST_u_bounds_3<Real>;
    get_current_lambda_bounds = ST_lambda_bounds_2_or_3<Real>;
    break;
  case 4:
    get_current_psi_bounds = ST_psi_bounds_3_or_4<Real>;
    get_current_u_bounds = ST_u_bounds_4<Real>;
    get_current_lambda_bounds = ST_lambda_bounds_4_or_5<Real>;
    break;
  case 5:
    get_current_psi_bounds = ST_psi_bounds_5<Real>;
    get_current_u_bounds = ST_u_bounds_0_or_5<Real>;
    get_current_lambda_bounds = ST_lambda_bounds_4_or_5<Real>;
    break;
  default:
    assert(false && "Invalid sub_integral_index!");
  }

  switch (sub_triangle_index) {
  case 0:
    map_u_and_v_coords = ST_map_u_and_v_coords_0<Real>;
    break;
  case 1:
    map_u_and_v_coords = ST_map_u_and_v_coords_1<Real>;
    break;
  case 2:
    map_u_and_v_coords = ST_map_u_and_v_coords_2<Real>;
    break;
  case 3:
    map_u_and_v_coords = ST_map_u_and_v_coords_3<Real>;
    break;
  }
}

template <class CoefficientType, class Real>
template <class TestComponent, class TrialComponent, class Integrand,
          class Kernel>
void DIRECTFN_ST2<DoFParent<2, 3, Real>, CoefficientType, Real>::
    integrate_singular(Integrand integrand, Kernel kernel,
                       TestComponent testComponent,
                       TrialComponent trialComponent,
                       const DoFParent<2, 3, Real> &cell_test,
                       const DoFParent<2, 3, Real> &cell_trial,
                       const DoFMask &mask_test, const DoFMask &mask_trial,
                       DenseSubMatrix<CoefficientType> *output) {
  //  std::ofstream uv_out("st_points_out_cpp.txt");
  //  uv_out << std::setprecision(16);
  //  unsigned int counter = 0;
  using ComponentsReturnType = typename std::result_of<TestComponent(
      const DoFParent<2, 3, Real> &, const unsigned int &,
      const Point<2, Real> &)>::type;

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
  // KahanContiguousMatrix<CoefficientType> temp_matrix(mask_test.size(),
  // mask_trial.size());
  auto cell_test_geom = cell_test.get_cell_on_mesh();
  auto cell_trial_geom = cell_trial.get_cell_on_mesh();

  std::complex<Real> complexj = {0., 1.};

  // For the self term, we have integrations over 4 subtriangles
  // These subtriangles then decompose into 6 4-D integrals
  for (unsigned int sub_triangle_index = 0;
       sub_triangle_index < n_sub_triangles; ++sub_triangle_index)
    for (unsigned int sub_integral_index = 0;
         sub_integral_index < n_sub_integrals; ++sub_integral_index) {
      // Set integral bounds functions
      this->set_integral_bounds_function_pointers(sub_triangle_index,
                                                  sub_integral_index);
      // With the sub_triangle_index and sub_integral_index set, we now have to
      // evaulate the 4-D integral that we have for this case

      // First, get the upper and lower bounds for the outer integral
      Real psi_lower, psi_upper;
      // this->get_psi_bounds(sub_integral_index, &psi_lower, &psi_upper);
      this->get_current_psi_bounds(psi_lower, psi_upper);

      // Get transformation to work with Gauss quadrature...
      const Real psi_average = (psi_lower + psi_upper) / (Real(2.));
      const Real psi_half_length = (psi_upper - psi_lower) / (Real(2.));
      // Get Gauss quadrature parameters
      std::vector<double> *xgl_psi, *wgl_psi;
      qcollection.get_weights_and_points(this->int_0_quadrature_index, wgl_psi,
                                         xgl_psi);

      for (unsigned int psi_index = 0; psi_index < wgl_psi->size();
           ++psi_index) {
        const Real jac0 = psi_half_length * wgl_psi->operator[](psi_index);
        // Get the value of psi
        const Real psi =
            psi_half_length * xgl_psi->operator[](psi_index) + psi_average;
        const Real tan_psi_less_pi_half = tan(psi - constants<Real>::PI / Real(2.0));
        const Real cos_psi = cos(psi);
        const Real sin_psi = sin(psi);

        // Get the upper and lower bounds for "u", the next integral
        Real u_lower, u_upper;
        this->get_current_u_bounds(tan_psi_less_pi_half, u_lower, u_upper);
        // Get transformation to work with Gauss quadrature...
        const Real u_average = (u_lower + u_upper) / (Real(2.));
        const Real u_half_length = (u_upper - u_lower) / (Real(2.));

        // Get the Gauss quadrature parameters
        std::vector<double> *xgl_u, *wgl_u;
        qcollection.get_weights_and_points(this->int_1_quadrature_index, wgl_u,
                                           xgl_u);

        for (unsigned int u_index = 0; u_index < wgl_u->size(); ++u_index) {
          const Real jac1 = jac0 * u_half_length * wgl_u->operator[](u_index);
          const Real u = u_half_length * xgl_u->operator[](u_index) + u_average;

          // Get the upper and lower bounds for \Lambda
          Real lambda_lower, lambda_upper;
          this->get_current_lambda_bounds(cos_psi, sin_psi, u, lambda_lower,
                                          lambda_upper, psi);
          // Get transformation to work with Gauss quadrature...
          const Real lambda_average =
              (lambda_lower + lambda_upper) / (Real(2.));
          const Real lambda_half_length =
              (lambda_upper - lambda_lower) / (Real(2.));

          // Get the Gauss quadrature parameters
          std::vector<double> *xgl_lambda, *wgl_lambda;
          qcollection.get_weights_and_points(this->int_2_quadrature_index,
                                             wgl_lambda, xgl_lambda);

          for (unsigned int lambda_index = 0; lambda_index < wgl_lambda->size();
               ++lambda_index) {
            const Real jac2 = jac1 * lambda_half_length *
                              wgl_lambda->operator[](lambda_index);
            const Real lambda =
                lambda_half_length * xgl_lambda->operator[](lambda_index) +
                lambda_average;

            Point<2, Real> uv_test_unmapped = {u, lambda * sin_psi - Real(1.)};
            // Convert uv_test
            const auto uv_test = this->map_u_and_v_coords(uv_test_unmapped);

            // Now that we've computed uv_test, we can precompute the values of
            // the shape functions
            ComponentsReturnType dof_test_memoize[mask_test.size()];
            for (unsigned int dof_test_index = 0;
                 dof_test_index < mask_test.size(); ++dof_test_index) {
              if (!test_has_integrals_to_perform[dof_test_index])
                continue;
              dof_test_memoize[dof_test_index] =
                  testComponent(cell_test, dof_test_index, uv_test);
            }
            const auto unitary_test_u =
                cell_test_geom->unitary_vector(uv_test, u_dir);
            const auto unitary_test_v =
                cell_test_geom->unitary_vector(uv_test, v_dir);
            const auto r_test = cell_test_geom->r(uv_test);

            // Get the upper and lower bounds for \rho
            const Real rho_lower = Real(0.0);
            const Real rho_upper = lambda;

            const Real rho_average = (rho_upper) / (Real(2.));
            const Real rho_half_length = (rho_upper) / (Real(2.));

            // Get the Gauss quadrature parameters
            std::vector<double> *xgl_rho, *wgl_rho;
            qcollection.get_weights_and_points(this->int_3_quadrature_index,
                                               wgl_rho, xgl_rho);

            for (unsigned int rho_index = 0; rho_index < wgl_rho->size();
                 ++rho_index) {

              const Real rho =
                  rho_half_length * xgl_rho->operator[](rho_index) +
                  rho_average;

              const Real transformation_jacobian = rho * sin_psi;
              const Real full_jacobian = transformation_jacobian * jac2 *
                                         rho_half_length *
                                         wgl_rho->operator[](rho_index);
              //++counter;

              Point<2, Real> uv_trial_unmapped = {
                  uv_test_unmapped(0) + rho * cos_psi,
                  -rho * sin_psi + uv_test_unmapped(1)};
              // Convert uv_trial
              const auto uv_trial = this->map_u_and_v_coords(uv_trial_unmapped);

              const auto unitary_trial_u =
                  cell_trial_geom->unitary_vector(uv_trial, u_dir);
              const auto unitary_trial_v =
                  cell_trial_geom->unitary_vector(uv_trial, v_dir);

              const auto r_trial = cell_trial_geom->r(uv_trial);
              const Real R = (r_test - r_trial).norm();
              const CoefficientType kernel_value = kernel(R);
              //  uv_out << uv_test(0) << " "<< uv_test(1)<< " " << uv_trial(0)
              //  << " " << uv_trial(1) << std::endl;
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
                // const auto test_component_value = testComponent(cell_test,
                // dof_test_index, uv_test);
                const auto &test_component_value =
                    dof_test_memoize[dof_test_index];

                for (unsigned int dof_trial_index = 0;
                     dof_trial_index < mask_trial.size(); ++dof_trial_index) {
                  if (is_filled(dof_test_index, dof_trial_index))
                    continue;
                  const DoFBase<cell_trial.cell_dim> &dof_trial =
                      cell_trial.get_dof(dof_trial_index);
                  // Now, we evaluate the integrand and then multiply by the
                  // integration weights...

                  const auto &trial_component_value =
                      dof_trial_memoize[dof_trial_index];
                  const auto value_to_add =
                      full_jacobian *
                      integrand(test_component_value, trial_component_value,
                                dof_test.direction == u_dir ? unitary_test_u
                                                            : unitary_test_v,
                                dof_trial.direction == u_dir ? unitary_trial_u
                                                             : unitary_trial_v,
                                kernel_value, R);
                  //   temp_matrix(dof_test_index, dof_trial_index) +=
                  //   value_to_add;
                  temp_matrix.accumulate(dof_test_index, dof_trial_index,
                                         value_to_add);
                }
              }
            }
          }
        }
      }
    }
  temp_matrix.apply_carry();
  // Now sync these integrals with the main matrix
  for (unsigned int i = 0; i < mask_test.size(); ++i) {
    if (!test_has_integrals_to_perform[i])
      continue;
    const auto &dof_test = cell_test.get_dof(i);
    for (unsigned int j = 0; j < mask_trial.size(); ++j) {
      if (is_filled(i,j))
        continue;
      const auto &dof_trial = cell_trial.get_dof(j);
      // Below is temp and does not work if the integrals are not symmetric!
      auto value_to_add = Real(0.5)*(temp_matrix(i, j) + temp_matrix(j,i));
      output->at(dof_test.global_index, dof_trial.global_index) += value_to_add;
         // temp_matrix(i, j);

      output->get_is_filled_at(dof_test.global_index, dof_trial.global_index) = true;
    }
  }
  //  uv_out.close();
}

template <class CoefficientType, class Real>
DIRECTFN_ST2<DoFParent<2, 3, Real>, CoefficientType, Real>::DIRECTFN_ST2(
    const unsigned int &N0, const unsigned int &N1, const unsigned int &N2,
    const unsigned int &N3) {
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
template <class TestComponent, class TrialComponent, class Integrand,
          class Kernel>
void DIRECTFN_ST2<DoFParent<2, 3, Real>, CoefficientType, Real>::
    integrate_singular_adaptive(
        Integrand integrand, Kernel kernel, TestComponent testComponent,
        TrialComponent trialComponent, const DoFParent<2, 3, Real> &cell_test,
        const DoFParent<2, 3, Real> &cell_trial, const DoFMask &mask_test,
        const DoFMask &mask_trial, DenseSubMatrix<CoefficientType> *output) {
  //  std::ofstream uv_out("st_points_out_cpp.txt");
  //  uv_out << std::setprecision(16);
  //  unsigned int counter = 0;
  using ComponentsReturnType = typename std::result_of<TestComponent(
      const DoFParent<2, 3, Real> &, const unsigned int &,
      const Point<2, Real> &)>::type;

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
  // KahanContiguousMatrix<CoefficientType> temp_matrix(mask_test.size(),
  // mask_trial.size());
  auto cell_test_geom = cell_test.get_cell_on_mesh();
  auto cell_trial_geom = cell_trial.get_cell_on_mesh();

  std::complex<Real> complexj = {0., 1.};

  // For the self term, we have integrations over 4 subtriangles
  // These subtriangles then decompose into 6 4-D integrals
  for (unsigned int sub_triangle_index = 0;
       sub_triangle_index < n_sub_triangles; ++sub_triangle_index)
    for (unsigned int sub_integral_index = 0;
         sub_integral_index < n_sub_integrals; ++sub_integral_index) {
      // Set integral bounds functions
      this->set_integral_bounds_function_pointers(sub_triangle_index,
                                                  sub_integral_index);
      // With the sub_triangle_index and sub_integral_index set, we now have to
      // evaulate the 4-D integral that we have for this case

      // First, get the upper and lower bounds for the outer integral
      Real psi_lower, psi_upper;
      // this->get_psi_bounds(sub_integral_index, &psi_lower, &psi_upper);
      this->get_current_psi_bounds(psi_lower, psi_upper);

      // Get transformation to work with Gauss quadrature...
      const Real psi_average = (psi_lower + psi_upper) / (Real(2.));
      const Real psi_half_length = (psi_upper - psi_lower) / (Real(2.));
      // Get Gauss quadrature parameters
      std::vector<double> *xgl_psi, *wgl_psi;
      qcollection.get_weights_and_points(this->int_0_quadrature_index, wgl_psi,
                                         xgl_psi);

      for (unsigned int psi_index = 0; psi_index < wgl_psi->size();
           ++psi_index) {
        const Real jac0 = psi_half_length * wgl_psi->operator[](psi_index);
        // Get the value of psi
        const Real psi =
            psi_half_length * xgl_psi->operator[](psi_index) + psi_average;
        const Real tan_psi_less_pi_half = tan(psi - constants<Real>::PI / Real(2.0));
        const Real cos_psi = cos(psi);
        const Real sin_psi = sin(psi);

        // Get the upper and lower bounds for "u", the next integral
        Real u_lower, u_upper;
        this->get_current_u_bounds(tan_psi_less_pi_half, u_lower, u_upper);
        // Get transformation to work with Gauss quadrature...
        const Real u_average = (u_lower + u_upper) / (Real(2.));
        const Real u_half_length = (u_upper - u_lower) / (Real(2.));

        // Get the Gauss quadrature parameters
        std::vector<double> *xgl_u, *wgl_u;
        qcollection.get_weights_and_points(this->int_1_quadrature_index, wgl_u,
                                           xgl_u);

        for (unsigned int u_index = 0; u_index < wgl_u->size(); ++u_index) {
          const Real jac1 = jac0 * u_half_length * wgl_u->operator[](u_index);
          const Real u = u_half_length * xgl_u->operator[](u_index) + u_average;

          // Get the upper and lower bounds for \Lambda
          Real lambda_lower, lambda_upper;
          this->get_current_lambda_bounds(cos_psi, sin_psi, u, lambda_lower,
                                          lambda_upper, psi);
          // Get transformation to work with Gauss quadrature...
          const Real lambda_average =
              (lambda_lower + lambda_upper) / (Real(2.));
          const Real lambda_half_length =
              (lambda_upper - lambda_lower) / (Real(2.));

          // Get the Gauss quadrature parameters
          std::vector<double> *xgl_lambda, *wgl_lambda;
          qcollection.get_weights_and_points(this->int_2_quadrature_index,
                                             wgl_lambda, xgl_lambda);

          for (unsigned int lambda_index = 0; lambda_index < wgl_lambda->size();
               ++lambda_index) {
            const Real jac2 = jac1 * lambda_half_length *
                              wgl_lambda->operator[](lambda_index);
            const Real lambda =
                lambda_half_length * xgl_lambda->operator[](lambda_index) +
                lambda_average;

            Point<2, Real> uv_test_unmapped = {u, lambda * sin_psi - Real(1.)};
            // Convert uv_test
            const auto uv_test = this->map_u_and_v_coords(uv_test_unmapped);

            // Now that we've computed uv_test, we can precompute the values of
            // the shape functions
            ComponentsReturnType dof_test_memoize[mask_test.size()];
            for (unsigned int dof_test_index = 0;
                 dof_test_index < mask_test.size(); ++dof_test_index) {
              if (!test_has_integrals_to_perform[dof_test_index])
                continue;
              dof_test_memoize[dof_test_index] =
                  testComponent(cell_test, dof_test_index, uv_test);
            }
            const auto unitary_test_u =
                cell_test_geom->unitary_vector(uv_test, u_dir);
            const auto unitary_test_v =
                cell_test_geom->unitary_vector(uv_test, v_dir);
            const auto r_test = cell_test_geom->r(uv_test);

            // Get the upper and lower bounds for \rho
            const Real rho_lower = Real(0.0);
            const Real rho_upper = lambda;

            const Real rho_average = (rho_upper) / (Real(2.));
            const Real rho_half_length = (rho_upper) / (Real(2.));

            // Get the Gauss quadrature parameters
            std::vector<double> *xgl_rho, *wgl_rho;
            qcollection.get_weights_and_points(this->int_3_quadrature_index,
                                               wgl_rho, xgl_rho);
            for (unsigned int dof_test_index = 0;
                 dof_test_index < mask_test.size(); ++dof_test_index) {
              if (!test_has_integrals_to_perform[dof_test_index])
                continue;
              const DoFBase<cell_test.cell_dim> &dof_test =
                  cell_test.get_dof(dof_test_index);
              // const auto test_component_value = testComponent(cell_test,
              // dof_test_index, uv_test);
              const auto &test_component_value =
                  dof_test_memoize[dof_test_index];
              for (unsigned int dof_trial_index = 0;
                   dof_trial_index < mask_trial.size(); ++dof_trial_index) {
                if (!mask_trial[dof_trial_index])
                  continue;
                const DoFBase<cell_trial.cell_dim> &dof_trial =
                    cell_trial.get_dof(dof_trial_index);
                // Now, we evaluate the integrand and then multiply by the
                // integration weights...

                auto integrand_function = [&](const Real& rho){
                  const Real transformation_jacobian = rho * sin_psi*jac2;
                  Point<2, Real> uv_trial_unmapped = {
                      uv_test_unmapped(0) + rho * cos_psi,
                      -rho * sin_psi + uv_test_unmapped(1)};
                  // Convert uv_trial
                  const auto uv_trial =
                      this->map_u_and_v_coords(uv_trial_unmapped);

                  const auto unitary_trial_u =
                      cell_trial_geom->unitary_vector(uv_trial, u_dir);
                  const auto unitary_trial_v =
                      cell_trial_geom->unitary_vector(uv_trial, v_dir);

                  const auto r_trial = cell_trial_geom->r(uv_trial);
                  const Real R = (r_test - r_trial).norm();
                  const CoefficientType kernel_value = kernel(R);

                  const auto &trial_component_value = trialComponent(cell_trial, dof_trial_index, uv_trial);

                  return
                      transformation_jacobian *
                      integrand(test_component_value, trial_component_value,
                                dof_test.direction == u_dir ? unitary_test_u
                                                            : unitary_test_v,
                                dof_trial.direction == u_dir ? unitary_trial_u
                                                             : unitary_trial_v,
                                kernel_value, R);

                };
                Real error = 0;
                CoefficientType value_to_add = boost::math::quadrature::gauss_kronrod<Real, 15>::integrate(integrand_function, rho_lower, rho_upper, 5, 1e-12, &error);

                temp_matrix.accumulate(dof_test_index, dof_trial_index,
                                       value_to_add);
              }
            }
          }
        }
      }
    }
  temp_matrix.apply_carry();
  // Now sync these integrals with the main matrix
  for (unsigned int i = 0; i < mask_test.size(); ++i) {
    if (!mask_test[i])
      continue;
    const auto &dof_test = cell_test.get_dof(i);
    for (unsigned int j = 0; j < mask_trial.size(); ++j) {
      if (!mask_trial[j])
        continue;
      const auto &dof_trial = cell_trial.get_dof(j);
      output->get_is_filled_at(dof_test.global_index, dof_trial.global_index) = true;
      output->at(dof_test.global_index, dof_trial.global_index) +=
          temp_matrix(i, j);
    }
  }
  //  uv_out.close();
}

} // namespace NumericalIntegration
DROMON_NAMESPACE_CLOSE
#endif // DROMON_DIRECTFN_SINGULAR_MOD_H