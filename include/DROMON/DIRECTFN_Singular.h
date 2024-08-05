//
// Created by Jake J. Harmon (jake.harmon@ieee.org) on 4/15/22.
//

#ifndef DROMON_DIRECTFN_SINGULAR_H
#define DROMON_DIRECTFN_SINGULAR_H

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
        template<class DoFCellType, class CoefficientType = std::complex<double>,
                class Real = double>
        class DIRECTFN_ST {
        public:
            DIRECTFN_ST() {
                assert(false && "The general implementation of the DIRECTFN_ST numerical "
                                "integrator is not yet implemented!");
            }

        private:
        };

        template<class CoefficientType, class Real>
        class DIRECTFN_ST<DoFParent<2, 3, Real>, CoefficientType, Real> {
        public:
            /**
             *
             * @param N0 - Order of numerical integration for outermost integral
             * @param N1 - See above, for next integral
             * @param N2 - See above, for next integral
             * @param N3 - See above, for innermost integral
             */
            DIRECTFN_ST(const unsigned int &N0, const unsigned int &N1,
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
            template<class Integrand, class Kernel>
            void integrate_singular(Integrand integrand, Kernel kernel,
                                    const DoFParent<2, 3, Real> &cell_test,
                                    const DoFParent<2, 3, Real> &cell_trial,
                                    const DoFMask &mask_test, const DoFMask &mask_trial,
                                    DenseSubMatrix<CoefficientType> *output);

        private:
            void get_psi_bounds(const unsigned int &sub_integral_index, Real *psi_lower,
                                Real *psi_upper) const;

            void get_u_bounds(const unsigned int &sub_integral_index, const Real &psi,
                              Real *u_lower, Real *u_upper) const;

            void get_lambda_bounds(const unsigned int &sub_integral_index,
                                   const Real &psi, const Real &u, Real *lambda_lower,
                                   Real *lambda_upper) const;

            Point<2, Real> map_u_and_v_coords(const unsigned int &sub_triangle_index,
                                              const Point<2, Real> &uv) const;

            const unsigned int n_sub_triangles = 4;
            const unsigned int n_sub_integrals = 6;

            unsigned int int_0_quadrature_index;
            unsigned int int_1_quadrature_index;
            unsigned int int_2_quadrature_index;
            unsigned int int_3_quadrature_index;

            Quadrature::QuadratureCollection<Quadrature::GaussQuadrature> qcollection;
        };

        template<class DoFCellType, class CoefficientType = std::complex<double>,
                class Real = double>
        class DIRECTFN_ET {
        public:
            DIRECTFN_ET() {
                assert(false && "The general implementation of the DIRECTFN_ST numerical "
                                "integrator is not yet implemented!");
            }

        private:
        };

        template<class CoefficientType, class Real>
        class DIRECTFN_ET<DoFParent<2, 3, Real>, CoefficientType, Real> {
        public:
            /**
             *
             * @param N0 - Order of numerical integration for outermost integral
             * @param N1 - See above, for next integral
             * @param N2 - See above, for next integral
             * @param N3 - See above, for innermost integral
             */
            DIRECTFN_ET(const unsigned int &N0, const unsigned int &N1,
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
            template<class Integrand, class Kernel>
            void integrate_singular(
                    Integrand integrand, Kernel kernel, const DoFParent<2, 3, Real> &cell_test,
                    const DoFParent<2, 3, Real> &cell_trial, const DoFMask &mask_test,
                    const DoFMask &mask_trial,
                    const MultiIndex<2, unsigned int> &local_adjacent_edge_indices,
                    DenseSubMatrix<CoefficientType> *output);

        private:
            std::ofstream output;

            void get_theta_bounds(const unsigned int &sub_integral_index,
                                  Real *theta_lower, Real *theta_upper) const;

            void get_psi_bounds(const unsigned int &sub_integral_index, const Real &theta,
                                Real *psi_lower, Real *psi_upper) const;

            void get_u_bounds(const unsigned int &sub_integral_index,
                              const bool &is_lower_region, const Real &theta,
                              const Real &psi, Real *u_lower, Real *u_upper) const;

            void get_lambda_bounds(const unsigned int &sub_integral_index,
                                   const unsigned int &part_index, const Real &theta,
                                   const Real &psi, const Real &u, Real *lambda_lower,
                                   Real *lambda_upper) const;

            void map_u_and_v_coords(const unsigned int &local_singular_edge_index,
                                    Point<2, Real> *uv) const;

            const bool has_double[8] = {false, true, true, true, true, true, false, true};
            const bool is_psi_minus_region[8][2] = {
                    {false, false},
                    {true,  false},
                    {true,  true},
                    {true,  false},
                    {true,  true},
                    {false, true},
                    {false, false},
                    {true,  false}};
            unsigned int int_0_quadrature_index;
            unsigned int int_1_quadrature_index;
            unsigned int int_2_quadrature_index;
            unsigned int int_3_quadrature_index;

            Quadrature::QuadratureCollection<Quadrature::GaussQuadrature> qcollection;

            const unsigned int n_sub_integrals = 8;
        };

        template<class CoefficientType, class Real>

        void DIRECTFN_ET<DoFParent<2, 3, Real>, CoefficientType,
                Real>::map_u_and_v_coords(const unsigned int
                                          &local_singular_edge_index,
                                          Point<2, Real> *uv) const {
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
            if (local_singular_edge_index == 1) {
                Real temp = (*uv)(0);
                (*uv)[0] = -(*uv)(1);
                (*uv)[1] = temp;
            } else if (local_singular_edge_index == 2) {
                Real temp = (*uv)(0);
                (*uv)[0] = (*uv)(1);
                (*uv)[1] = -temp;
            } else if (local_singular_edge_index == 3) {
                (*uv)[0] *= -1;
                (*uv)[1] *= -1;
            }
            // else if (local_singular_edge_index == 0) do nothing;
        }

        template<class CoefficientType, class Real>
        void DIRECTFN_ET<DoFParent<2, 3, Real>, CoefficientType, Real>::
        get_lambda_bounds(const unsigned int &sub_integral_index,
                          const unsigned int &part_index, const Real &theta,
                          const Real &psi, const Real &u, Real *lambda_lower,
                          Real *lambda_upper) const {
            *lambda_lower = Real(0.);
            if (sub_integral_index == 0) {
                *lambda_upper = (u + Real(1.)) / (cos(theta) * cos(psi));
                return;
            } else if (sub_integral_index == 6) {
                *lambda_upper = (u - Real(1.)) / (cos(theta) * cos(psi));
                return;
            } else {
                assert(part_index == 0 ||
                       part_index == 1 && "Invalid part_index! Must be either 0 or 1!");
                if (part_index == 0)
                    switch (sub_integral_index) {
                        case 1:
                        case 2:
                        case 3: {
                            const Real test = cos(theta);
                            const Real test2 = cos(psi);
                            *lambda_upper = (u + Real(1.0)) / (cos(theta) * cos(psi));
                        }
                            break;
                        case 4:
                            *lambda_upper = Real(2.) / (sin(theta) * cos(psi));
                            break;
                        case 5:
                            *lambda_upper = Real(2.) / sin(psi);
                            break;
                        case 7:
                            *lambda_upper = (u - Real(1.)) / (cos(theta) * cos(psi));
                            break;
                        default:
                            assert(false && "Sub_integra_index is invalid for this region!");
                    }
                else
                    switch (sub_integral_index) {
                        case 1:
                        case 3:
                        case 7:
                            *lambda_upper = Real(2.) / sin(psi);
                            break;
                        case 2:
                            *lambda_upper = Real(2.) / (sin(theta) * cos(psi));
                            break;
                        case 4:
                        case 5:
                            *lambda_upper = (u - Real(1.)) / (cos(theta) * cos(psi));
                            break;
                        default:
                            assert(false && "Sub_integral_index is invalid for this region!");
                    }
            }
        }

        template<class CoefficientType, class Real>
        void DIRECTFN_ET<DoFParent<2, 3, Real>, CoefficientType, Real>::get_u_bounds(
                const unsigned int &sub_integral_index, const bool &is_lower_region,
                const Real &theta, const Real &psi, Real *u_lower, Real *u_upper) const {
            if (sub_integral_index == 0 || sub_integral_index == 6) {
                *u_lower = Real(-1.);
                *u_upper = Real(1.);
            } else {
                if (is_lower_region) {
                    *u_lower = Real(-1.);
                    switch (sub_integral_index) {
                        case 1:
                            *u_upper = Real(2.) * cos(theta) / tan(psi) - Real(1.);
                            break;
                        case 2:
                            *u_upper = Real(2.0) * tan(constants<Real>::PI / Real(2.0) - theta) -
                                       Real(1.0); // u_1_theta;
                            break;
                        case 3:
                            *u_upper = Real(2.) * cos(theta) / tan(psi) - Real(1.); // u_1_theta
                            break;
                        case 4:
                            *u_upper = Real(2.) / tan(theta) + Real(1.); // u_2_theta;
                            break;
                        case 5:
                        case 7:
                            *u_upper = Real(2.) * cos(theta) / tan(psi) + Real(1.); // u_2_psi;
                            break;
                        default:
                            assert(false && "Invalid sub_integral_index for specified region!");
                    }
                } else {
                    *u_upper = Real(1.);
                    switch (sub_integral_index) {
                        case 1:
                            *u_lower = Real(2.) * cos(theta) / tan(psi) - Real(1.);
                            break;
                        case 2:
                            *u_lower = Real(2.) / tan(theta) - Real(1.); // u_1_theta;
                            break;
                        case 3:
                            *u_lower = Real(2.) * cos(theta) / tan(psi) - Real(1.); // u_1_psi;
                            break;
                        case 4:
                            *u_lower = Real(2.) / tan(theta) + Real(1.);
                            break;
                        case 5:
                        case 7:
                            *u_lower = Real(2.) * cos(theta) / tan(psi) + Real(1.);
                            break;
                        default:
                            assert(false && "Invalid sub_integral_index for specified region!");
                    }
                }
            }
        }

        template<class CoefficientType, class Real>
        void DIRECTFN_ET<DoFParent<2, 3, Real>, CoefficientType, Real>::get_psi_bounds(
                const unsigned int &sub_integral_index, const Real &theta, Real *psi_lower,
                Real *psi_upper) const {
            switch (sub_integral_index) {
                case 0:
                    *psi_lower = Real(0.);
                    *psi_upper = atan(cos(theta));
                    break;
                case 1:
                    *psi_lower = atan(cos(theta));
                    *psi_upper = constants<Real>::PI / Real(2.);
                    break;
                case 2:
                    *psi_lower = Real(0.);
                    *psi_upper = atan(sin(theta));
                    break;
                case 3:
                    *psi_lower = atan(sin(theta));
                    *psi_upper = constants<Real>::PI / Real(2.);
                    break;
                case 4:
                    *psi_lower = Real(0.);
                    *psi_upper = atan(sin(theta));
                    break;
                case 5:
                    *psi_lower = atan(sin(theta));
                    *psi_upper = constants<Real>::PI / Real(2.);
                    break;
                case 6:
                    *psi_lower = Real(0.);
                    *psi_upper = -atan(cos(theta));
                    break;
                case 7:
                    *psi_lower = -atan(cos(theta));
                    *psi_upper = constants<Real>::PI / Real(2.);
                    break;
                default:
                    assert(false && "Sub_integral_index is invalid!");
            }
        }

        template<class CoefficientType, class Real>
        void DIRECTFN_ET<DoFParent<2, 3, Real>, CoefficientType,
                Real>::get_theta_bounds(const unsigned int &sub_integral_index,
                                        Real *theta_lower,
                                        Real *theta_upper) const {
            switch (sub_integral_index) {
                case 0:
                case 1:
                    *theta_lower = Real(0.);
                    *theta_upper = Real(constants<Real>::PI) / Real(4.);
                    break;
                case 2:
                case 3:
                    *theta_lower = Real(constants<Real>::PI) / Real(4.);
                    *theta_upper = Real(constants<Real>::PI) / Real(2.);
                    break;
                case 4:
                case 5:
                    *theta_lower = Real(constants<Real>::PI) / Real(2.);
                    *theta_upper = Real(3.) * constants<Real>::PI / Real(4.);
                    break;
                case 6:
                case 7:
                    *theta_lower = Real(3.) * constants<Real>::PI / Real(4.);
                    *theta_upper = constants<Real>::PI;
                    break;
            }
        }

        template<class CoefficientType, class Real>
        DIRECTFN_ET<DoFParent<2, 3, Real>, CoefficientType, Real>::DIRECTFN_ET(
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

        template<class CoefficientType, class Real>
        template<class Integrand, class Kernel>
        void DIRECTFN_ST<DoFParent<2, 3, Real>, CoefficientType, Real>::
        integrate_singular(Integrand integrand, Kernel kernel,
                           const DoFParent<2, 3, Real> &cell_test,
                           const DoFParent<2, 3, Real> &cell_trial,
                           const DoFMask &mask_test, const DoFMask &mask_trial,
                           DenseSubMatrix<CoefficientType> *output) {
          auto cell_test_geom = cell_test.get_cell_on_mesh();
          auto cell_trial_geom = cell_trial.get_cell_on_mesh();
  //          unsigned int counter = 0;
//            std::ofstream uv_out("uv_points_out_ST.txt");
//            uv_out << std::setprecision(16);
            std::complex<Real> complexj = {0., 1.};

            // For the self term, we have integrations over 4 subtriangles
            // These subtriangles then decompose into 6 4-D integrals
            for (unsigned int sub_triangle_index = 0;
                 sub_triangle_index < n_sub_triangles; ++sub_triangle_index)
                for (unsigned int sub_integral_index = 0;
                     sub_integral_index < n_sub_integrals; ++sub_integral_index) {
                    // With the sub_triangle_index and sub_integral_index set, we now have to
                    // evaulate the 4-D integral that we have for this case

                    // First, get the upper and lower bounds for the outer integral
                    Real psi_lower, psi_upper;
                    this->get_psi_bounds(sub_integral_index, &psi_lower, &psi_upper);
                    // Get transformation to work with Gauss quadrature...
                    const Real psi_average = (psi_lower + psi_upper) / (Real(2.));
                    const Real psi_half_length = (psi_upper - psi_lower) / (Real(2.));
                    // Get Gauss quadrature parameters
                    std::vector<double> *xgl_psi, *wgl_psi;
                    qcollection.get_weights_and_points(this->int_0_quadrature_index, wgl_psi,
                                                       xgl_psi);

                    for (unsigned int psi_index = 0; psi_index < wgl_psi->size();
                         ++psi_index) {
                        // Get the value of psi
                        const Real psi = psi_half_length * xgl_psi->operator[](psi_index) + psi_average;

                        // Get the upper and lower bounds for "u", the next integral
                        Real u_lower, u_upper;
                        this->get_u_bounds(sub_integral_index, psi, &u_lower, &u_upper);
                        // Get transformation to work with Gauss quadrature...
                        const Real u_average = (u_lower + u_upper) / (Real(2.));
                        const Real u_half_length = (u_upper - u_lower) / (Real(2.));

                        // Get the Gauss quadrature parameters
                        std::vector<double> *xgl_u, *wgl_u;
                        qcollection.get_weights_and_points(this->int_1_quadrature_index, wgl_u,
                                                           xgl_u);

                        for (unsigned int u_index = 0; u_index < wgl_u->size(); ++u_index) {
                          if (sub_triangle_index == 1 && sub_integral_index == 4 && psi_index == 0 && u_index == 8)
                          {
                            int pause = 0;
                          }
                            const Real u = u_half_length * xgl_u->operator[](u_index) + u_average;

                            // Get the upper and lower bounds for \Lambda
                            Real lambda_lower, lambda_upper;
                            this->get_lambda_bounds(sub_integral_index, psi, u, &lambda_lower,
                                                    &lambda_upper);
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
                              if (sub_triangle_index == 1 && sub_integral_index == 4 && psi_index == 0 && u_index == 8 && lambda_index == 8)
                              {
                                int pause = 0;
                              }
                                const Real lambda =
                                        lambda_half_length * xgl_lambda->operator[](lambda_index) +
                                        lambda_average;

                                Point<2, Real> uv_test_unmapped = {u, lambda * sin(psi) - Real(1.)};
                                // Convert uv_test
                                const auto uv_test =
                                        this->map_u_and_v_coords(sub_triangle_index, uv_test_unmapped);
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
                                            rho_half_length * xgl_rho->operator[](rho_index) + rho_average;
                                    auto test_cos_psi = cos(psi);
                                    auto test_sin_psi = sin(psi);
                                    const Real transformation_jacobian = rho * sin(psi);
                                    Point<2, Real> uv_trial_unmapped = {
                                            uv_test_unmapped(0) + rho * cos(psi),
                                            -rho * sin(psi) + uv_test_unmapped(1)};
                                    // Convert uv_trial
                                    const auto uv_trial = this->map_u_and_v_coords(sub_triangle_index,
                                                                                   uv_trial_unmapped);

                                    const auto unitary_trial_u =
                                        cell_trial_geom->unitary_vector(uv_trial, u_dir);
                                    const auto unitary_trial_v =
                                        cell_trial_geom->unitary_vector(uv_trial, v_dir);

                                    const auto r_trial = cell_trial_geom->r(uv_trial);
                                    const Real R = (r_test - r_trial).norm();
                                    const CoefficientType kernel_value = kernel(R);
                                    for (unsigned int dof_test_index = 0;
                                         dof_test_index < mask_test.size(); ++dof_test_index) {
                                        if (!mask_test[dof_test_index])
                                            continue;
                                        const DoFBase<cell_test.cell_dim> &dof_test =
                                                cell_test.get_dof(dof_test_index);

                                        for (unsigned int dof_trial_index = 0;
                                             dof_trial_index < mask_trial.size(); ++dof_trial_index) {
                                            if (!mask_trial[dof_trial_index])
                                                continue;
                                            const DoFBase<cell_trial.cell_dim> &dof_trial =
                                                    cell_trial.get_dof(dof_trial_index);
                                            std::cout << "this does happen in the code\n";
                                            // Now, we evaluate the integrand and then multiply by the
                                            // integration weights...
                                            output->at(dof_test.global_index, dof_trial.global_index) +=
                                                    transformation_jacobian * psi_half_length *
                                                    u_half_length * lambda_half_length * rho_half_length *
                                                    wgl_psi->operator[](psi_index) * wgl_u->operator[](u_index) *
                                                    wgl_lambda->operator[](lambda_index) * wgl_rho->operator[](rho_index) *
                                                    integrand(cell_test, cell_trial, dof_test_index,
                                                              dof_trial_index, uv_test, uv_trial, unitary_test_u, unitary_test_v, unitary_trial_u, unitary_trial_v, kernel_value, R);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
        }

        template<class CoefficientType, class Real>
        template<class Integrand, class Kernel>
        void DIRECTFN_ET<DoFParent<2, 3, Real>, CoefficientType, Real>::
        integrate_singular(
                Integrand integrand, Kernel kernel, const DoFParent<2, 3, Real> &cell_test,
                const DoFParent<2, 3, Real> &cell_trial, const DoFMask &mask_test,
                const DoFMask &mask_trial,
                const MultiIndex<2, unsigned int> &local_adjacent_edge_indices,
                DenseSubMatrix<CoefficientType> *output) {

          // First, we create a temporary Matrix, as this is more efficient than accessing the DenseSubMatrix frequently
          // While the mask might have disabled entries, given that the number of DoFs on cells will be fairly low,
          // the difference is not significant
          ContiguousMatrix<CoefficientType> temp_matrix(mask_test.size(), mask_trial.size());


          auto cell_test_geom = cell_test.get_cell_on_mesh();
          auto cell_trial_geom = cell_trial.get_cell_on_mesh();
            std::complex<Real> complexj = {0., 1.};
            //  std::ofstream uv_out("uv_points_out.txt");
            //  uv_out << std::setprecision(16);
            // We have the 8 primary sub integrals (some of which split into 2 additional
            // integrals)
            unsigned int counter = 0;
            for (unsigned sub_integral_index = 0; sub_integral_index < n_sub_integrals;
                 ++sub_integral_index) {

                // We have the outer integral over \theta
                Real theta_lower, theta_upper;
                this->get_theta_bounds(sub_integral_index, &theta_lower, &theta_upper);
                // Get transformation to work with Gauss quadrature...
                const Real theta_average = (theta_lower + theta_upper) / (Real(2.));
                const Real theta_half_length = (theta_upper - theta_lower) / (Real(2.));
                // Get Gauss quadrature parameters
                std::vector<double> *xgl_theta, *wgl_theta;
                qcollection.get_weights_and_points(this->int_0_quadrature_index, wgl_theta,
                                                   xgl_theta);
                for (unsigned int theta_index = 0; theta_index < wgl_theta->size();
                     ++theta_index) {
                    const Real jac0 = theta_half_length*wgl_theta->operator[](theta_index);
                    // Get the value of psi
                    const Real theta =
                            theta_half_length * xgl_theta->operator[](theta_index) + theta_average;

                    // Now we have the integral over \Psi
                    Real psi_lower, psi_upper;
                    this->get_psi_bounds(sub_integral_index, theta, &psi_lower, &psi_upper);
                    const Real psi_average = (psi_lower + psi_upper) / (Real(2.));
                    const Real psi_half_length = (psi_upper - psi_lower) / (Real(2.));
                    // Get Gauss quadrature parameters
                    std::vector<double> *xgl_psi, *wgl_psi;
                    qcollection.get_weights_and_points(this->int_1_quadrature_index, wgl_psi,
                                                       xgl_psi);

                    for (unsigned int psi_index = 0; psi_index < wgl_psi->size();
                         ++psi_index) {
                        const Real jac1 = psi_half_length*wgl_psi->operator[](psi_index)*jac0;
                        const Real psi = psi_half_length * xgl_psi->operator[](psi_index) + psi_average;

                        // Now we have the integral over u, which, for certain sub integrals,
                        // actually splits into a pair of integrals as determined by is_lower
                        // region...
                        Real u_lower, u_upper;
                        const unsigned int n_parts =
                                this->has_double[sub_integral_index] ? 2 : 1;
                        for (unsigned part_index = 0; part_index < n_parts; ++part_index) {
                            bool is_lower_region = true;
                            if (sub_integral_index == 7) {
                                if (part_index == 0)
                                    is_lower_region = false;
                            } else if (part_index > 0)
                                is_lower_region = false;
                            this->get_u_bounds(sub_integral_index, is_lower_region, theta, psi,
                                               &u_lower, &u_upper);
                            const Real u_average = (u_lower + u_upper) / (Real(2.));
                            const Real u_half_length = (u_upper - u_lower) / (Real(2.));
                            // Get Gauss quadrature parameters
                            std::vector<double> *xgl_u, *wgl_u;
                            qcollection.get_weights_and_points(this->int_2_quadrature_index,
                                                               wgl_u, xgl_u);

                            for (unsigned int u_index = 0; u_index < wgl_u->size(); ++u_index) {
                                const Real jac2 = u_half_length*wgl_u->operator[](u_index)*jac1;
                                const Real u = u_half_length * xgl_u->operator[](u_index) + u_average;
                                Real lambda_lower, lambda_upper;

                                this->get_lambda_bounds(sub_integral_index, part_index, theta, psi,
                                                        u, &lambda_lower, &lambda_upper);
                                const Real lambda_average =
                                        (lambda_lower + lambda_upper) / (Real(2.));
                                const Real lambda_half_length =
                                        (lambda_upper - lambda_lower) / (Real(2.));
                                // Get Gauss quadrature parameters
                                std::vector<double> *xgl_lambda, *wgl_lambda;
                                qcollection.get_weights_and_points(this->int_3_quadrature_index,
                                                                   wgl_lambda, xgl_lambda);
                                for (unsigned int lambda_index = 0;
                                     lambda_index < wgl_lambda->size(); ++lambda_index) {
                                    const Real lambda =
                                            lambda_half_length * xgl_lambda->operator[](lambda_index) +
                                            lambda_average;
                                    // Now we need to get the actual uv coordinates
                                    auto tester = sin(psi);
                                    Point<2, Real> uv_test = {u, lambda * sin(psi) - Real(1.0)};
                                    Point<2, Real> uv_trial = {lambda * cos(psi) * cos(theta) - u,
                                                               lambda * cos(psi) * sin(theta) -
                                                               Real(1.)};
                                    const Real transformation_jacobian = lambda * lambda * cos(psi);
                                    const Real full_jacobian = transformation_jacobian*wgl_lambda->operator[](lambda_index)*lambda_half_length*jac2;
                                    // However, in constructing these points, we have assumed a
                                    // particular orientation
                                    // of the edge adjacency
                                    // Since this adjacency is not satisfied in most cases, we must
                                    // introduce the correct transformation

                                    this->map_u_and_v_coords(local_adjacent_edge_indices(0),
                                                             &uv_test);
                                    this->map_u_and_v_coords(local_adjacent_edge_indices(1),
                                                             &uv_trial);

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
                                    //  uv_out << uv_test(0) << " " << uv_test(1) << " " <<
                                    //  uv_trial(0) << " " << uv_trial(1) << std::endl;
                                    for (unsigned int dof_test_index = 0;
                                         dof_test_index < mask_test.size(); ++dof_test_index) {
                                        if (!mask_test[dof_test_index])
                                            continue;
                                        const DoFBase<cell_test.cell_dim> &dof_test =
                                                cell_test.get_dof(dof_test_index);

                                        for (unsigned int dof_trial_index = 0;
                                             dof_trial_index < mask_trial.size(); ++dof_trial_index) {
                                            if (!mask_trial[dof_trial_index])
                                                continue;
                                            const DoFBase<cell_trial.cell_dim> &dof_trial =
                                                    cell_trial.get_dof(dof_trial_index);
                                            // Now, we evaluate the integrand and then multiply by the
                                            // integration weights...
//                                            output->at(dof_test.global_index,
//                                                       dof_trial.global_index) +=
//                                                    transformation_jacobian *
//                                                    theta_half_length * psi_half_length
//                                                    * u_half_length * lambda_half_length
//                                                    * wgl_theta->at(theta_index) *
//                                                    wgl_psi->at(psi_index) *
//                                                    wgl_u->at(u_index) *
//                                                    wgl_lambda->at(lambda_index) *
//                                                    integrand(cell_test, cell_trial,
//                                                              dof_test_index,
//                                                              dof_trial_index, uv_test,
//                                                              uv_trial, unitary_test_u, unitary_test_v, unitary_trial_u, unitary_trial_v, kernel_value);
                                            temp_matrix(dof_test_index, dof_trial_index) += full_jacobian*integrand(cell_test, cell_trial,dof_test_index,dof_trial_index, uv_test,uv_trial, unitary_test_u, unitary_test_v, unitary_trial_u, unitary_trial_v, kernel_value, R);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            for (unsigned int i =0; i < mask_test.size(); ++i)
            {
              if (!mask_test[i])
                continue;
              const auto &dof_test = cell_test.get_dof(i);
              for (unsigned int j = 0; j < mask_trial.size(); ++j)
              {
                if (!mask_trial[j])
                  continue;
                const auto &dof_trial = cell_trial.get_dof(j);
                output->at(dof_test.global_index, dof_trial.global_index) += temp_matrix(i,j);
              }
            }
        }

        template<class CoefficientType, class Real>
        void DIRECTFN_ST<DoFParent<2, 3, Real>, CoefficientType, Real>::get_psi_bounds(
                const unsigned int &sub_integral_index, Real *psi_lower,
                Real *psi_upper) const {
            switch (sub_integral_index) {
                case 0:
                    *psi_lower = 0;
                    *psi_upper = constants<Real>::PI / Real(4.);
                    break;
                case 1:
                case 2:
                    *psi_lower = constants<Real>::PI / Real(4.);
                    *psi_upper = constants<Real>::PI / Real(2.);
                    break;
                case 3:
                case 4:
                    *psi_lower = constants<Real>::PI / Real(2.);
                    *psi_upper = Real(3.) * constants<Real>::PI / Real(4.);
                    break;
                case 5:
                    *psi_lower = Real(3.) * constants<Real>::PI / Real(4.);
                    *psi_upper = constants<Real>::PI;
                    break;
            }
        }

        template<class CoefficientType, class Real>
        void DIRECTFN_ST<DoFParent<2, 3, Real>, CoefficientType, Real>::get_u_bounds(
                const unsigned int &sub_integral_index, const Real &psi, Real *u_lower,
                Real *u_upper) const {
            switch (sub_integral_index) {
                case 0:
                case 5:
                    *u_lower = Real(-1.);
                    *u_upper = Real(1.);
                    break;
                case 1:
                    *u_lower = Real(2.) * tan(psi - constants<Real>::PI / Real(2.)) + 1;
                    *u_upper = Real(1.);
                    break;
                case 2:
                    *u_lower = Real(-1.);
                    *u_upper = Real(2.) * tan(psi - constants<Real>::PI / Real(2.)) + 1;
                    break;
                case 3:
                    *u_lower = Real(2.) * tan(psi - constants<Real>::PI / Real(2.)) - 1;
                    *u_upper = Real(1.);
                    break;
                case 4:
                    *u_lower = Real(-1.);
                    *u_upper = Real(2.) * tan(psi - constants<Real>::PI / Real(2.)) - 1;
                    break;
                default:
                    assert(false && "Sub_integral_index is invalid!");
            }
        }

        template<class CoefficientType, class Real>
        void DIRECTFN_ST<DoFParent<2, 3, Real>, CoefficientType, Real>::
        get_lambda_bounds(const unsigned int &sub_integral_index, const Real &psi,
                          const Real &u, Real *lambda_lower,
                          Real *lambda_upper) const {
            *lambda_lower = Real(0.0);
            switch (sub_integral_index) {
                case 0:
                case 1:
                    *lambda_upper = (Real(1.) - u) / cos(psi);
                    break;
                case 2:
                case 3:
                    *lambda_upper = (Real(2.)) / sin(psi);
                    break;
                case 4:
                case 5:
                 // *lambda_upper = -(u + Real(1.0)) / cos(psi);
                  *lambda_upper = (u + Real(1.0)) / cos(constants<Real>::PI - psi);
                    break;
                default:
                    assert(false && "Sub_integral_index not valid!");
            }
        }

        template<class CoefficientType, class Real>
        Point<2, Real>
        DIRECTFN_ST<DoFParent<2, 3, Real>, CoefficientType, Real>::map_u_and_v_coords(
                const unsigned int &sub_triangle_index, const Point<2, Real> &uv) const {
            Point<2, Real> out;
            switch (sub_triangle_index) {
                case 0:
                    out = uv;
                    break;
                case 1:
                    out[0] = Real(-1.) * uv(1);
                    out[1] = uv(0);
                    break;
                case 2:
                    out[0] = Real(-1.) * uv(0);
                    out[1] = Real(-1.) * uv(1);
                    break;
                case 3:
                    out[0] = uv(1);
                    out[1] = Real(-1.) * uv(0);
                    break;
                default:
                    assert(false && "Invalid sub_triangle_index!");
            }
            return out;
        }

        template<class CoefficientType, class Real>
        DIRECTFN_ST<DoFParent<2, 3, Real>, CoefficientType, Real>::DIRECTFN_ST(
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

        template<class DoFCellType, class CoefficientType = std::complex<double>,
                class Real = double>
        class DIRECTFN_VT {
        public:
            DIRECTFN_VT() {
                assert(false && "The general implementation of the DIRECTFN_VT numerical "
                                "integrator is not yet implemented!");
            }

        private:
        };

        template<class CoefficientType, class Real>
        class DIRECTFN_VT<DoFParent<2, 3, Real>, CoefficientType, Real> {
        public:
            /**
             *
             * @param N0 - Order of numerical integration for outermost integral
             * @param N1 - See above, for next integral
             * @param N2 - See above, for next integral
             * @param N3 - See above, for innermost integral
             */
            DIRECTFN_VT(const unsigned int &N0, const unsigned int &N1,
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
            template<class Integrand, class Kernel>
            void integrate_singular(Integrand integrand, Kernel kernel,
                                    const DoFParent<2, 3, Real> &cell_test,
                                    const DoFParent<2, 3, Real> &cell_trial,
                                    const DoFMask &mask_test, const DoFMask &mask_trial,
                                    const MultiIndex<2, unsigned int> &local_adjacent_vertex_indices,
                                    DenseSubMatrix<CoefficientType> *output);

        private:
            static constexpr unsigned int n_sub_ranges = 4;

            void get_theta_bounds_test(const unsigned int &sub_range_index,
                                       Real *theta_lower, Real *theta_upper) const;

            void get_theta_bounds_trial(const unsigned int &sub_range_index,
                                        Real *theta_lower, Real *theta_upper) const;

            void get_L_test(const unsigned int &sub_range_index, const Real &theta_test,
                            Real *lp) const;

            void get_L_trial(const unsigned int &sub_range_index, const Real &theta_trial,
                             Real *lq) const;

            void get_psi_bounds_I(const Real &atan_Lq_over_Lp, Real *psi_lower, Real *psi_upper);

            void get_psi_bounds_II(const Real &atan_Lq_over_Lp, Real *psi_lower, Real *psi_upper);

            void map_u_and_v_coords(const unsigned int &local_singular_vertex_index,
                                    Point<2, Real> *uv) const;


            unsigned int int_0_quadrature_index;
            unsigned int int_1_quadrature_index;
            unsigned int int_2_quadrature_index;
            unsigned int int_3_quadrature_index;

            Quadrature::QuadratureCollection<Quadrature::GaussQuadrature> qcollection;
        };

        template<class CoefficientType, class Real>
        template<class Integrand, class Kernel>
        void DIRECTFN_VT<DoFParent<2, 3, Real>, CoefficientType, Real>::
        integrate_singular(Integrand integrand, Kernel kernel,
                           const DoFParent<2, 3, Real> &cell_test,
                           const DoFParent<2, 3, Real> &cell_trial,
                           const DoFMask &mask_test, const DoFMask &mask_trial,
                           const MultiIndex<2, unsigned int> &local_adjacent_vertex_indices,
                           DenseSubMatrix<CoefficientType> *output) {
          auto cell_test_geom = cell_test.get_cell_on_mesh();
          auto cell_trial_geom = cell_trial.get_cell_on_mesh();

            for (unsigned int sub_range_index = 0; sub_range_index < this->n_sub_ranges; ++sub_range_index) {

                Real theta_test_lower, theta_test_upper;
                this->get_theta_bounds_test(sub_range_index, &theta_test_lower, &theta_test_upper);
                const Real theta_test_average = (theta_test_lower + theta_test_upper) / (Real(2.));
                const Real theta_test_half_length = (theta_test_upper - theta_test_lower) / (Real(2.));
                // Get Gauss quadrature parameters
                std::vector<double> *xgl_theta_test, *wgl_theta_test;
                qcollection.get_weights_and_points(this->int_0_quadrature_index, wgl_theta_test,
                                                   xgl_theta_test);
                for (unsigned int theta_test_index = 0; theta_test_index < wgl_theta_test->size(); ++theta_test_index) {
                  if (theta_test_index == 2)
                  {
                    int pause = 0;
                  }
                    const Real theta_test =
                            theta_test_half_length * xgl_theta_test->operator[](theta_test_index) + theta_test_average;
                    Real L_test;
                    this->get_L_test(sub_range_index, theta_test, &L_test);
                    Real theta_trial_lower, theta_trial_upper;
                    this->get_theta_bounds_trial(sub_range_index, &theta_trial_lower, &theta_trial_upper);
                    const Real theta_trial_average = (theta_trial_lower + theta_trial_upper) / (Real(2.));
                    const Real theta_trial_half_length = (theta_trial_upper - theta_trial_lower) / (Real(2.));
                    // Get Gauss quadrature parameters
                    std::vector<double> *xgl_theta_trial, *wgl_theta_trial;
                    qcollection.get_weights_and_points(this->int_1_quadrature_index, wgl_theta_trial,
                                                       xgl_theta_trial);
                    for (unsigned int theta_trial_index = 0;
                         theta_trial_index < wgl_theta_trial->size(); ++theta_trial_index) {
                      if (theta_trial_index == 6)
                      {
                        int pause = 0;
                      }
                        const Real theta_trial =
                                theta_trial_half_length * xgl_theta_trial->operator[](theta_trial_index) + theta_trial_average;
                        Real L_trial;
                        this->get_L_trial(sub_range_index, theta_trial, &L_trial);
                        const Real atan_Ltrial_over_Ltest = atan(L_trial / L_test);

                        std::vector<double> *xgl_psi, *wgl_psi;
                        qcollection.get_weights_and_points(this->int_2_quadrature_index, wgl_psi,
                                                           xgl_psi);

                        Real psi_lower, psi_upper, psi_average, psi_half_length;
                        // For the first integral part
                        this->get_psi_bounds_I(atan_Ltrial_over_Ltest, &psi_lower, &psi_upper);
                        psi_average = (psi_upper + psi_lower)/Real(2.0);
                        psi_half_length = (psi_upper - psi_lower)/Real(2.0);
                        for (unsigned int psi_index = 0; psi_index < wgl_psi->size(); ++psi_index)
                        {
                            const Real psi = psi_half_length*xgl_psi->operator[](psi_index) + psi_average;
                            const Real lambda_lower = Real(0.0);
                            const Real lambda_upper = L_test/cos(psi);
                            const Real lambda_average = lambda_upper/Real(2.0);
                            const Real lambda_half_length = lambda_upper/Real(2.0);
                            std::vector<double> *xgl_lambda, *wgl_lambda;
                            qcollection.get_weights_and_points(this->int_3_quadrature_index, wgl_lambda,
                                                               xgl_lambda);
                            for (unsigned int lambda_index = 0; lambda_index < wgl_lambda->size(); ++lambda_index)
                            {
                                const Real lambda = lambda_half_length*(xgl_lambda->operator[](lambda_index) + Real(1.0));

                                Point<2, Real> uv_test = {lambda*cos(psi)*cos(theta_test) - Real(1.0), lambda*cos(psi)*sin(theta_test) - Real(1.0)};
                                Point<2, Real> uv_trial = {lambda*sin(psi)*cos(theta_trial) - Real(1.0), lambda*sin(psi)*sin(theta_trial) - Real(1.0)};

                                const Real transformation_jacobian = lambda * lambda * lambda * cos(psi)*sin(psi);

                                // Map these uv_coordinates to the necessary ones for the given cell orientation...
                                this->map_u_and_v_coords(local_adjacent_vertex_indices(0), &uv_test);
                                this->map_u_and_v_coords(local_adjacent_vertex_indices(1), &uv_trial);
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
                                for (unsigned int dof_test_index = 0;
                                     dof_test_index < mask_test.size(); ++dof_test_index) {
                                    if (!mask_test[dof_test_index])
                                        continue;
                                    const DoFBase<cell_test.cell_dim> &dof_test =
                                            cell_test.get_dof(dof_test_index);

                                    for (unsigned int dof_trial_index = 0;
                                         dof_trial_index < mask_trial.size(); ++dof_trial_index) {
                                        if (!mask_trial[dof_trial_index])
                                            continue;
                                        const DoFBase<cell_trial.cell_dim> &dof_trial =
                                                cell_trial.get_dof(dof_trial_index);
                                        // Now, we evaluate the integrand and then multiply by the
                                        // integration weights...
                                        output->at(dof_test.global_index,
                                                   dof_trial.global_index) +=
                                                transformation_jacobian *
                                               theta_test_half_length*theta_trial_half_length*psi_half_length*lambda_half_length
                                                * wgl_theta_test->operator[](theta_test_index) * wgl_theta_trial->operator[](theta_trial_index) * wgl_psi->operator[](psi_index)
                                                * wgl_lambda->operator[](lambda_index) *
                                                integrand(cell_test, cell_trial,
                                                          dof_test_index,
                                                          dof_trial_index, uv_test,
                                                          uv_trial, unitary_test_u, unitary_test_v, unitary_trial_u, unitary_trial_v, kernel_value, R);
                                    }
                                }
                            }
                        }

                        // For the second integral part
                        this->get_psi_bounds_II(atan_Ltrial_over_Ltest, &psi_lower, &psi_upper);
                        psi_average = (psi_upper + psi_lower)/Real(2.0);
                        psi_half_length = (psi_upper - psi_lower)/Real(2.0);
                        for (unsigned int psi_index = 0; psi_index < wgl_psi->size(); ++psi_index)
                        {

                            const Real psi = psi_half_length*xgl_psi->at(psi_index) + psi_average;
                            const Real lambda_lower = Real(0.0);
                            const Real lambda_upper = L_trial/sin(psi);
                            const Real lambda_average = lambda_upper/Real(2.0);
                            const Real lambda_half_length = lambda_upper/Real(2.0);
                            std::vector<double> *xgl_lambda, *wgl_lambda;
                            qcollection.get_weights_and_points(this->int_3_quadrature_index, wgl_lambda,
                                                               xgl_lambda);
                            for (unsigned int lambda_index = 0; lambda_index < wgl_lambda->size(); ++lambda_index)
                            {

                                const Real lambda = lambda_half_length*(xgl_lambda->at(lambda_index) + Real(1.0));

                                Point<2, Real> uv_test = {lambda*cos(psi)*cos(theta_test) - Real(1.0), lambda*cos(psi)*sin(theta_test) - Real(1.0)};
                                Point<2, Real> uv_trial = {lambda*sin(psi)*cos(theta_trial) - Real(1.0), lambda*sin(psi)*sin(theta_trial) - Real(1.0)};

                                const Real transformation_jacobian = lambda * lambda * lambda * cos(psi)*sin(psi);
                               const Real full_jacobian = transformation_jacobian *
                                        theta_test_half_length*theta_trial_half_length*psi_half_length*lambda_half_length
                                    * wgl_theta_test->at(theta_test_index) * wgl_theta_trial->at(theta_trial_index) * wgl_psi->at(psi_index)
                                    * wgl_lambda->at(lambda_index);

                                // Map these uv_coordinates to the necessary ones for the given cell orientation...
                                this->map_u_and_v_coords(local_adjacent_vertex_indices(0), &uv_test);
                                this->map_u_and_v_coords(local_adjacent_vertex_indices(1), &uv_trial);
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
                                for (unsigned int dof_test_index = 0;
                                     dof_test_index < mask_test.size(); ++dof_test_index) {
                                    if (!mask_test[dof_test_index])
                                        continue;
                                    const DoFBase<cell_test.cell_dim> &dof_test =
                                            cell_test.get_dof(dof_test_index);

                                    for (unsigned int dof_trial_index = 0;
                                         dof_trial_index < mask_trial.size(); ++dof_trial_index) {
                                        if (!mask_trial[dof_trial_index])
                                            continue;
                                        const DoFBase<cell_trial.cell_dim> &dof_trial =
                                                cell_trial.get_dof(dof_trial_index);
                                        // Now, we evaluate the integrand and then multiply by the
                                        // integration weights...
                                        output->at(dof_test.global_index,
                                                   dof_trial.global_index) +=
                                            full_jacobian *
                                            integrand(cell_test, cell_trial,
                                                      dof_test_index,
                                                      dof_trial_index, uv_test,
                                                      uv_trial, unitary_test_u, unitary_test_v, unitary_trial_u, unitary_trial_v, R);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        template<class CoefficientType, class Real>
        void DIRECTFN_VT<DoFParent<2, 3, Real>, CoefficientType,
                Real>::map_u_and_v_coords(const unsigned int
                                          &local_singular_vertex_index,
                                          Point<2, Real> *uv) const
                                          {
                    switch (local_singular_vertex_index)
                    {
                        case 0:
                            return;
                        case 1: {
                          Real temp = (*uv)(0);
                          (*uv)[0] = -(*uv)(1);
                          (*uv)[1] = temp;
                        }
                            break;

                        case 2: {
                          Real temp = (*uv)(0);
                          (*uv)[0] = (*uv)(1);
                          (*uv)[1] = -temp;
                        }
                            break;
                        case 3:
                            (*uv)[0] *= Real(-1.0);
                            (*uv)[1] *= Real(-1.0);
                            break;
                        default:
                          assert(false && "Local_singular_vertex_index not valid!");
                    }

        }

        template<class CoefficientType, class Real>
        void DIRECTFN_VT<DoFParent<2, 3, Real>, CoefficientType,
                Real>::get_psi_bounds_II(const Real &atan_Lq_over_Lp,
                                         Real *psi_lower, Real *psi_upper) {
            *psi_lower = atan_Lq_over_Lp;
            *psi_upper = constants<Real>::PI / Real(2.0);
        }

        template<class CoefficientType, class Real>
        void DIRECTFN_VT<DoFParent<2, 3, Real>, CoefficientType,
                Real>::get_psi_bounds_I(const Real &atan_Lq_over_Lp,
                                        Real *psi_lower, Real *psi_upper) {
            *psi_lower = Real(0.0);
            *psi_upper = atan_Lq_over_Lp;
        }

        template<class CoefficientType, class Real>
        void DIRECTFN_VT<DoFParent<2, 3, Real>, CoefficientType, Real>::get_L_test(
                const unsigned int &sub_range_index, const Real &theta_test,
                Real *lp) const {
            switch (sub_range_index) {
                case 0:
                case 1:
                    *lp = Real(2.0) / cos(theta_test);
                    break;
                case 2:
                case 3:
                    *lp = Real(2.0) / sin(theta_test);
                    break;
                default:
                    assert(false && "Invalid sub_range_index!");
            }
        }

        template<class CoefficientType, class Real>
        void DIRECTFN_VT<DoFParent<2, 3, Real>, CoefficientType, Real>::get_L_trial(
                const unsigned int &sub_range_index, const Real &theta_trial,
                Real *lq) const {
            switch (sub_range_index) {
                case 0:
                case 2:
                    *lq = Real(2.0) / cos(theta_trial);
                    break;
                case 1:
                case 3:
                    *lq = Real(2.0) / sin(theta_trial);
                    break;
                default:
                    assert(false && "Invalid sub_range_index!");
            }
        }

        template<class CoefficientType, class Real>
        void DIRECTFN_VT<DoFParent<2, 3, Real>, CoefficientType, Real>::
        get_theta_bounds_test(const unsigned int &sub_range_index,
                              Real *theta_lower, Real *theta_upper) const {
            switch (sub_range_index) {
                case 0:
                case 1:
                    *theta_lower = Real(0.0);
                    *theta_upper = constants<Real>::PI / Real(4.0);
                    break;
                case 2:
                case 3:
                    *theta_lower = constants<Real>::PI / Real(4.0);
                    *theta_upper = constants<Real>::PI / Real(2.0);
                    break;
            }
        }

        template<class CoefficientType, class Real>
        void DIRECTFN_VT<DoFParent<2, 3, Real>, CoefficientType, Real>::
        get_theta_bounds_trial(const unsigned int &sub_range_index,
                               Real *theta_lower, Real *theta_upper) const {
            switch (sub_range_index) {
                case 0:
                case 2:
                    *theta_lower = Real(0.0);
                    *theta_upper = constants<Real>::PI / Real(4.0);
                    break;
                case 1:
                case 3:
                    *theta_lower = constants<Real>::PI / Real(4.0);
                    *theta_upper = constants<Real>::PI / Real(2.0);
                    break;
            }
        }

        template<class CoefficientType, class Real>
        DIRECTFN_VT<DoFParent<2, 3, Real>, CoefficientType, Real>::DIRECTFN_VT(
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

    } // namespace NumericalIntegration

DROMON_NAMESPACE_CLOSE

#endif // DROMON_DIRECTFN_SINGULAR_H
