//
// Created by Jake J. Harmon (jake.harmon@ieee.org) on 5/1/22.
//

#ifndef DROMON_POSTPROCESSING_H
#define DROMON_POSTPROCESSING_H
#include "config.h"
#include "Kernels.h"
#include "DoFHandler.h"
#include "Excitations.h"
#include "QuadratureCollection.h"

DROMON_NAMESPACE_OPEN

namespace PostProcessing
{

template <unsigned int spacedim, class Real, class CoefficientType, class Excitation>
Point<spacedim, CoefficientType> compute_E_scattered(DoFHandler<2, spacedim, CoefficientType, Real>* dof_handler, const std::vector<CoefficientType>& solution_coefficients,
                                                     const MaterialData<Real>& materialData, const Excitation &excitation, const Real& theta_sc, const Real& phi_sc, const unsigned int& ngl, const double& R_dist_scalar = 200000.0)
{
  Quadrature::QuadratureCollection<Quadrature::GaussQuadrature> qcollection;
  qcollection.push_back(ngl);

  const auto& omega = excitation.omega;

  Point<spacedim, CoefficientType> out;
  const unsigned int dim = 2;
  const Real wavelength = 2*constants<Real>::PI/excitation.gamma.imag();
  const Real R_dist_field = R_dist_scalar*wavelength;
  Point<spacedim, Real> r_field = {R_dist_field*sin(theta_sc)*cos(phi_sc), R_dist_field*sin(theta_sc)*sin(phi_sc), R_dist_field*cos(theta_sc)};


  // std::cout << "wavelength: " << wavelength << std::endl;
  //   std::cout << "R_dist_scalar: " << R_dist_scalar << std::endl;
  // std::cout << "R_dist_field: " << R_dist_field << std::endl;
  // std::cout << "theta_sc: " << theta_sc << " phi_sc: " << phi_sc << std::endl;
  //std::cout << "r_field: " << r_field << std::endl;

  // Convert spherical to cartesian
  // We first must compute the scattered electric field
  // E^sc = (-j\omega\mu*\int_S J_S g + k^-2*nabla_S\cdot J_s \nabla g dS + int_S M_S\times\nabla g dS)
  // The RCS is then approximated by 4*pi*r^2*|E^sc|^2/|E^inc|^2

  // Loop through all the cells, integrate over the shape functions
  // Depending on whether the shape function is for an electric DoF or a Magnetic one,
  // integrate the correct component for computing the RCS

  // Create the vector where we will accumulate the contributions from every DoF
  Point<spacedim, CoefficientType> contributions[solution_coefficients.size()];

  for (auto& cell_trial : (dof_handler)->get_dof_parents())
  {
    std::vector<double> *xgl, *wgl;
    qcollection.get_weights_and_points(ngl, wgl, xgl);
    const std::complex<Real> complexj = {0., 1.};
    const unsigned int n_dofs = cell_trial.n_dofs();
    auto cell_trial_geom = cell_trial.get_cell_on_mesh();

    const auto exterior_material =
        (materialData.get_material_domain(
             cell_trial_geom->material_domain_id()))
            .get_exterior();
    const auto& eps = exterior_material.eps;
    const auto& mu = exterior_material.mu;
    const auto& epsr = exterior_material.epsr;
    const auto& mur = exterior_material.mur;
    const CoefficientType wavenumber = omega*sqrt(eps*mu);

    for (unsigned int i_trial = 0; i_trial < wgl->size(); ++i_trial) {
      for (unsigned int j_trial = 0; j_trial < wgl->size(); ++j_trial) {
        
        const Point<2, Real> uv_trial(xgl->at(i_trial), xgl->at(j_trial));
        // Now we must loop through all the degrees of freedom
        // For those degrees of freedom that are masked, we skip them
        auto unitary_trial_u = cell_trial_geom->unitary_vector(uv_trial, u_dir);
        auto unitary_trial_v = cell_trial_geom->unitary_vector(uv_trial, v_dir);

        const auto r_trial = cell_trial_geom->r(uv_trial);
        const auto R_vec = r_field - r_trial;
        const Real wavelength = 2*constants<Real>::PI/excitation.gamma.imag();


        const auto R = R_vec.norm();
        // We need both g and grad_g to compute the scattered field
        // In the case of a very large R, simplifications may be made (i.e.,
        // neglecting the R^2 terms etc.), but since the amount of time is not significant,
        // we compute the full equation
        const auto g = kernels::g(R, omega, epsr, mur);
        const auto g_prime = kernels::grad_g(R_vec, R, omega, epsr, mur);

        // const auto g = kernels::gScaled(R, omega, epsr, mur, wavelength);
        // const auto g_prime = kernels::grad_gScaled(R_vec, omega, epsr, mur, wavelength);

        const Real jacobian = wgl->at(i_trial) * wgl->at(j_trial);
        for (unsigned int dof_trial_index = 0; dof_trial_index < n_dofs;
             ++dof_trial_index)
        {
          // Turn the index into the DoF object
          const DoFBase<dim> &dof_trial = cell_trial.get_dof(dof_trial_index);
          // If the DoF is not active, continue
          if (!dof_trial.is_active)
            continue;
          const auto dof_value = dof_trial.evaluate_shape_function(uv_trial);
          const auto& dof_vector_value = dof_trial.direction == u_dir ? unitary_trial_u : unitary_trial_v;
          // Now, based on whether we have an Electric or Magntic DoF, we must decide
          // how to accumulate the contributions
          if (dof_trial.cur_type == Electric) {
            // In this case, we must also compute the divergence
            Real dof_div_value =
                dof_trial.evaluate_shape_function_divergence(uv_trial);
            // First we have -j\omega\mu*J_S*g and then we have
            // the second part of k^-2\nabla_S\cdotJ_s \nabla g
            const auto temp = jacobian*(-complexj*omega*mu*dof_value*dof_vector_value*g + g_prime*(dof_div_value/(wavenumber*wavenumber)));
            contributions[dof_trial.active_index] += jacobian*(-complexj*omega*mu*dof_value*dof_vector_value*g + g_prime*(dof_div_value/(wavenumber*wavenumber)));
          }
          else if (dof_trial.cur_type == Magnetic)
          {
            // The term to compute here is M_s\times\nabla g
            contributions[dof_trial.active_index] += jacobian*dof_vector_value.dot(g_prime)*dof_value;
          }
        }
      }
    }
  }
  // Now that we've computed all the contributions, to actually evaluate the scattered electric field
  // all we need to do is multiply by the solution coefficients at each entry and then sum the results up
  for (unsigned int i = 0; i < solution_coefficients.size(); ++i)
  {
    out += contributions[i]*solution_coefficients[i];
  }
  return out;
}

template <unsigned int spacedim, class Real, class CoefficientType, class Excitation>
Real compute_RCS(DoFHandler<2, spacedim, CoefficientType, Real>* dof_handler, const std::vector<CoefficientType>& solution_coefficients,
                 const MaterialData<Real>& materialData, const Excitation &excitation, const Real& theta_sc, const Real& phi_sc, const unsigned int& ngl, const double& R_dist_scalar = 200000.0)
{

  const Real wavelength = 2*constants<Real>::PI/excitation.gamma.imag();
  const Real R_dist_field = R_dist_scalar;//*wavelength;


  const auto Esc = compute_E_scattered(dof_handler, solution_coefficients, materialData, excitation, theta_sc, phi_sc, ngl);
  const auto Esc_norm = Esc.norm();

  return Esc_norm*Esc_norm*4.0*constants<Real>::PI*R_dist_field*R_dist_field/wavelength/wavelength;
}

}

DROMON_NAMESPACE_CLOSE
#endif // DROMON_POSTPROCESSING_H