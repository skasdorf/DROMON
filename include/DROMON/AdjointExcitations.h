//
// Created by Jake J. Harmon (jake.harmon@ieee.org) on 5/14/22.
//

#ifndef DROMON_ADJOINTEXCITATIONS_H
#define DROMON_ADJOINTEXCITATIONS_H
#include "DoFBase.h"
#include "DoFMask.h"
#include "DoFHandler.h"
#include "Materials.h"
#include "SubMatrix.h"
#include "config.h"
#include "QuadratureCollection.h"
#include <complex>
#include "Kernels.h"
DROMON_NAMESPACE_OPEN

namespace Excitations {
// In contrast to the excitation for the forward problem, which is largely
// agnostic to the properties of the discretization, as the excitation
// for the adjoint problem amounts to processing based on this discretization
// we must include additional information
template <class DoFCellType, class CoefficientType, class Real = double, class MaterialField = double>
struct AdjointExcitation {
  explicit AdjointExcitation(const Real &frequency);
  AdjointExcitation(const Real &frequency,
                    const Material<MaterialField> &exterior_material);

  // const Real freq;
  const std::complex<Real> gamma;
  const Real omega;
  const Real frequency;
  MaterialField wavenumber;
  const Material<MaterialField> exterior_material;
  Real wavelength;


  virtual void fill_excitation(const DoFCellType &cell_test,
                                  const DoFMask &mask_test,
                                  const unsigned int &ngl,
                                  DenseSubVector<CoefficientType> *output);
};
template <class DoFCellType, class CoefficientType, class Real, class MaterialField>
void AdjointExcitation<DoFCellType, CoefficientType, Real,MaterialField>::fill_excitation(
    const DoFCellType &cell_test, const DoFMask &mask_test,
    const unsigned int &ngl, DenseSubVector<CoefficientType> *output) {
  return;
}
template <class DoFCellType, class CoefficientType, class Real, class MaterialField>
AdjointExcitation<DoFCellType, CoefficientType, Real,MaterialField>::AdjointExcitation(
    const Real &frequency)
    : frequency(frequency), omega(2 * constants<Real>::PI * frequency),
      gamma(2 * constants<Real>::PI * frequency * constants<Real>::ROOT_EPS0MU0_ *
            std::complex<Real>(0, 1.0)) {
  wavelength = (2.0) * constants<Real>::PI / gamma.imag();
  std::cout << "wavelength in AE:" << wavelength << std::endl;
  this->wavenumber = omega*sqrt(exterior_material.eps*exterior_material.mu);
}

template <class DoFCellType, class CoefficientType, class Real, class MaterialField>
AdjointExcitation<DoFCellType, CoefficientType, Real,MaterialField>::AdjointExcitation(
    const Real &frequency,
    const Material<MaterialField> &exterior_material)
    : frequency(frequency), omega(2 * constants<double>::PI * frequency),
      gamma(2 * constants<double>::PI * frequency *
            std::sqrt(exterior_material.eps * exterior_material.mu) *
            std::complex<Real>(0, 1.0)), exterior_material(exterior_material) {
  wavelength = (2.0) * constants<double>::PI / gamma.imag();
}

template <class DoFCellType, class CoefficientType, class Real, class MaterialField>
struct AdjointScatteredFieldExcitation
    : public AdjointExcitation<DoFCellType, CoefficientType, Real,MaterialField> {
  explicit AdjointScatteredFieldExcitation() = default;
};

template <class CoefficientType, class Real, class MaterialField>
struct AdjointScatteredFieldExcitation<DoFParent<2, 3, Real>, CoefficientType,
                                       Real, MaterialField>
    : public AdjointExcitation<DoFParent<2, 3, Real>, CoefficientType, Real,MaterialField> {

  AdjointScatteredFieldExcitation(const Real &frequency,
                                           const Real &theta_sc,
                                           const Real &phi_sc,
                                           const Point<3, Real>& isolatation_direction,
                                           const Real &R_dist_field_scalar = 200000.0);

  AdjointScatteredFieldExcitation(const Real &frequency,
                                  const Material<MaterialField> &exterior_material,
                                  const Real &theta_sc, const Real &phi_sc,
                                  const Point<3, Real>& isolatation_direction,
                                  const Real &R_dist_field_scalar = 200000.0);


  void fill_excitation(const DoFParent<2, 3, Real> &cell_test,
                          const DoFMask &mask_test, const unsigned int &ngl,
                          DenseSubVector<CoefficientType> *output) override;

private:
  Real theta_sc;
  Real phi_sc;
  Real R_dist_field;
  Point<3, Real> r_field;
  Point<3, Real> isolation_direction;
};
template <class CoefficientType, class Real, class MaterialField>
void AdjointScatteredFieldExcitation<
    DoFParent<2, 3, Real>, CoefficientType,
    Real,MaterialField>::fill_excitation(const DoFParent<2, 3, Real> &cell_test,
                              const DoFMask &mask_test, const unsigned int &ngl,
                              DenseSubVector<CoefficientType> *output) {
  Quadrature::QuadratureCollection<Quadrature::GaussQuadrature> qcollection;
  qcollection.push_back(ngl);
  std::vector<double> *xgl, *wgl;
  qcollection.get_weights_and_points(ngl, wgl, xgl);
  // For this particular problem, we simply intregrate over this single patch
  // and compute the necessary integrals related to the DoFs
  // These integrals are the exact same as those performed for in PostProcessing
  auto cell_test_geom = cell_test.get_cell_on_mesh();


  const auto& eps = this->exterior_material.eps;
  const auto& mu = this->exterior_material.mu;
  const auto& epsr = this->exterior_material.epsr;
  const auto& mur = this->exterior_material.mur;
  for (unsigned int i_trial = 0; i_trial < wgl->size(); ++i_trial) {
    for (unsigned int j_trial = 0; j_trial < wgl->size(); ++j_trial) {
      const Point<2, Real> uv_test(xgl->at(i_trial), xgl->at(j_trial));
      // Now we must loop through all the degrees of freedom
      // For those degrees of freedom that are masked, we skip them
      auto unitary_trial_u = cell_test_geom->unitary_vector(uv_test, u_dir);
      auto unitary_trial_v = cell_test_geom->unitary_vector(uv_test, v_dir);

      auto r_trial = cell_test_geom->r(uv_test);
      auto R_vec = r_field - r_trial;//*this->wavelength;

      // R_vec = R_vec * this->wavelength;
      const auto R = R_vec.norm();
      // std::cout << "R.norm(): " << R << std::endl;

      // std::cout << "r_field: " << this->wavelength << std::endl;

      // We need both g and grad_g to compute the scattered field
      // In the case of a very large R, simplifications may be made (i.e.,
      // neglecting the R^2 terms etc.), but since the amount of time is not significant,
      // we compute the full equation

      //  const auto g = kernels::gScaled(R, this->omega, epsr, mur, this->wavelength);
      //  const auto g_prime = kernels::grad_gScaled(R_vec, this->omega, epsr, mur, this->wavelength);

      auto g = kernels::g(R, this->omega, epsr, mur);
      auto g_prime = kernels::grad_g(R_vec, R, this->omega, epsr, mur);

      // std::cout << "gamma: " << std::complex<Real>(0.0,1.0)*this->omega*sqrt(epsr*mur)*constants<Real>::ROOT_EPS0MU0_ * this->wavelength << std::endl;


      // auto g1 = kernels::gScaled(R, this->omega, epsr, mur, this->wavelength);
      // std::cout << "gscaled: " << g1 << "   g: " << g << "  g*exp(jlambda): " << g*std::exp(-constants<Real>::complexj*this->wavelength) <<std::endl;


      const Real jacobian = wgl->at(i_trial) * wgl->at(j_trial);
      for (unsigned int dof_test_index = 0; dof_test_index < mask_test.size();
           ++dof_test_index)
      {
        // Turn the index into the DoF object
        const auto &dof_test = cell_test.get_dof(dof_test_index);
        if (!mask_test[dof_test_index])
          continue;
        const auto dof_value = dof_test.evaluate_shape_function(uv_test);
        const auto& dof_vector_value = dof_test.direction == u_dir ? unitary_trial_u : unitary_trial_v;
        // Now, based on whether we have an Electric or Magntic DoF, we must decide
        // how to accumulate the contributions
        if (dof_test.cur_type == Electric) {
          // In this case, we must also compute the divergence
          Real dof_div_value =
              dof_test.evaluate_shape_function_divergence(uv_test);
          // First we have -j\omega\mu*J_S*g and then we have
          // the second part of k^-2\nabla_S\cdotJ_s \nabla g
          //contributions[dof_trial.active_index] += jacobian*(-constants<Real>::complexj*this->omega*mu*dof_value*dof_vector_value*g + g_prime*(dof_div_value/(this->wavenumber*this->wavenumber)));
          
          auto value = jacobian*(-constants<Real>::complexj*this->omega*mu*dof_value*dof_vector_value*g + g_prime*(dof_div_value/(this->wavenumber*this->wavenumber)));

          output->at(dof_test.global_index) += std::conj(this->isolation_direction.dot(value));
          // output->at(dof_test.global_index) += std::conj(this->isolation_direction.dot(value));
        }
        else if (dof_test.cur_type == Magnetic)
        {
          std::cout << "does it ever make it to magnetic??.........................\n";
          // The term to compute here is M_s\times\nabla g
         // contributions[dof_trial.active_index] += jacobian*dof_vector_value.dot(g_prime)*dof_value;
          const auto value = jacobian*dof_vector_value.dot(g_prime)*dof_value;
          output->at(dof_test.global_index) += std::conj(this->isolation_direction.dot(value));
        }
      }
    }
  }
}

template <class CoefficientType, class Real, class MaterialField>
AdjointScatteredFieldExcitation<DoFParent<2, 3, Real>, CoefficientType, Real, MaterialField>::
    AdjointScatteredFieldExcitation(const Real &frequency,
                                    const Material<MaterialField> &exterior_material,
                                    const Real &theta_sc, const Real &phi_sc,
                                    const Point<3, Real>& isolatation_direction,
                                    const Real &R_dist_field_scalar)
    : AdjointExcitation<DoFParent<2, 3, Real>, CoefficientType, Real>(
          frequency, exterior_material) {


    this->R_dist_field = R_dist_field_scalar*this->wavelength;

  this->r_field = {this->R_dist_field*sin(theta_sc)*cos(phi_sc), this->R_dist_field*sin(theta_sc)*sin(phi_sc), this->R_dist_field*cos(theta_sc)};
  this->isolation_direction = isolatation_direction;
}
template <class CoefficientType, class Real, class MaterialField>
AdjointScatteredFieldExcitation<DoFParent<2, 3, Real>, CoefficientType, Real, MaterialField>::
    AdjointScatteredFieldExcitation(const Real &frequency, const Real &theta_sc,
                                    const Real &phi_sc,
                                    const Point<3, Real>& isolatation_direction,
                                    const Real &R_dist_field_scalar)
    : AdjointExcitation<DoFParent<2, 3, Real>, CoefficientType, Real>(
          frequency) {

    this->R_dist_field = R_dist_field_scalar*this->wavelength;

  this->r_field = {this->R_dist_field*sin(theta_sc)*cos(phi_sc), this->R_dist_field*sin(theta_sc)*sin(phi_sc), this->R_dist_field*cos(theta_sc)};
  this->isolation_direction = isolatation_direction;

}


template <class DoFCellType, class CoefficientType, class Real, class MaterialField>
struct AdjointRCSExcitation
    : public AdjointExcitation<DoFCellType, CoefficientType, Real,MaterialField> {
  explicit AdjointRCSExcitation() = default;
};

template <class CoefficientType, class Real, class MaterialField>
struct AdjointRCSExcitation<DoFParent<2, 3, Real>, CoefficientType,
                                       Real, MaterialField>
    : public AdjointExcitation<DoFParent<2, 3, Real>, CoefficientType, Real,MaterialField> {
  explicit AdjointRCSExcitation(const Real &frequency,
                                           const Real &theta_sc,
                                           const Real &phi_sc,
                                           const Point<3, std::complex<Real>>& scattered_field,
                                           const Real &R_dist_field_scalar = 200000.0);
  AdjointRCSExcitation(const Real &frequency,
                                  const Material<MaterialField> &exterior_material,
                                  const Real &theta_sc, const Real &phi_sc,
                                  const Point<3, std::complex<Real>>& scattered_field,
                                  const Real &R_dist_field_scalar = 200000.0);


  void fill_excitation(const DoFParent<2, 3, Real> &cell_test,
                          const DoFMask &mask_test, const unsigned int &ngl,
                          DenseSubVector<CoefficientType> *output) override;

private:
  Real theta_sc;
  Real phi_sc;
  Real R_dist_field;
  Point<3, Real> r_field;
  Point<3, std::complex<Real>> scattered_field;
  Point<3, std::complex<Real>> conj_scattered_field;
};
template <class CoefficientType, class Real, class MaterialField>
void AdjointRCSExcitation<
    DoFParent<2, 3, Real>, CoefficientType,
    Real,MaterialField>::fill_excitation(const DoFParent<2, 3, Real> &cell_test,
                                             const DoFMask &mask_test, const unsigned int &ngl,
                                             DenseSubVector<CoefficientType> *output) {
  Quadrature::QuadratureCollection<Quadrature::GaussQuadrature> qcollection;
  qcollection.push_back(ngl);
  std::vector<double> *xgl, *wgl;
  qcollection.get_weights_and_points(ngl, wgl, xgl);
  // For this particular problem, we simply intregrate over this single patch
  // and compute the necessary integrals related to the DoFs
  // These integrals are the exact same as those performed for in PostProcessing
  auto cell_test_geom = cell_test.get_cell_on_mesh();

  // std::cout << "RCS is being used\n";


  const auto& eps = this->exterior_material.eps;
  const auto& mu = this->exterior_material.mu;
  const auto& epsr = this->exterior_material.epsr;
  const auto& mur = this->exterior_material.mur;
  for (unsigned int i_trial = 0; i_trial < wgl->size(); ++i_trial) {
    for (unsigned int j_trial = 0; j_trial < wgl->size(); ++j_trial) {
      const Point<2, Real> uv_test(xgl->at(i_trial), xgl->at(j_trial));
      // Now we must loop through all the degrees of freedom
      // For those degrees of freedom that are masked, we skip them
      auto unitary_trial_u = cell_test_geom->unitary_vector(uv_test, u_dir);
      auto unitary_trial_v = cell_test_geom->unitary_vector(uv_test, v_dir);

      const auto r_trial = cell_test_geom->r(uv_test);
      const auto R_vec = r_field - r_trial;
      const auto R = R_vec.norm();
      // We need both g and grad_g to compute the scattered field
      // In the case of a very large R, simplifications may be made (i.e.,
      // neglecting the R^2 terms etc.), but since the amount of time is not significant,
      // we compute the full equation

      const auto g = kernels::g(R, this->omega, epsr, mur);
      const auto g_prime = kernels::grad_g(R_vec, R, this->omega, epsr, mur);

      // const auto g = kernels::gScaled(R, this->omega, epsr, mur, this->wavelength);
      // const auto g_prime = kernels::grad_gScaled(R_vec, this->omega, epsr, mur, this->wavelength);

      const Real jacobian = wgl->at(i_trial) * wgl->at(j_trial);
      for (unsigned int dof_test_index = 0; dof_test_index < mask_test.size();
           ++dof_test_index)
      {
        // Turn the index into the DoF object
        const auto &dof_test = cell_test.get_dof(dof_test_index);
        if (!mask_test[dof_test_index])
          continue;
        const auto dof_value = dof_test.evaluate_shape_function(uv_test);
        const auto& dof_vector_value = dof_test.direction == u_dir ? unitary_trial_u : unitary_trial_v;
        // Now, based on whether we have an Electric or Magntic DoF, we must decide
        // how to accumulate the contributions
        if (dof_test.cur_type == Electric) {
          // In this case, we must also compute the divergence
          Real dof_div_value =
              dof_test.evaluate_shape_function_divergence(uv_test);
          // First we have -j\omega\mu*J_S*g and then we have
          // the second part of k^-2\nabla_S\cdotJ_s \nabla g
          //contributions[dof_trial.active_index] += jacobian*(-constants<Real>::complexj*this->omega*mu*dof_value*dof_vector_value*g + g_prime*(dof_div_value/(this->wavenumber*this->wavenumber)));
          const auto value = jacobian*(-constants<Real>::complexj*this->omega*mu*dof_value*dof_vector_value*g + g_prime*(dof_div_value/(this->wavenumber*this->wavenumber)));
       //   output->at(dof_test.global_index) += value.dot(this->isolation_direction);
          output->at(dof_test.global_index) += Real(8.0*constants<Real>::PI*R_dist_field*R_dist_field)*std::conj(conj_scattered_field.dot(value));
         // (4.0D0*PI*RRR**2*CONJG(2.0D0*CONJG(REFERENCE_SCATTERED_FIELD(KOMPK))*(CC_SIE(NUNK+I)*(-Diel_Dom(1)%SIGM)))/(ALAMBDA**2))
        }
        else if (dof_test.cur_type == Magnetic)
        {
          // The term to compute here is M_s\times\nabla g
          // contributions[dof_trial.active_index] += jacobian*dof_vector_value.dot(g_prime)*dof_value;
          const auto value = jacobian*dof_vector_value.dot(g_prime)*dof_value;
          output->at(dof_test.global_index) += Real(8.0*constants<Real>::PI*R_dist_field*R_dist_field)*std::conj(conj_scattered_field.dot(value));
        }
      }
    }
  }
}

template <class CoefficientType, class Real, class MaterialField>
AdjointRCSExcitation<DoFParent<2, 3, Real>, CoefficientType, Real, MaterialField>::
    AdjointRCSExcitation(const Real &frequency,
                                    const Material<MaterialField> &exterior_material,
                                    const Real &theta_sc, const Real &phi_sc,
                                    const Point<3, std::complex<Real>>& scattered_field,
                                    const Real &R_dist_field_scalar)
    : AdjointExcitation<DoFParent<2, 3, Real>, CoefficientType, Real>(
          frequency, exterior_material) {

    this->R_dist_field = R_dist_field_scalar*this->wavelength;

  this->r_field = {this->R_dist_field*sin(theta_sc)*cos(phi_sc), this->R_dist_field*sin(theta_sc)*sin(phi_sc), this->R_dist_field*cos(theta_sc)};
 this->scattered_field = scattered_field;
 this->conj_scattered_field = scattered_field.conj();
}
template <class CoefficientType, class Real, class MaterialField>
AdjointRCSExcitation<DoFParent<2, 3, Real>, CoefficientType, Real, MaterialField>::
    AdjointRCSExcitation(const Real &frequency, const Real &theta_sc,
                                    const Real &phi_sc,
                                    const Point<3, std::complex<Real>>& scattered_field,
                                    const Real &R_dist_field_scalar)
    : AdjointExcitation<DoFParent<2, 3, Real>, CoefficientType, Real>(
          frequency) {

    this->R_dist_field = R_dist_field_scalar*this->wavelength;

  this->r_field = {this->R_dist_field*sin(theta_sc)*cos(phi_sc), this->R_dist_field*sin(theta_sc)*sin(phi_sc), this->R_dist_field*cos(theta_sc)};
  this->scattered_field = scattered_field;
  this->conj_scattered_field = scattered_field.conj();
}


} // namespace Excitations

DROMON_NAMESPACE_CLOSE
#endif // DROMON_ADJOINTEXCITATIONS_H
