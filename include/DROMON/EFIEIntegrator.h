//
// Created by Jake J. Harmon (jake.harmon@ieee.org) on 4/1/22.
//

#ifndef DROMON_EFIEINTEGRATOR_H
#define DROMON_EFIEINTEGRATOR_H

#include "DIRECTFN_Singular.h"
#include "DIRECTFN_ST_Singular.h"
#include "DIRECTFN_ET_Singular.h"
#include "DIRECTFN_VT_Singular.h"
#include "DoFBase.h"
#include "DoFMask.h"
#include "IntegratorBase.h"
#include "Kernels.h"
#include "Materials.h"
#include "MultiIndex.h"
#include "Point.h"
#include "QuadratureCollection.h"
#include "SubMatrix.h"
#include "config.h"

DROMON_NAMESPACE_OPEN
template <class DoFCellType, class CoefficientType, class Real = double>
class EFIEIntegrator
    : public IntegratorBase<DoFCellType, CoefficientType, Real> {
public:
  explicit EFIEIntegrator() = default;

  virtual void integrate(
      const DoFCellType &cell_test, const DoFCellType &cell_trial,
      const DoFMask &mask_test, const DoFMask &mask_trial,
      const MaterialData<Real> &material_data,
      const Excitations::Excitation<DoFCellType::space_dim, Real> &excitation,
      DenseSubMatrix<CoefficientType> *output, const unsigned int &ngl_order = 0, const unsigned int& ngl_vertex_order = 0, const unsigned int& ngl_edge_order = 0, const unsigned int& ngl_self_order = 0) override;

  virtual void integrate_HOPS(
      const DoFCellType &cell_test, const DoFCellType &cell_trial,
      const DoFMask &mask_test, const DoFMask &mask_trial,
      const MaterialData<Real> &material_data,
      const Excitations::Excitation<DoFCellType::space_dim, Real> &excitation,
      DenseSubMatrix<CoefficientType> *output, const unsigned int &ngl_order = 0, const unsigned int& ngl_vertex_order = 0, const unsigned int& ngl_edge_order = 0, const unsigned int& ngl_self_order = 0);
  
  virtual void integrate_perturb(
      const DoFCellType &cell_test, const DoFCellType &cell_trial,
      const DoFMask &mask_test, const DoFMask &mask_trial,
      const MaterialData<Real> &material_data,
      const Excitations::Excitation<DoFCellType::space_dim, Real> &excitation,
      DenseSubMatrix<CoefficientType> *output, const unsigned int &ngl_order = 0, const unsigned int& ngl_vertex_order = 0, const unsigned int& ngl_edge_order = 0, const unsigned int& ngl_self_order = 0);
  
  virtual void integrate(const DoFCellType &cell_test,
                         const DoFCellType &cell_trial,
                         const DoFMask &mask_test, const DoFMask &mask_trial,
                         DenseSubMatrix<CoefficientType> *output) override;

  // virtual void integrate_HOPS(const DoFCellType &cell_test,
  //                        const DoFCellType &cell_trial,
  //                        const DoFMask &mask_test, const DoFMask &mask_trial,
  //                        DenseSubMatrix<CoefficientType> *output);

  virtual void fill_excitation(
      const DoFCellType &cell_test, const DoFMask &mask_test,
      const Excitations::Excitation<DoFCellType::space_dim, Real> &excitation,
      const unsigned int& ngl, DenseSubVector<CoefficientType> *output) override;

  virtual void fill_excitation_HOPS(
      const DoFCellType &cell_test, const DoFMask &mask_test,
      const Excitations::Excitation<DoFCellType::space_dim, Real> &excitation,
      const unsigned int& ngl, const MaterialData<Real> &material_data, DenseSubVector<CoefficientType> *output);

    virtual void fill_excitation_perturb(
      const DoFCellType &cell_test, const DoFMask &mask_test,
      const Excitations::Excitation<DoFCellType::space_dim, Real> &excitation,
      const unsigned int& ngl, const MaterialData<Real> &material_data, DenseSubVector<CoefficientType> *output);


  virtual bool is_symmetric() const override;
  virtual const std::string get_integrator_name() const override;

private:

  //   Pointer to the FE, FECollection, or FECollectionCollector...
  //  FE *fe;
  //  std::unique_ptr<
  //      NumericalIntegration::DIRECTFN_ST<DoFCellType, CoefficientType, Real>>
  //      directfnSt;
  //
  //  std::unique_ptr<
  //      NumericalIntegration::DIRECTFN_ET<DoFCellType, CoefficientType, Real>>
  //      directfnEt;
  //
  //  std::unique_ptr<
  //      NumericalIntegration::DIRECTFN_VT<DoFCellType, CoefficientType, Real>>
  //      directfnVt;

  template <class TestComponent, class TrialComponent, class Integrand,
            class Kernel>
  void integrate_regular(
      Integrand integrand, Kernel kernel, TestComponent testComponent,
      TrialComponent trialComponent, const DoFCellType &cell_test,
      const DoFCellType &cell_trial, const DoFMask &mask_test,
      const DoFMask &mask_trial, const MaterialData<Real> &material_data,
      const Excitations::Excitation<DoFCellType::space_dim, Real> &excitation,const unsigned int& ngl,
      DenseSubMatrix<CoefficientType> *output);

  void find_regularity_type(
      const DoFCellType &cell_test, const DoFCellType &cell_trial,
      RegularityType *regularity,
      MultiIndex<2, unsigned int> *regularity_local_index) const;
  Quadrature::QuadratureCollection<Quadrature::GaussQuadrature> qcollection;

  //  //
  //  RegularityType regularity;
  //  // In the event that there is a (near-)singularity, @p regularity
  //  local_index
  //  // contains the local indices (e.g., the edge index) of the singularity
  //  for
  //  // the two cells being tested
  //  MultiIndex<2, unsigned int> regularity_local_index;

};
/**
 * As the EFIE is symmetric with respect to the trial and testing cells,
 * the integrals only need to be evaluated once...
 * @tparam DoFCellType
 * @tparam CoefficientType
 * @return
 */
template <class DoFCellType, class CoefficientType, class Real>
bool EFIEIntegrator<DoFCellType, CoefficientType, Real>::is_symmetric() const {
  return true;
}
template <class DoFCellType, class CoefficientType, class Real>
const std::string
EFIEIntegrator<DoFCellType, CoefficientType, Real>::get_integrator_name()
    const {
  return "EFIEIntegrator";
}
/**
 * Here, we integrate the three possible combinations of depending on the
 * orientation of the cells with respect to each other, namely their distance
 * which engenders the (near-)singular behavior
 * @tparam DoFCellType
 * @tparam CoefficientType
 * @param cell_test
 * @param cell_trial
 * @param mask_test
 * @param mask_trial
 * @param output
 */

template <class DoFCellType, class CoefficientType, class Real>
void EFIEIntegrator<DoFCellType, CoefficientType, Real>::integrate(
    const DoFCellType &cell_test, const DoFCellType &cell_trial,
    const DoFMask &mask_test, const DoFMask &mask_trial,
    const MaterialData<Real> &material_data,
    const Excitations::Excitation<DoFCellType::space_dim, Real> &excitation,
    DenseSubMatrix<CoefficientType> *output, const unsigned int &ngl_order, const unsigned int& ngl_vertex_order, const unsigned int& ngl_edge_order, const unsigned int& ngl_self_order) {

  // Setup the quadrature orders for the given problem...
  const unsigned int quadrature_order_regular = ngl_order == 0 ? 5 : ngl_order;
  const unsigned int sing_quadrature_order_face = ngl_order == 0 ? 5 : ngl_edge_order;
  const unsigned int sing_quadrature_order_self = ngl_order == 0 ? 5 : ngl_self_order;
  const unsigned int sing_quadrature_order_point = ngl_order == 0 ? 5 : ngl_vertex_order;


  const std::complex<Real> complexj = {0., 1.};
  const auto cell_test_geom = cell_test.get_cell_on_mesh();


  const Material<Real> &exterior_material =
      (material_data.get_material_domain(cell_test_geom->material_domain_id()))
          .get_exterior();

  auto kernel = [&excitation, &exterior_material](const Real &R) {
    return kernels::g(R, excitation.omega, exterior_material.epsr,
                      exterior_material.mur);
    // return kernels::gScaled(R, excitation.omega, exterior_material.epsr,
    //                   exterior_material.mur, 2*constants<Real>::PI/excitation.gamma.imag());
  };
  auto test_component = [](const DoFCellType &test_cell,
                           const unsigned int &test_dof_index,
                           const Point<DoFCellType::cell_dim, Real> &uv_test) {
    const DoFBase<DoFCellType::cell_dim> &dof_test =
        test_cell.get_dof(test_dof_index);

    Point<2, Real> out = {dof_test.evaluate_shape_function(uv_test),
                          dof_test.evaluate_shape_function_divergence(uv_test)};
    return out;
  };
  auto trial_component = [](const DoFCellType &trial_cell,
                                      const unsigned int &trial_dof_index,
                                      const Point<DoFCellType::cell_dim, Real> &uv_trial){
    const DoFBase<DoFCellType::cell_dim> &dof_trial =
        trial_cell.get_dof(trial_dof_index);

    Point<2, Real> out = {dof_trial.evaluate_shape_function(uv_trial),
                          dof_trial.evaluate_shape_function_divergence(uv_trial)};
    return out;
  };
  // Let's generate a lambda function for the integrand
  ;

  // Test_values and trial_values are those returned by test_component and trial_component, respectively
  auto combined_integrand = [&excitation, &exterior_material,&complexj](const Point<2, Real>& test_values,
                                                              const Point<2, Real>& trial_values,
                               const Point<DoFCellType::space_dim, Real> &unitary_test,
                               const Point<DoFCellType::space_dim, Real> &unitary_trial,
                               const CoefficientType &kernel_value, Real R)
  {

        // std::cout << "R value in combined_integrand: " << R << std::endl;

        const Real interm_comp_0 = -test_values(0)*trial_values(0)*unitary_test.dot(unitary_trial)*excitation.omega;
        const Real interm_comp_1 = test_values(1)*trial_values(1)/excitation.omega;
       // return (complexj*kernel_value)*(interm_comp_0*exterior_material.mu + interm_comp_1/exterior_material.eps);
        const auto temp_out = kernel_value*(interm_comp_0*exterior_material.mu + interm_comp_1/exterior_material.eps);
        return std::complex<Real>(-temp_out.imag(), temp_out.real());
  };

  // First, we need to find the regularity type
  RegularityType regularity;
  MultiIndex<2, unsigned int> regularity_local_index;
  this->find_regularity_type(cell_test, cell_trial, &regularity,
                             &regularity_local_index);

  // Then, based on the regularity type, we execute the possible options...
  switch (regularity) {
  case Regular:
    this->integrate_regular(combined_integrand, kernel, test_component, trial_component, cell_test, cell_trial, mask_test,
                            mask_trial, material_data, excitation, quadrature_order_regular, output);
    break;
  case EdgeSingular:
  case FaceSingular: {
    // Make the DIRECTFN object
    auto directfnEt = std::make_unique<
        NumericalIntegration::DIRECTFN_ET2<DoFCellType, CoefficientType, Real>>(
        NumericalIntegration::DIRECTFN_ET2<DoFCellType, CoefficientType, Real>(
            sing_quadrature_order_face, sing_quadrature_order_face, ceil(2.0*sing_quadrature_order_face/3.0),
            ceil(2.0*sing_quadrature_order_face/3.0)));

    // Now, actually perform integration
    (directfnEt)
        ->integrate_singular(combined_integrand, kernel, test_component, trial_component, cell_test, cell_trial,
                             mask_test, mask_trial, regularity_local_index, output);
  } break;
  case PointSingular: {
    // Make the DIRECTFN object
   // unsigned int sing_quadrature_order = cell_test.active_degree(0) * 10;
    auto directfnVt = std::make_unique<
        NumericalIntegration::DIRECTFN_VT2<DoFCellType, CoefficientType, Real>>(
        NumericalIntegration::DIRECTFN_VT2<DoFCellType, CoefficientType, Real>(
            sing_quadrature_order_point, sing_quadrature_order_point, sing_quadrature_order_point,
            sing_quadrature_order_point));

    // Now, actually perform integration
    (directfnVt)
        ->integrate_singular(combined_integrand, kernel, test_component, trial_component, cell_test, cell_trial,
                             mask_test, mask_trial, regularity_local_index, output);
  } break;
  case SelfSingular:

    {
      auto directfnSt = std::make_unique<NumericalIntegration::DIRECTFN_ST2<DoFCellType, CoefficientType, Real>>(NumericalIntegration::DIRECTFN_ST2<DoFCellType, CoefficientType, Real>(
                       sing_quadrature_order_self, ceil(2.*sing_quadrature_order_self/3.),
          ceil(2.*sing_quadrature_order_self/3.), ceil(2.*sing_quadrature_order_self/3.)));
      // Now, actually perform integration
      (directfnSt)
          ->integrate_singular(combined_integrand, kernel, test_component, trial_component, cell_test, cell_trial,
                               mask_test, mask_trial, output);
    }

    // assert(false && "Other regularity cases not currently implemented!");
    break;
  case EdgeNearSingular:
  case PointNearSingular:
    assert(false && "Other regularity cases not currently implemented!");
    break;
  default:
    assert(false && "Other regularity cases not currently implemented!");
  }
}


template <class DoFCellType, class CoefficientType, class Real>
void EFIEIntegrator<DoFCellType, CoefficientType, Real>::integrate_HOPS(
    const DoFCellType &cell_test, const DoFCellType &cell_trial,
    const DoFMask &mask_test, const DoFMask &mask_trial,
    const MaterialData<Real> &material_data,
    const Excitations::Excitation<DoFCellType::space_dim, Real> &excitation,
    DenseSubMatrix<CoefficientType> *output, const unsigned int &ngl_order, const unsigned int& ngl_vertex_order, const unsigned int& ngl_edge_order, const unsigned int& ngl_self_order) {

  // Setup the quadrature orders for the given problem...
  const unsigned int quadrature_order_regular = ngl_order == 0 ? 5 : ngl_order;
  const unsigned int sing_quadrature_order_face = ngl_order == 0 ? 5 : ngl_edge_order;
  const unsigned int sing_quadrature_order_self = ngl_order == 0 ? 5 : ngl_self_order;
  const unsigned int sing_quadrature_order_point = ngl_order == 0 ? 5 : ngl_vertex_order;


  const std::complex<Real> complexj = {0., 1.};
  const auto cell_test_geom = cell_test.get_cell_on_mesh();

  const Material<Real> &exterior_material =
      (material_data.get_material_domain(cell_test_geom->material_domain_id()))
          .get_exterior();

  auto kernel = [&excitation, &exterior_material](const Real &R) {
    return kernels::g(R, excitation.omega, exterior_material.epsr,
                      exterior_material.mur);
    // return kernels::gScaled(R, excitation.omega, exterior_material.epsr,
    //                   exterior_material.mur, 2*constants<Real>::PI/excitation.gamma.imag());
  };
  auto test_component = [](const DoFCellType &test_cell,
                           const unsigned int &test_dof_index,
                           const Point<DoFCellType::cell_dim, Real> &uv_test) {
    const DoFBase<DoFCellType::cell_dim> &dof_test =
        test_cell.get_dof(test_dof_index);

    Point<2, Real> out = {dof_test.evaluate_shape_function(uv_test),
                          dof_test.evaluate_shape_function_divergence(uv_test)};
    return out;
  };
  auto trial_component = [](const DoFCellType &trial_cell,
                                      const unsigned int &trial_dof_index,
                                      const Point<DoFCellType::cell_dim, Real> &uv_trial){
    const DoFBase<DoFCellType::cell_dim> &dof_trial =
        trial_cell.get_dof(trial_dof_index);

    Point<2, Real> out = {dof_trial.evaluate_shape_function(uv_trial),
                          dof_trial.evaluate_shape_function_divergence(uv_trial)};
    return out;
  };
  // Let's generate a lambda function for the integrand
  ;

  // Test_values and trial_values are those returned by test_component and trial_component, respectively
  auto combined_integrand_HOPS = [&excitation, &exterior_material,&complexj, this](const Point<2, Real>& test_values,
                                                              const Point<2, Real>& trial_values,
                               const Point<DoFCellType::space_dim, Real> &unitary_test,
                               const Point<DoFCellType::space_dim, Real> &unitary_trial,
                               const CoefficientType &kernel_value, Real R)
  {

        // std::cout << "sqrt(epsmu): " << sqrt(exterior_material.eps*exterior_material.mu) << std::endl;

      // // ////omega version------------------
        // const std::complex<double> d1 (1.0/excitation.omega, -R*sqrt(8.8541878128e-12 * 12.566370614359e-7));
        // const std::complex<double> d2 (-1.0/excitation.omega, -R*sqrt(8.8541878128e-12 * 12.566370614359e-7));

      ////// epsilon-----------------------------------
      // const std::complex<double> d1 (0.0, -exterior_material.mu*excitation.omega*R/2.0/sqrt(exterior_material.eps*exterior_material.mu));
      // const std::complex<double> d2 (-1./exterior_material.eps, -excitation.omega*R*exterior_material.mu/sqrt(exterior_material.eps*exterior_material.mu)/2.0);

      const std::complex<double> d1 = -complexj*exterior_material.mu*excitation.omega*R/2.0/sqrt(exterior_material.eps*exterior_material.mu);
      const std::complex<double> d2 = -1./exterior_material.eps - complexj*excitation.omega*R*exterior_material.mu/sqrt(exterior_material.eps*exterior_material.mu)/2.0;

        const std::complex<double> interm_comp_0 = -complexj*test_values(0)*trial_values(0)*unitary_test.dot(unitary_trial)*excitation.omega*exterior_material.mu;
        const std::complex<double> interm_comp_1 = complexj*test_values(1)*trial_values(1)/excitation.omega/exterior_material.eps;
       // return (complexj*kernel_value)*(interm_comp_0*exterior_material.mu + interm_comp_1/exterior_material.eps);
        const auto temp_out = kernel_value*(interm_comp_0*d1 + interm_comp_1*d2);
        // return std::complex<Real>(-temp_out.imag(), temp_out.real());
        return temp_out;


      // // // ////omega version, pulled omega from comps and into d1/d2 grad terms-------------------
      //   const std::complex<double> interm_comp_0 = -complexj*exterior_material.mu*test_values(0)*trial_values(0)*unitary_test.dot(unitary_trial);
      //   const std::complex<double> interm_comp_1 = complexj/exterior_material.eps*test_values(1)*trial_values(1);

      //   // const auto d1 = 1/excitation.omega + kernels::dgdw(R, excitation.omega, exterior_material.epsr,
      //   //               exterior_material.mur);
        
      //   const auto d1 = 1.0 - complexj*excitation.omega*R*sqrt(8.8541878128e-12 * 12.566370614359e-7);
      //   // const std::complex<Real> d1(1/excitation.omega, -sqrt(8.8541878128e-12 * 12.566370614359e-7)*R);


      //   const auto d2 = (-1.0 - complexj*excitation.omega*R*sqrt(8.8541878128e-12 * 12.566370614359e-7))/excitation.omega/excitation.omega;
      //  // return (complexj*kernel_value)*(interm_comp_0*exterior_material.mu + interm_comp_1/exterior_material.eps);
      //   // const auto temp_out = kernel_value*(interm_comp_0*exterior_material.mu * (1/excitation.omega -complexj*R*sqrt(8.8541878128e-12 * 12.566370614359e-7))
      //   //  + interm_comp_1/exterior_material.eps * (-1/excitation.omega -complexj*R*sqrt(8.8541878128e-12 * 12.566370614359e-7)));
      //   //  + interm_comp_1/exterior_material.eps * std::complex<Real>(-1/excitation.omega, -R*sqrt(8.8541878128e-12 * 12.566370614359e-7)));
      //    const auto temp_out = kernel_value*(interm_comp_0*d1 + interm_comp_1*d2);

      /////new version---------------------
      //   const Real interm_comp_0 = -test_values(0)*trial_values(0)*unitary_test.dot(unitary_trial);
      //   const Real interm_comp_1 = test_values(1)*trial_values(1);
      //  // return (complexj*kernel_value)*(interm_comp_0*exterior_material.mu + interm_comp_1/exterior_material.eps);
      //   const auto temp_out = kernel_value*(interm_comp_0*exterior_material.mu * (1.0 - excitation.omega*complexj*R*sqrt(8.8541878128e-12 * 12.566370614359e-7))
      //    + interm_comp_1/exterior_material.eps * (-1.0 - excitation.omega*complexj*R*sqrt(8.8541878128e-12 * 12.566370614359e-7))/excitation.omega/excitation.omega);
      //   //  + interm_comp_1/exterior_material.eps * std::complex<Real>(-1/excitation.omega, -R*sqrt(8.8541878128e-12 * 12.566370614359e-7)));
      //   //  const auto temp_out = kernel_value*(interm_comp_0*exterior_material.mu*R + interm_comp_1/exterior_material.eps*R);

        // return std::complex<Real>(-temp_out.imag(), temp_out.real());
        

  };


    auto combined_integrand = [&excitation, &exterior_material,&complexj](const Point<2, Real>& test_values,
                                                              const Point<2, Real>& trial_values,
                               const Point<DoFCellType::space_dim, Real> &unitary_test,
                               const Point<DoFCellType::space_dim, Real> &unitary_trial,
                               const CoefficientType &kernel_value, Real R)
  {

        // std::cout << "R value in combined_integrand: " << R << std::endl;

        const Real interm_comp_0 = -test_values(0)*trial_values(0)*unitary_test.dot(unitary_trial)*excitation.omega;
        const Real interm_comp_1 = test_values(1)*trial_values(1)/excitation.omega;
       // return (complexj*kernel_value)*(interm_comp_0*exterior_material.mu + interm_comp_1/exterior_material.eps);
        const auto temp_out = kernel_value*(interm_comp_0*exterior_material.mu + interm_comp_1/exterior_material.eps);
        return std::complex<Real>(-temp_out.imag(), temp_out.real());
  };

  // First, we need to find the regularity type
  RegularityType regularity;
  MultiIndex<2, unsigned int> regularity_local_index;
  this->find_regularity_type(cell_test, cell_trial, &regularity,
                             &regularity_local_index);

  // Then, based on the regularity type, we execute the possible options...
  switch (regularity) {
  case Regular:
    this->integrate_regular(combined_integrand_HOPS, kernel, test_component, trial_component, cell_test, cell_trial, mask_test,
                            mask_trial, material_data, excitation, quadrature_order_regular, output);
    break;
  case EdgeSingular:
  case FaceSingular: {
    // Make the DIRECTFN object
    auto directfnEt = std::make_unique<
        NumericalIntegration::DIRECTFN_ET2<DoFCellType, CoefficientType, Real>>(
        NumericalIntegration::DIRECTFN_ET2<DoFCellType, CoefficientType, Real>(
            sing_quadrature_order_face, sing_quadrature_order_face, ceil(2.0*sing_quadrature_order_face/3.0),
            ceil(2.0*sing_quadrature_order_face/3.0)));

    // Now, actually perform integration
    (directfnEt)
        ->integrate_singular(combined_integrand_HOPS, kernel, test_component, trial_component, cell_test, cell_trial,
                             mask_test, mask_trial, regularity_local_index, output);
  } break;
  case PointSingular: {
    // Make the DIRECTFN object
   // unsigned int sing_quadrature_order = cell_test.active_degree(0) * 10;
    auto directfnVt = std::make_unique<
        NumericalIntegration::DIRECTFN_VT2<DoFCellType, CoefficientType, Real>>(
        NumericalIntegration::DIRECTFN_VT2<DoFCellType, CoefficientType, Real>(
            sing_quadrature_order_point, sing_quadrature_order_point, sing_quadrature_order_point,
            sing_quadrature_order_point));

    // Now, actually perform integration
    (directfnVt)
        ->integrate_singular(combined_integrand_HOPS, kernel, test_component, trial_component, cell_test, cell_trial,
                             mask_test, mask_trial, regularity_local_index, output);
  } break;
  case SelfSingular:

    {

      auto directfnSt = std::make_unique<NumericalIntegration::DIRECTFN_ST2<DoFCellType, CoefficientType, Real>>(NumericalIntegration::DIRECTFN_ST2<DoFCellType, CoefficientType, Real>(
                       sing_quadrature_order_self, ceil(2.*sing_quadrature_order_self/3.),
          ceil(2.*sing_quadrature_order_self/3.), ceil(2.*sing_quadrature_order_self/3.)));
      // Now, actually perform integration
      (directfnSt)
          ->integrate_singular(combined_integrand_HOPS, kernel, test_component, trial_component, cell_test, cell_trial,
                               mask_test, mask_trial, output);
    }

    // assert(false && "Other regularity cases not currently implemented!");
    break;
  case EdgeNearSingular:
  case PointNearSingular:
    assert(false && "Other regularity cases not currently implemented!");
    break;
  default:
    assert(false && "Other regularity cases not currently implemented!");
  }
}

template <class DoFCellType, class CoefficientType, class Real>
void EFIEIntegrator<DoFCellType, CoefficientType, Real>::integrate_perturb(
    const DoFCellType &cell_test, const DoFCellType &cell_trial,
    const DoFMask &mask_test, const DoFMask &mask_trial,
    const MaterialData<Real> &material_data,
    const Excitations::Excitation<DoFCellType::space_dim, Real> &excitation,
    DenseSubMatrix<CoefficientType> *output, const unsigned int &ngl_order, const unsigned int& ngl_vertex_order, const unsigned int& ngl_edge_order, const unsigned int& ngl_self_order) {

  // Setup the quadrature orders for the given problem...
  const unsigned int quadrature_order_regular = ngl_order == 0 ? 5 : ngl_order;
  const unsigned int sing_quadrature_order_face = ngl_order == 0 ? 5 : ngl_edge_order;
  const unsigned int sing_quadrature_order_self = ngl_order == 0 ? 5 : ngl_self_order;
  const unsigned int sing_quadrature_order_point = ngl_order == 0 ? 5 : ngl_vertex_order;


  const std::complex<Real> complexj = {0., 1.};
  const auto cell_test_geom = cell_test.get_cell_on_mesh();

  const Material<Real> &exterior_material =
      (material_data.get_material_domain(cell_test_geom->material_domain_id()))
          .get_exterior();

  auto kernel = [&excitation, &exterior_material](const Real &R) {

//----------------------------------------------------------------------------------------------------
    //shouldn't need to perturb here
    //original
    return kernels::g(R, excitation.omega, exterior_material.epsr,
                      exterior_material.mur);

    //

      //// perturbed material 
    // double epsr_perturb = exterior_material.epsr*1.01;
    // double eps_perturb = exterior_material.eps*1.01;
    // return kernels::g(R, excitation.omega, epsr_perturb,
    //                   exterior_material.mur);
//----------------------------------------------------------------------------------------------------
  };
  auto test_component = [](const DoFCellType &test_cell,
                           const unsigned int &test_dof_index,
                           const Point<DoFCellType::cell_dim, Real> &uv_test) {
    const DoFBase<DoFCellType::cell_dim> &dof_test =
        test_cell.get_dof(test_dof_index);

    Point<2, Real> out = {dof_test.evaluate_shape_function(uv_test),
                          dof_test.evaluate_shape_function_divergence(uv_test)};
    return out;
  };
  auto trial_component = [](const DoFCellType &trial_cell,
                                      const unsigned int &trial_dof_index,
                                      const Point<DoFCellType::cell_dim, Real> &uv_trial){
    const DoFBase<DoFCellType::cell_dim> &dof_trial =
        trial_cell.get_dof(trial_dof_index);

    Point<2, Real> out = {dof_trial.evaluate_shape_function(uv_trial),
                          dof_trial.evaluate_shape_function_divergence(uv_trial)};
    return out;
  };
  // Let's generate a lambda function for the integrand
  ;

    auto combined_integrand_perturb = [&excitation, &exterior_material,&complexj](const Point<2, Real>& test_values,
                                                              const Point<2, Real>& trial_values,
                               const Point<DoFCellType::space_dim, Real> &unitary_test,
                               const Point<DoFCellType::space_dim, Real> &unitary_trial,
                               const CoefficientType &kernel_value, Real R)
  {

//-------------------------------------------------------------------------------------------------------------------
          //// perturbed material 
        // double epsr_perturb = exterior_material.epsr*1.01;
        // double eps_perturb = exterior_material.eps*1.01;
        // std::complex<double> g = kernels::g(R*1.01, excitation.omega, epsr_perturb, exterior_material.mur);


        // perturbed radius
        std::complex<double> g = kernels::g(R*1.01, excitation.omega, exterior_material.epsr, exterior_material.mur);



        const Real interm_comp_0 = -test_values(0)*trial_values(0)*unitary_test.dot(unitary_trial)*excitation.omega;
        const Real interm_comp_1 = test_values(1)*trial_values(1)/excitation.omega;

        //original
        // const auto temp_out = kernel_value*(interm_comp_0*exterior_material.mu + interm_comp_1/exterior_material.eps);
        //perturbed material
        // const auto temp_out = kernel_value*(interm_comp_0*exterior_material.mu + interm_comp_1/eps_perturb);
        //perturbed radius
        const auto temp_out = g*(interm_comp_0*exterior_material.mu + interm_comp_1/exterior_material.eps);

//--------------------------------------------------------------------------------------------------------------------

        return std::complex<Real>(-temp_out.imag(), temp_out.real());
  };

  // First, we need to find the regularity type
  RegularityType regularity;
  MultiIndex<2, unsigned int> regularity_local_index;
  this->find_regularity_type(cell_test, cell_trial, &regularity,
                             &regularity_local_index);

  // Then, based on the regularity type, we execute the possible options...
  switch (regularity) {
  case Regular:
    this->integrate_regular(combined_integrand_perturb, kernel, test_component, trial_component, cell_test, cell_trial, mask_test,
                            mask_trial, material_data, excitation, quadrature_order_regular, output);
    break;
  case EdgeSingular:
  case FaceSingular: {
    // Make the DIRECTFN object
    auto directfnEt = std::make_unique<
        NumericalIntegration::DIRECTFN_ET2<DoFCellType, CoefficientType, Real>>(
        NumericalIntegration::DIRECTFN_ET2<DoFCellType, CoefficientType, Real>(
            sing_quadrature_order_face, sing_quadrature_order_face, ceil(2.0*sing_quadrature_order_face/3.0),
            ceil(2.0*sing_quadrature_order_face/3.0)));

    // Now, actually perform integration
    (directfnEt)
        ->integrate_singular(combined_integrand_perturb, kernel, test_component, trial_component, cell_test, cell_trial,
                             mask_test, mask_trial, regularity_local_index, output);
  } break;
  case PointSingular: {
    // Make the DIRECTFN object
   // unsigned int sing_quadrature_order = cell_test.active_degree(0) * 10;
    auto directfnVt = std::make_unique<
        NumericalIntegration::DIRECTFN_VT2<DoFCellType, CoefficientType, Real>>(
        NumericalIntegration::DIRECTFN_VT2<DoFCellType, CoefficientType, Real>(
            sing_quadrature_order_point, sing_quadrature_order_point, sing_quadrature_order_point,
            sing_quadrature_order_point));

    // Now, actually perform integration
    (directfnVt)
        ->integrate_singular(combined_integrand_perturb, kernel, test_component, trial_component, cell_test, cell_trial,
                             mask_test, mask_trial, regularity_local_index, output);
  } break;
  case SelfSingular:

    {

      auto directfnSt = std::make_unique<NumericalIntegration::DIRECTFN_ST2<DoFCellType, CoefficientType, Real>>(NumericalIntegration::DIRECTFN_ST2<DoFCellType, CoefficientType, Real>(
                       sing_quadrature_order_self, ceil(2.*sing_quadrature_order_self/3.),
          ceil(2.*sing_quadrature_order_self/3.), ceil(2.*sing_quadrature_order_self/3.)));
      // Now, actually perform integration
      (directfnSt)
          ->integrate_singular(combined_integrand_perturb, kernel, test_component, trial_component, cell_test, cell_trial,
                               mask_test, mask_trial, output);
    }

    // assert(false && "Other regularity cases not currently implemented!");
    break;
  case EdgeNearSingular:
  case PointNearSingular:
    assert(false && "Other regularity cases not currently implemented!");
    break;
  default:
    assert(false && "Other regularity cases not currently implemented!");
  }
}

/**
 *
 * @tparam DoFCellType
 * @tparam CoefficientType
 * @param cell_test
 * @param cell_trial
 * @param regularity The type of regularity (e.g., EdgeSingularity,
 * NodeSingularity)
 * @param regularity_local_index In the event that there is a
 * (near-)singularity, @p regularity_local_index contains the local indices
 * (e.g., the edge index) of the singularity for the two cells being tested
 */
template <class DoFCellType, class CoefficientType, class Real>
void EFIEIntegrator<DoFCellType, CoefficientType, Real>::find_regularity_type(
    const DoFCellType &cell_test, const DoFCellType &cell_trial,
    RegularityType *regularity,
    MultiIndex<2, unsigned int> *regularity_local_index) const {

  auto mesh_cell_test = cell_test.get_cell_on_mesh();
  auto mesh_cell_trial = cell_trial.get_cell_on_mesh();
  AdjacencyType adjacencyType;
  mesh_cell_test->check_if_adjacent(mesh_cell_trial, &adjacencyType,
                                    regularity_local_index);

  if (adjacencyType == Disjoint)
    *regularity = Regular;
  else if (adjacencyType == Self)
    *regularity = SelfSingular;
  else if (adjacencyType == Edge)
    *regularity = EdgeSingular;
  else if (adjacencyType == Vertex)
    *regularity = PointSingular;
  else if (adjacencyType == Face)
    *regularity = FaceSingular;
}
template <class DoFCellType, class CoefficientType, class Real>
template <class TestComponent, class TrialComponent, class Integrand,
          class Kernel>
void EFIEIntegrator<DoFCellType, CoefficientType, Real>::integrate_regular(
    Integrand integrand, Kernel kernel, TestComponent testComponent,
    TrialComponent trialComponent, const DoFCellType &cell_test,
    const DoFCellType &cell_trial, const DoFMask &mask_test,
    const DoFMask &mask_trial, const MaterialData<Real> &material_data,
    const Excitations::Excitation<DoFCellType::space_dim, Real> &excitation, const unsigned int& ngl,
    DenseSubMatrix<CoefficientType> *output) {
    using ComponentsReturnType = Point<2, Real>;

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
  // In this case, we simply perform standard dim-D numerical integration

  // For the first component we have the integral
  // I_EE_1 = -jwu int_S_{test} int_S_{trial} J_{test}\cdot J_{trial} g
  // dS_{trial} dS_{test} For the second part we have the integral I_EE_2 =
  // j/w/e int_S{test} int_S{trial} (div(J_{test})*(div(J_{trial})) g dS_{trial}
  // dS_{test}
  std::vector<double> *xgl, *wgl;
  //const unsigned int ngl = cell_test.active_degree(0) * 10;
  qcollection.push_back(ngl);
  qcollection.get_weights_and_points(ngl, wgl, xgl);
  const std::complex<Real> complexj = {0., 1.};

  const unsigned dim = cell_test.cell_dim;
  switch (dim) {
  case 2: {
    auto cell_test_geom = cell_test.get_cell_on_mesh();
    auto cell_trial_geom = cell_trial.get_cell_on_mesh();
    const Material<Real> exterior_material =
        (material_data.get_material_domain(
             cell_test_geom->material_domain_id()))
            .get_exterior();

    for (unsigned int i_test = 0; i_test < wgl->size(); ++i_test) {
      const Real jac0 = wgl->operator[](i_test);
      for (unsigned int j_test = 0; j_test < wgl->size(); ++j_test) {
        const Real jac1 = jac0 * wgl->operator[](j_test);
        const Point<dim, Real> uv_test(xgl->operator[](i_test),
                                       xgl->operator[](j_test));
        // Now we must loop through all the degrees of freedom
        // For those degrees of freedom that are masked, we skip them
        const auto unitary_test_u =
            cell_test_geom->unitary_vector(uv_test, u_dir);
        const auto unitary_test_v =
            cell_test_geom->unitary_vector(uv_test, v_dir);

        auto r_test = cell_test_geom->r(uv_test);

        // Now that we've computed uv_test, we can precompute the values of the shape functions
        ComponentsReturnType dof_test_memoize[mask_test.size()];
        for (unsigned int dof_test_index = 0;
             dof_test_index < mask_test.size(); ++dof_test_index) {
          if (!test_has_integrals_to_perform[dof_test_index])
            continue;

          dof_test_memoize[dof_test_index] = testComponent(cell_test, dof_test_index, uv_test);
        }

        for (unsigned int i_trial = 0; i_trial < wgl->size(); ++i_trial) {
          const Real jac2 = jac1 * wgl->operator[](i_trial);
          for (unsigned int j_trial = 0; j_trial < wgl->size(); ++j_trial) {
            const Real jacobian = jac2 * wgl->operator[](j_trial);
            const Point<dim, Real> uv_trial(xgl->operator[](i_trial),
                                            xgl->operator[](j_trial));
            const auto unitary_trial_u =
                cell_trial_geom->unitary_vector(uv_trial, u_dir);
            const auto unitary_trial_v =
                cell_trial_geom->unitary_vector(uv_trial, v_dir);
            auto r_trial = cell_trial_geom->r(uv_trial);

            Real R = (r_test - r_trial).norm();
            // std::cout << "R in method: " << R << std::endl;
            // std::cout << "r_test: " << r_test[0] << ", " << r_test[1] << ", " << r_test[2] << std::endl;
            // std::cout << "r_trial: " << r_trial[0] << ", " << r_trial[1] << ", " << r_trial[2] << std::endl;

            const auto kernel_value = kernel(R);

            // std::cout << "R after kernel: " << R << std::endl;

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

              // Turn the index into the DoF object
              const DoFBase<dim> &dof_test = cell_test.get_dof(dof_test_index);
              const auto& test_component_value = dof_test_memoize[dof_test_index];

              for (unsigned int dof_trial_index = 0;
                   dof_trial_index < mask_trial.size(); ++dof_trial_index) {
                if (is_filled(dof_test_index, dof_trial_index) == true)
                  continue;

                // Turn the index into the DoF object
                const DoFBase<dim> &dof_trial =
                    cell_trial.get_dof(dof_trial_index);

                const auto &trial_component_value =
                    dof_trial_memoize[dof_trial_index];

                // if(dof_trial_index == 0 && dof_test_index == 0)
                // stdLL

                const auto value_to_add = jacobian *
                                          integrand(test_component_value, trial_component_value,
                                                    dof_test.direction == u_dir ? unitary_test_u
                                                                                : unitary_test_v,
                                                    dof_trial.direction == u_dir ? unitary_trial_u
                                                                                 : unitary_trial_v,
                                                    kernel_value, R);
                temp_matrix.accumulate(dof_test_index, dof_trial_index,value_to_add);

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
      const DoFBase<dim> &dof_test = cell_test.get_dof(i);
      for (unsigned int j = 0; j < mask_trial.size(); ++j) {
        if (is_filled(i,j))
          continue;
        const DoFBase<dim> &dof_trial = cell_trial.get_dof(j);

        output->at(dof_test.global_index, dof_trial.global_index) +=
            temp_matrix(i, j);
        output->get_is_filled_at(dof_test.global_index, dof_trial.global_index) = true;
      }
    }
  } break;
  default:
    assert(false && "Not implemented!");
  }
}
// template <class DoFCellType, class CoefficientType,  class Real>
// EFIEIntegrator<DoFCellType, CoefficientType,  Real>::EFIEIntegrator(FE *fe) {
//   this->fe = fe;
//   // Now we have to go through and decide based on what type of FE we have...
// }
template <class DoFCellType, class CoefficientType, class Real>
void EFIEIntegrator<DoFCellType, CoefficientType, Real>::integrate(
    const DoFCellType &cell_test, const DoFCellType &cell_trial,
    const DoFMask &mask_test, const DoFMask &mask_trial,
    DenseSubMatrix<CoefficientType> *output) {
  assert(false && "Not implemented for this integrator!");
}
template <class DoFCellType, class CoefficientType, class Real>
void EFIEIntegrator<DoFCellType, CoefficientType, Real>::fill_excitation(
    const DoFCellType &cell_test, const DoFMask &mask_test,
    const Excitations::Excitation<DoFCellType::space_dim, Real> &excitation,const unsigned int& ngl,
    DenseSubVector<CoefficientType> *output) {
  std::vector<double> *xgl, *wgl;
  qcollection.push_back(ngl);
  qcollection.get_weights_and_points(ngl, wgl, xgl);
  const std::complex<Real> complexj = {0., 1.};

  const unsigned dim = cell_test.cell_dim;
  switch (dim) {
  case 2: {
    auto cell_test_geom = cell_test.get_cell_on_mesh();
    for (unsigned int i_test = 0; i_test < wgl->size(); ++i_test) {
      for (unsigned int j_test = 0; j_test < wgl->size(); ++j_test) {
        const Point<dim, Real> uv_test(xgl->at(i_test), xgl->at(j_test));
        // Now we must loop through all the degrees of freedom
        // For those degrees of freedom that are masked, we skip them
        auto unitary_test_u = cell_test_geom->unitary_vector(uv_test, u_dir);
        auto unitary_test_v = cell_test_geom->unitary_vector(uv_test, v_dir);

        const auto r_test = cell_test_geom->r(uv_test);

        const auto excitation_dot_au =
            excitation.evaluate_excitation_in_direction(r_test, unitary_test_u);
        const auto excitation_dot_av =
            excitation.evaluate_excitation_in_direction(r_test, unitary_test_v);

        const Real jacobian = wgl->at(i_test) * wgl->at(j_test);
        for (unsigned int dof_test_index = 0; dof_test_index < mask_test.size();
             ++dof_test_index) {
          if (!mask_test[dof_test_index])
            continue;
          // Turn the index into the DoF object
          const DoFBase<dim> &dof_test = cell_test.get_dof(dof_test_index);
          // Access the pointer to the FiniteElement from the DoF
          //          auto fe_pointer_test =
          //              cell_test.get_active_fe_pointer(dof_test_index);
          const auto dof_value = dof_test.evaluate_shape_function(uv_test);

          // The negative sign is to account for the excitation being on the RHS
          // of the governing equation...
          output->at(dof_test.global_index) -=
              jacobian * dof_value *
              (dof_test.direction == u_dir ? excitation_dot_au
                                           : excitation_dot_av);
        }
      }
    }

  } break;
  default:
    assert(false && "Not implemented!");
  }
}

template <class DoFCellType, class CoefficientType, class Real>
void EFIEIntegrator<DoFCellType, CoefficientType, Real>::fill_excitation_HOPS(
    const DoFCellType &cell_test, const DoFMask &mask_test,
    const Excitations::Excitation<DoFCellType::space_dim, Real> &excitation,const unsigned int& ngl,
    const MaterialData<Real> &material_data, DenseSubVector<CoefficientType> *output) {
  std::vector<double> *xgl, *wgl;
  qcollection.push_back(ngl);
  qcollection.get_weights_and_points(ngl, wgl, xgl);
  const std::complex<Real> complexj = {0., 1.};
  const unsigned dim = cell_test.cell_dim;
  switch (dim) {
  case 2: {
    auto cell_test_geom = cell_test.get_cell_on_mesh();
    for (unsigned int i_test = 0; i_test < wgl->size(); ++i_test) {
      for (unsigned int j_test = 0; j_test < wgl->size(); ++j_test) {
        const Point<dim, Real> uv_test(xgl->at(i_test), xgl->at(j_test));
        // Now we must loop through all the degrees of freedom
        // For those degrees of freedom that are masked, we skip them
        auto unitary_test_u = cell_test_geom->unitary_vector(uv_test, u_dir);
        auto unitary_test_v = cell_test_geom->unitary_vector(uv_test, v_dir);

        const Real wavelength = 2*constants<Real>::PI/excitation.gamma.imag();
        const Material<Real> &exterior_material =
      (material_data.get_material_domain(cell_test_geom->material_domain_id()))
          .get_exterior();


        auto r_test = cell_test_geom->r(uv_test);
        const auto r_test_dot_n = r_test[0] * -1.0;
        const auto excitation_dot_au =
            excitation.evaluate_excitation_in_direction(r_test, unitary_test_u);
        const auto excitation_dot_av =
            excitation.evaluate_excitation_in_direction(r_test, unitary_test_v);

        // std::cout << "full excitation: " << excitationTest[0] << " " << excitationTest[1] << " " << excitationTest[2] << std::endl;
        // std::cout << "  dotExcitation: " << excitation_dot_au << std::endl;

        const Real jacobian = wgl->at(i_test) * wgl->at(j_test);
        for (unsigned int dof_test_index = 0; dof_test_index < mask_test.size();
             ++dof_test_index) {
          if (!mask_test[dof_test_index])
            continue;
          // Turn the index into the DoF object
          const DoFBase<dim> &dof_test = cell_test.get_dof(dof_test_index);
          // Access the pointer to the FiniteElement from the DoF
          //          auto fe_pointer_test =
          //              cell_test.get_active_fe_pointer(dof_test_index);
          const auto dof_value = dof_test.evaluate_shape_function(uv_test);

          // The negative sign is to account for the excitation being on the RHS
          // of the governing equation...

          output->at(dof_test.global_index) -=
              //omega perturbation
              jacobian * dof_value  * 

              //omega perturbation
              // (-complexj * r_test_dot_n * sqrt(8.8541878128e-12 * 12.566370614359e-7)) *
              //freq perturbation
              // (-complexj * constants<Real>::PI*2.0 * r_test_dot_n * sqrt(8.8541878128e-12 * 12.566370614359e-7)) *
              //eps perturbation
              -complexj*exterior_material.mu*excitation.omega*r_test_dot_n/2.0/sqrt(exterior_material.eps*exterior_material.mu) *

              (dof_test.direction == u_dir ? excitation_dot_au
                                           : excitation_dot_av);
        }
      }
    }

  } break;
  default:
    assert(false && "Not implemented!");
  }
}

template <class DoFCellType, class CoefficientType, class Real>
void EFIEIntegrator<DoFCellType, CoefficientType, Real>::fill_excitation_perturb(
    const DoFCellType &cell_test, const DoFMask &mask_test,
    const Excitations::Excitation<DoFCellType::space_dim, Real> &excitation,const unsigned int& ngl,
    const MaterialData<Real> &material_data, DenseSubVector<CoefficientType> *output) {
  std::vector<double> *xgl, *wgl;
  qcollection.push_back(ngl);
  qcollection.get_weights_and_points(ngl, wgl, xgl);
  const std::complex<Real> complexj = {0., 1.};

  const unsigned dim = cell_test.cell_dim;
  switch (dim) {
  case 2: {
    auto cell_test_geom = cell_test.get_cell_on_mesh();
    for (unsigned int i_test = 0; i_test < wgl->size(); ++i_test) {
      for (unsigned int j_test = 0; j_test < wgl->size(); ++j_test) {
        const Point<dim, Real> uv_test(xgl->at(i_test), xgl->at(j_test));
        // Now we must loop through all the degrees of freedom
        // For those degrees of freedom that are masked, we skip them
        auto unitary_test_u = cell_test_geom->unitary_vector(uv_test, u_dir);
        auto unitary_test_v = cell_test_geom->unitary_vector(uv_test, v_dir);

        auto r_test = cell_test_geom->r(uv_test);
////------------------------------------------------------------------------------------------
        //perturbed radius
        r_test *= 1.01;
        const auto excitation_dot_au =
            excitation.evaluate_excitation_in_direction_perturb(r_test, unitary_test_u, 1.0);
        const auto excitation_dot_av =
            excitation.evaluate_excitation_in_direction_perturb(r_test, unitary_test_v, 1.0);

        //perturbed material
        // const auto excitation_dot_au =
        //     excitation.evaluate_excitation_in_direction_perturb(r_test, unitary_test_u, 1.01);
        // const auto excitation_dot_av =
        //     excitation.evaluate_excitation_in_direction_perturb(r_test, unitary_test_v, 1.01);

//----------------------------------------------------------------------------------------------

        const Real jacobian = wgl->at(i_test) * wgl->at(j_test);
        for (unsigned int dof_test_index = 0; dof_test_index < mask_test.size();
             ++dof_test_index) {
          if (!mask_test[dof_test_index])
            continue;
          // Turn the index into the DoF object
          const DoFBase<dim> &dof_test = cell_test.get_dof(dof_test_index);
          // Access the pointer to the FiniteElement from the DoF
          //          auto fe_pointer_test =
          //              cell_test.get_active_fe_pointer(dof_test_index);
          const auto dof_value = dof_test.evaluate_shape_function(uv_test);

          // The negative sign is to account for the excitation being on the RHS
          // of the governing equation...
          output->at(dof_test.global_index) -=
              jacobian * dof_value *
              (dof_test.direction == u_dir ? excitation_dot_au
                                           : excitation_dot_av);
        }
      }
    }

  } break;
  default:
    assert(false && "Not implemented!");
  }
}

// template <class DoFCellType, class CoefficientType, class Real>
// EFIEIntegrator<DoFCellType, CoefficientType, Real>::EFIEIntegrator()
//{
//
//   unsigned int sing_quadrature_order = cell_test.active_degree(0) * 10;
//   this->directfnEt = std::make_unique<
//       NumericalIntegration::DIRECTFN_ET<DoFCellType, CoefficientType, Real>>(
//       NumericalIntegration::DIRECTFN_ET<DoFCellType, CoefficientType, Real>(
//           sing_quadrature_order, sing_quadrature_order,
//           sing_quadrature_order, sing_quadrature_order));
//
//   unsigned int sing_quadrature_order = cell_test.active_degree(0) * 10;
//   this->directfnVt = std::make_unique<
//       NumericalIntegration::DIRECTFN_VT<DoFCellType, CoefficientType, Real>>(
//       NumericalIntegration::DIRECTFN_VT<DoFCellType, CoefficientType, Real>(
//           sing_quadrature_order, sing_quadrature_order,
//           sing_quadrature_order, sing_quadrature_order));
//
//     unsigned int sing_quadrature_order = cell_test.active_degree(0) * 10;
//     this->directfnSt = std::make_unique<NumericalIntegration::DIRECTFN_ST<
//         DoFCellType, CoefficientType, Real>>(
//         NumericalIntegration::DIRECTFN_ST<DoFCellType, CoefficientType,
//         Real>(
//             sing_quadrature_order, sing_quadrature_order,
//             sing_quadrature_order, sing_quadrature_order));
// }

DROMON_NAMESPACE_CLOSE

#endif // DROMON_EFIEINTEGRATOR_H