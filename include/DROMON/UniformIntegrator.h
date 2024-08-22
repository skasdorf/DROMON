//
// Created by Jake J. Harmon (jake.harmon@ieee.org) on 4/3/22.
//

#ifndef DROMON_UNIFORMINTEGRATOR_H
#define DROMON_UNIFORMINTEGRATOR_H

#include "DoFMask.h"
#include "IntegratorBase.h"
#include "QuadratureCollection.h"
#include "SubMatrix.h"
#include "config.h"

DROMON_NAMESPACE_OPEN

template <class DoFCellType, class CoefficientType, class Real = double>
class UniformIntegrator
    : public IntegratorBase<DoFCellType, CoefficientType, Real> {
public:
  UniformIntegrator() = default;
  // virtual void attach_quadrature_collection() override;
  virtual bool is_symmetric() const override;
  virtual const std::string get_integrator_name() const override;
  virtual void integrate(const DoFCellType &cell_test,
                         const DoFCellType &cell_trial,
                         const DoFMask &mask_test, const DoFMask &mask_trial,
                         DenseSubMatrix<CoefficientType> *output) override;
  virtual void integrate(const DoFCellType &cell_test,
                         const DoFCellType &cell_trial,
                         const DoFMask &mask_test, const DoFMask &mask_trial,
                         const MaterialData<Real> &material_data,
                         const Excitations::Excitation<DoFCellType::space_dim, Real> &excitation,
                         DenseSubMatrix<CoefficientType> *output, const unsigned int& ngl_order = 0) override;
  virtual void fill_excitation(const DoFCellType &cell_test, const DoFMask &mask_test,  const Excitations::Excitation<DoFCellType::space_dim, Real> &excitation,const unsigned int& ngl, DenseSubVector<CoefficientType> *output) override;

private:
  void integrate_near_singular(
      const DoFCellType &cell_test, const DoFCellType &cell_trial,
      const DoFMask &mask_test, const DoFMask &mask_trial,
      const RegularityType &regularity,
      const MultiIndex<2, unsigned int> &regularity_local_index,
      DenseSubMatrix<CoefficientType> *output);
  void
  integrate_singular(const DoFCellType &cell_test,
                     const DoFCellType &cell_trial, const DoFMask &mask_test,
                     const DoFMask &mask_trial,
                     const RegularityType &regularity,
                     const MultiIndex<2, unsigned int> &regularity_local_index,
                     DenseSubMatrix<CoefficientType> *output);
  void
  integrate_regular(const DoFCellType &cell_test, const DoFCellType &cell_trial,
                    const DoFMask &mask_test, const DoFMask &mask_trial,
                    DenseSubMatrix<CoefficientType> *output);

private:
  Quadrature::QuadratureCollection<Quadrature::GaussQuadrature> qcollection;
};

template <class DoFCellType, class CoefficientType, class Real>
bool UniformIntegrator<DoFCellType, CoefficientType, Real>::is_symmetric()
    const {
  return false;
}
template <class DoFCellType, class CoefficientType, class Real>
const std::string
UniformIntegrator<DoFCellType, CoefficientType, Real>::get_integrator_name()
    const {
  return "UniformIntegrator";
}
template <class DoFCellType, class CoefficientType, class Real>
void UniformIntegrator<DoFCellType, CoefficientType, Real>::integrate_regular(
    const DoFCellType &cell_test, const DoFCellType &cell_trial,
    const DoFMask &mask_test, const DoFMask &mask_trial,
    DenseSubMatrix<CoefficientType> *output) {
  // For this integrator, we simply ignore the second cell and compute the
  // area of the first cell.
  if (cell_test != cell_trial)
    return;

  // Get the Quadrature rule for this cell...
  std::vector<double> *xgl, *wgl;
  unsigned int ngl = cell_test.active_degree(0) * 3;
  qcollection.get_weights_and_points(ngl, wgl, xgl);

  auto cell_test_geom = cell_test.get_cell_on_mesh();

  if (cell_test.cell_dim == 2) {
    // Since we don't actually care about the DoFs for this test...
    unsigned int dof_index = 0;
    for (unsigned int i = 0; i < wgl->size(); ++i) {
      double interm = 0.;
      for (unsigned int j = 0; j < wgl->size(); ++j) {
        // Evaluate the Jacobian...
        std::cout << "does it make it to uniform integrator?\n";
        auto jac_value = cell_test_geom->jacobian({xgl->at(i), xgl->at(j)});
        interm += jac_value * wgl->at(j);
      }
      output->at(dof_index, dof_index) += wgl->at(i) * interm;
    }
  } else
    assert(false && "Not implemented for dim != 2");
}
template <class DoFCellType, class CoefficientType, class Real>
void UniformIntegrator<DoFCellType, CoefficientType, Real>::integrate(
    const DoFCellType &cell_test, const DoFCellType &cell_trial,
    const DoFMask &mask_test, const DoFMask &mask_trial,
    DenseSubMatrix<CoefficientType> *output) {
  integrate_regular(cell_test, cell_trial, mask_test, mask_trial, output);
}
template <class DoFCellType, class CoefficientType, class Real>
void UniformIntegrator<DoFCellType, CoefficientType, Real>::
    integrate_near_singular(
        const DoFCellType &cell_test, const DoFCellType &cell_trial,
        const DoFMask &mask_test, const DoFMask &mask_trial,
        const RegularityType &regularity,
        const MultiIndex<2, unsigned int> &regularity_local_index,
        DenseSubMatrix<CoefficientType> *output) {
  assert(false && "Not implemented for this integrator!");
}
template <class DoFCellType, class CoefficientType, class Real>
void UniformIntegrator<DoFCellType, CoefficientType, Real>::integrate_singular(
    const DoFCellType &cell_test, const DoFCellType &cell_trial,
    const DoFMask &mask_test, const DoFMask &mask_trial,
    const RegularityType &regularity,
    const MultiIndex<2, unsigned int> &regularity_local_index,
    DenseSubMatrix<CoefficientType> *output) {
  assert(false && "Not implemented for this integrator!");
}
template <class DoFCellType, class CoefficientType, class Real>
void UniformIntegrator<DoFCellType, CoefficientType, Real>::integrate(
    const DoFCellType &cell_test, const DoFCellType &cell_trial,
    const DoFMask &mask_test, const DoFMask &mask_trial,
    const MaterialData<Real> &material_data,
    const Excitations::Excitation<DoFCellType::space_dim, Real> &excitation,
    DenseSubMatrix<CoefficientType> *output, const unsigned int& ngl_order)
{
  assert(false && "Not implemented for this integrator!");
}
template <class DoFCellType, class CoefficientType, class Real>
void UniformIntegrator<DoFCellType, CoefficientType, Real>::fill_excitation(
    const DoFCellType &cell_test, const DoFMask &mask_test,
    const Excitations::Excitation<DoFCellType::space_dim, Real> &excitation,
    const unsigned int& ngl, DenseSubVector<CoefficientType> *output)
{
    // Integrate over the surface of the cell
        assert(false && "Not implemented for this integrator!");
}

DROMON_NAMESPACE_CLOSE

#endif // DROMON_UNIFORMINTEGRATOR_H
