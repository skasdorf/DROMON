//
// Created by Jake J. Harmon (jake.harmon@ieee.org) on 3/30/22.
//

#ifndef DROMON_INTEGRATORBASE_H
#define DROMON_INTEGRATORBASE_H

#include "DoFMask.h"
#include "Excitations.h"
#include "MultiIndex.h"
#include "SubMatrix.h"
#include "config.h"

DROMON_NAMESPACE_OPEN

template <class DoFCellType, class CoefficientType, class Real = double>
class IntegratorBase {
public:
  IntegratorBase() = default;
  // IntegratorBase(FE* fe) = default;
  virtual bool is_symmetric() const = 0;
  virtual const std::string get_integrator_name() const = 0;
  // virtual void attach_quadrature_collection() = 0;
  virtual void integrate(const DoFCellType &cell_test,
                         const DoFCellType &cell_trial,
                         const DoFMask &mask_test, const DoFMask &mask_trial,
                         DenseSubMatrix<CoefficientType> *output) = 0;
  virtual void integrate(
      const DoFCellType &cell_test, const DoFCellType &cell_trial,
      const DoFMask &mask_test, const DoFMask &mask_trial,
      const MaterialData<Real> &material_data,
      const Excitations::Excitation<DoFCellType::space_dim, Real> &excitation,
      DenseSubMatrix<CoefficientType> *output,
      const unsigned int &ngl_order = 0, const unsigned int& ngl_vertex_order = 0, const unsigned int& ngl_edge_order = 0, const unsigned int& ngl_self_order = 0) = 0;
  virtual void fill_excitation(
      const DoFCellType &cell_test, const DoFMask &mask_test,
      const Excitations::Excitation<DoFCellType::space_dim, Real> &excitation,
      const unsigned int &ngl, DenseSubVector<CoefficientType> *output) = 0;

private:
  //  virtual void integrate_near_singular(const DoFCellType& cell_test, const
  //  DoFCellType& cell_trial, const DoFMask& mask_test, const DoFMask&
  //  mask_trial, const RegularityType& regularity, const MultiIndex<2, unsigned
  //  int>& regularity_local_index, DenseSubMatrix<CoefficientType>* output) =
  //  0; virtual void integrate_singular(const DoFCellType& cell_test, const
  //  DoFCellType& cell_trial, const DoFMask& mask_test, const DoFMask&
  //  mask_trial, const RegularityType& regularity, const MultiIndex<2, unsigned
  //  int>& regularity_local_index, DenseSubMatrix<CoefficientType>* output) =
  //  0; virtual void integrate_regular(const DoFCellType& cell_test, const
  //  DoFCellType& cell_trial, const DoFMask& mask_test, const DoFMask&
  //  mask_trial, DenseSubMatrix<CoefficientType>* output) = 0;
};

DROMON_NAMESPACE_CLOSE

#endif // DROMON_INTEGRATORBASE_H
