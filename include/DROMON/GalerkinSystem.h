//
// Created by Jake J. Harmon (jake.harmon@ieee.org) on 3/30/22.
//

#ifndef DROMON_GALERKINSYSTEM_H
#define DROMON_GALERKINSYSTEM_H

#include "DoFHandler.h"
#include "DoFMask.h"
#include "SubMatrix.h"
#include "config.h"

DROMON_NAMESPACE_OPEN

template <unsigned int dim, unsigned int spacedim,
          class CoefficientType = std::complex<double>, class Real = double>
class GalerkinSystem {
public:
  using DoFHandlerType = DoFHandler<dim, spacedim, CoefficientType, Real>;

  explicit GalerkinSystem(
      DoFHandler<dim, spacedim, CoefficientType, Real> *dof_handler);
  void update_system_size();
  void replace_dof_handler(
      DoFHandler<dim, spacedim, CoefficientType, Real> *dof_handler);

  template <class Integrator>
  void fill_system(Integrator &integrator,
                   const MaterialData<Real> &materialData,
                   const Excitations::Excitation<spacedim, Real> &excitation,
                   const unsigned int &ngl_order = 0, const unsigned int& ngl_vertex_order = 0, const unsigned int& ngl_edge_order = 0, const unsigned int& ngl_self_order = 0);

  template <class Integrator>
  void fill_system_HOPS(Integrator &integrator,
                   const MaterialData<Real> &materialData,
                   const Excitations::Excitation<spacedim, Real> &excitation,
                   const unsigned int &ngl_order = 0, const unsigned int& ngl_vertex_order = 0, const unsigned int& ngl_edge_order = 0, const unsigned int& ngl_self_order = 0);

  template <class Integrator>
  void fill_system_perturb(Integrator &integrator,
                   const MaterialData<Real> &materialData,
                   const Excitations::Excitation<spacedim, Real> &excitation, int perturbVar, double perturbSize,
                   const unsigned int &ngl_order = 0, const unsigned int& ngl_vertex_order = 0, const unsigned int& ngl_edge_order = 0, const unsigned int& ngl_self_order = 0);


  template <class Integrator> void fill_system(Integrator &integrator);
  template <class Integrator>
  void fill_forward_excitation(
      Integrator &integrator,
      const Excitations::Excitation<spacedim, Real> &excitation,
      const MaterialData<Real> &materialData, const unsigned int &ngl_order = 5);

  template <class Integrator>
  void fill_forward_excitation_HOPS(
      Integrator &integrator,
      const Excitations::Excitation<spacedim, Real> &excitation,
      const MaterialData<Real> &materialData, const unsigned int &ngl_order = 5);

  template <class Integrator>
  void fill_forward_excitation_perturb(
      Integrator &integrator,
      const Excitations::Excitation<spacedim, Real> &excitation,
      const MaterialData<Real> &materialData, int perturbVar, double perturbSize, const unsigned int &ngl_order = 5);

  template <class AdjointExcitation>
  void fill_adjoint_excitation(AdjointExcitation &adjointExcitation,
                               const unsigned int &ngl_order = 5);

  //  template <class Excitation>
  //  void fill_excitation(Excitation& excitation);

  DenseSubMatrix<CoefficientType> &subsystem_at(const unsigned int &i);
  DenseSubVector<CoefficientType> &subexcitation_at(const unsigned int &i);
  DenseSubVector<CoefficientType> &
  adjoint_subexcitation_at(const unsigned int &i);

  const DenseSubMatrix<CoefficientType> &
  const_subsystem_at(const unsigned int &i) const;
  const DenseSubVector<CoefficientType> &
  const_subexcitation_at(const unsigned int &i) const;

  DoFHandler<dim, spacedim, CoefficientType, Real> *get_dof_handler();
  bool is_symmetric_system() const;
  unsigned int size() const;

  void set_mask_from_dof_activity();

private:
  bool is_symmetric = false;
  std::vector<std::unique_ptr<DenseSubMatrix<CoefficientType>>> cell_subsystem;
  std::vector<std::unique_ptr<DenseSubVector<CoefficientType>>>
      cell_subexcitation;
  std::vector<std::unique_ptr<DenseSubVector<CoefficientType>>>
      cell_adjoint_subexcitation;
  std::vector<std::unique_ptr<DoFMask>> cell_dofmask;

  DoFHandler<dim, spacedim, CoefficientType, Real> *dof_handler;
};
template <unsigned int dim, unsigned int spacedim, class CoefficientType,
          class Real>
GalerkinSystem<dim, spacedim, CoefficientType, Real>::GalerkinSystem(
    DoFHandler<dim, spacedim, CoefficientType, Real> *dof_handler) {
  this->dof_handler = dof_handler;
  for (auto &cell : (this->dof_handler)->get_dof_parents()) {
    cell_subsystem.push_back(std::make_unique<DenseSubMatrix<CoefficientType>>(
        DenseSubMatrix<CoefficientType>(
            this->dof_handler->get_n_dofs_on_cell(cell.index),
            this->dof_handler->get_n_dofs())));

    cell_subexcitation.push_back(
        std::make_unique<DenseSubVector<CoefficientType>>(
            DenseSubVector<CoefficientType>(
                this->dof_handler->get_n_dofs_on_cell(cell.index))));

    cell_adjoint_subexcitation.push_back(
        std::make_unique<DenseSubVector<CoefficientType>>(
            DenseSubVector<CoefficientType>(
                this->dof_handler->get_n_dofs_on_cell(cell.index))));

    cell_dofmask.push_back(std::make_unique<DoFMask>(
        this->dof_handler->get_n_dofs_on_cell(cell.index)));
  }
}

template <unsigned int dim, unsigned int spacedim, class CoefficientType,
          class Real>
void GalerkinSystem<dim, spacedim, CoefficientType,
                    Real>::update_system_size() {
  //  for (auto& cell : cell_subsystem)
  //    cell->resize(this->dof_handler->n_dofs_on_cell(cell.index),
  //    this->dof_handler->n_dofs());
  //
  //  for (auto& cell : cell_dofmask)
  //    cell->resize(this->dof_handler->n_dofs_on_cell(cell.index));

  for (unsigned int i = 0; i < cell_subsystem.size(); ++i) {
    const auto n_dofs_on_cell =
        std::max(this->dof_handler->get_n_dofs_on_cell(i),
                 this->cell_subsystem[i]->n_rows());
    const auto n_total_dofs = std::max(this->dof_handler->get_n_dofs(),
                                       this->cell_subsystem[i]->n_cols());

    this->cell_subsystem[i]->resize(n_dofs_on_cell, n_total_dofs);
    this->cell_dofmask[i]->resize(n_dofs_on_cell);
    this->cell_subexcitation[i]->resize(n_dofs_on_cell);
    this->cell_adjoint_subexcitation[i]->resize(n_dofs_on_cell);
  }
}

template <unsigned int dim, unsigned int spacedim, class CoefficientType,
          class Real>
template <class Integrator>
void GalerkinSystem<dim, spacedim, CoefficientType, Real>::fill_system(
    Integrator &integrator, const MaterialData<Real> &materialData,
    const Excitations::Excitation<spacedim, Real> &excitation,
    const unsigned int &ngl_order, const unsigned int& ngl_vertex_order, const unsigned int& ngl_edge_order, const unsigned int& ngl_self_order) {
  const unsigned int uniform_ngl = ngl_order == 0 ? 5 : ngl_order;
  is_symmetric = integrator.is_symmetric();

  //std::cout << "this fill system is being called\n";
// Originally version integrating all cells simultaneously
  // However, the different interactions take different amounts of time, so we
  // should take care of that
//#pragma omp parallel for //collapse(2)
//  for (auto &cell_test : (this->dof_handler)->get_dof_parents()) {
//    for (auto &cell_trial : (this->dof_handler)->get_dof_parents()) {
//      //      ////////////////////////////
//      //      // Temporary: for debugging!
//      //      if (cell_trial.index == cell_test.index)
//      //        continue;
//      //      // Temporary: for debugging!
//      //      ////////////////////////////
//
//      if (integrator.is_symmetric() && cell_trial.index < cell_test.index)
//        continue;
//
//      integrator.integrate(
//          cell_test, cell_trial, *cell_dofmask[cell_test.index],
//          *cell_dofmask[cell_trial.index], materialData, excitation,
//          cell_subsystem[cell_test.index].get(), ngl_order);
//    }
//  }

  // Integrate self terms
#pragma omp parallel for
  for (auto &cell_test : (this->dof_handler)->get_dof_parents()) {
    auto & cell_trial = cell_test;
      //      ////////////////////////////
      //      // Temporary: for debugging!
      //      if (cell_trial.index == cell_test.index)
      //        continue;
      //      // Temporary: for debugging!
      //      ////////////////////////////
      integrator.integrate(
          cell_test, cell_trial, *cell_dofmask[cell_test.index],
          *cell_dofmask[cell_trial.index], materialData, excitation,
          cell_subsystem[cell_test.index].get(), ngl_order, ngl_vertex_order, ngl_edge_order, ngl_self_order);
  }

#pragma omp parallel for collapse(2) schedule(dynamic)
  for (auto &cell_test : (this->dof_handler)->get_dof_parents()) {
    for (auto &cell_trial : (this->dof_handler)->get_dof_parents()) {

      if (cell_trial.index == cell_test.index || (integrator.is_symmetric() && cell_trial.index < cell_test.index))
        continue;

      integrator.integrate(
          cell_test, cell_trial, *cell_dofmask[cell_test.index],
          *cell_dofmask[cell_trial.index], materialData, excitation,
          cell_subsystem[cell_test.index].get(), ngl_order, ngl_vertex_order, ngl_edge_order, ngl_self_order);
    }
  }


}

template <unsigned int dim, unsigned int spacedim, class CoefficientType,
          class Real>
template <class Integrator>
void GalerkinSystem<dim, spacedim, CoefficientType, Real>::fill_system_HOPS(
    Integrator &integrator, const MaterialData<Real> &materialData,
    const Excitations::Excitation<spacedim, Real> &excitation,
    const unsigned int &ngl_order, const unsigned int& ngl_vertex_order, const unsigned int& ngl_edge_order, const unsigned int& ngl_self_order) {
  const unsigned int uniform_ngl = ngl_order == 0 ? 5 : ngl_order;
  is_symmetric = integrator.is_symmetric();

// Originally version integrating all cells simultaneously
  // However, the different interactions take different amounts of time, so we
  // should take care of that
//#pragma omp parallel for //collapse(2)
//  for (auto &cell_test : (this->dof_handler)->get_dof_parents()) {
//    for (auto &cell_trial : (this->dof_handler)->get_dof_parents()) {
//      //      ////////////////////////////
//      //      // Temporary: for debugging!
//      //      if (cell_trial.index == cell_test.index)
//      //        continue;
//      //      // Temporary: for debugging!
//      //      ////////////////////////////
//
//      if (integrator.is_symmetric() && cell_trial.index < cell_test.index)
//        continue;
//
//      integrator.integrate(
//          cell_test, cell_trial, *cell_dofmask[cell_test.index],
//          *cell_dofmask[cell_trial.index], materialData, excitation,
//          cell_subsystem[cell_test.index].get(), ngl_order);
//    }
//  }

  // Integrate self terms
#pragma omp parallel for
  for (auto &cell_test : (this->dof_handler)->get_dof_parents()) {
    auto & cell_trial = cell_test;
      //      ////////////////////////////
      //      // Temporary: for debugging!
      //      if (cell_trial.index == cell_test.index)
      //        continue;
      //      // Temporary: for debugging!
      //      ////////////////////////////
      integrator.integrate_HOPS(
          cell_test, cell_trial, *cell_dofmask[cell_test.index],
          *cell_dofmask[cell_trial.index], materialData, excitation,
          cell_subsystem[cell_test.index].get(), ngl_order, ngl_vertex_order, ngl_edge_order, ngl_self_order);
  }

#pragma omp parallel for collapse(2) schedule(dynamic)
  for (auto &cell_test : (this->dof_handler)->get_dof_parents()) {
    for (auto &cell_trial : (this->dof_handler)->get_dof_parents()) {

      if (cell_trial.index == cell_test.index || (integrator.is_symmetric() && cell_trial.index < cell_test.index))
        continue;

      integrator.integrate_HOPS(
          cell_test, cell_trial, *cell_dofmask[cell_test.index],
          *cell_dofmask[cell_trial.index], materialData, excitation,
          cell_subsystem[cell_test.index].get(), ngl_order, ngl_vertex_order, ngl_edge_order, ngl_self_order);
    }
  }


}

template <unsigned int dim, unsigned int spacedim, class CoefficientType,
          class Real>
template <class Integrator>
void GalerkinSystem<dim, spacedim, CoefficientType, Real>::fill_system_perturb(
    Integrator &integrator, const MaterialData<Real> &materialData,
    const Excitations::Excitation<spacedim, Real> &excitation,  int perturbVar, double perturbSize,
    const unsigned int &ngl_order, const unsigned int& ngl_vertex_order, const unsigned int& ngl_edge_order, const unsigned int& ngl_self_order) {
  const unsigned int uniform_ngl = ngl_order == 0 ? 5 : ngl_order;
  is_symmetric = integrator.is_symmetric();

  // Integrate self terms
#pragma omp parallel for
  for (auto &cell_test : (this->dof_handler)->get_dof_parents()) {
    auto & cell_trial = cell_test;
      //      ////////////////////////////
      //      // Temporary: for debugging!
      //      if (cell_trial.index == cell_test.index)
      //        continue;
      //      // Temporary: for debugging!
      //      ////////////////////////////
      integrator.integrate_perturb(
          cell_test, cell_trial, *cell_dofmask[cell_test.index],
          *cell_dofmask[cell_trial.index], materialData, excitation, perturbVar, perturbSize,
          cell_subsystem[cell_test.index].get(), ngl_order, ngl_vertex_order, ngl_edge_order, ngl_self_order);
  }

#pragma omp parallel for collapse(2) schedule(dynamic)
  for (auto &cell_test : (this->dof_handler)->get_dof_parents()) {
    for (auto &cell_trial : (this->dof_handler)->get_dof_parents()) {

      if (cell_trial.index == cell_test.index || (integrator.is_symmetric() && cell_trial.index < cell_test.index))
        continue;

      integrator.integrate_perturb(
          cell_test, cell_trial, *cell_dofmask[cell_test.index],
          *cell_dofmask[cell_trial.index], materialData, excitation, perturbVar, perturbSize,
          cell_subsystem[cell_test.index].get(), ngl_order, ngl_vertex_order, ngl_edge_order, ngl_self_order);
    }
  }


}


template <unsigned int dim, unsigned int spacedim, class CoefficientType,
          class Real>
template <class Integrator>
void GalerkinSystem<dim, spacedim, CoefficientType, Real>::fill_system(
    Integrator &integrator) {

  std::cout << "this happens now\n";
  for (auto &cell_test : (this->dof_handler)->get_dof_parents())
    for (auto &cell_trial : (this->dof_handler)->get_dof_parents()) {
      if (integrator.is_symmetric() && cell_trial.index < cell_test.index)
        break;

      integrator.integrate(cell_test, cell_trial,
                           *cell_dofmask[cell_test.index],
                           *cell_dofmask[cell_trial.index],
                           cell_subsystem[cell_test.index].get());
    }
}



template <unsigned int dim, unsigned int spacedim, class CoefficientType,
          class Real>
DenseSubMatrix<CoefficientType> &
GalerkinSystem<dim, spacedim, CoefficientType, Real>::subsystem_at(
    const unsigned int &i) {
  assert(i < cell_subsystem.size() && "Index i out of range!");

  return *cell_subsystem[i];
}
template <unsigned int dim, unsigned int spacedim, class CoefficientType,
          class Real>
DenseSubVector<CoefficientType> &
GalerkinSystem<dim, spacedim, CoefficientType, Real>::subexcitation_at(
    const unsigned int &i) {
  assert(i < cell_subexcitation.size() && "Index i out of range!");

  return *cell_subexcitation[i];
}
template <unsigned int dim, unsigned int spacedim, class CoefficientType,
          class Real>
unsigned int
GalerkinSystem<dim, spacedim, CoefficientType, Real>::size() const {
  return this->cell_subsystem.size();
}
template <unsigned int dim, unsigned int spacedim, class CoefficientType,
          class Real>
DoFHandler<dim, spacedim, CoefficientType, Real> *
GalerkinSystem<dim, spacedim, CoefficientType, Real>::get_dof_handler() {
  return dof_handler;
}
template <unsigned int dim, unsigned int spacedim, class CoefficientType,
          class Real>
bool GalerkinSystem<dim, spacedim, CoefficientType, Real>::is_symmetric_system()
    const {
  return is_symmetric;
}
template <unsigned int dim, unsigned int spacedim, class CoefficientType,
          class Real>
const DenseSubMatrix<CoefficientType> &
GalerkinSystem<dim, spacedim, CoefficientType, Real>::const_subsystem_at(
    const unsigned int &i) const {
  assert(i < cell_subsystem.size() && "Index i out of range!");

  return *cell_subsystem[i];
}
template <unsigned int dim, unsigned int spacedim, class CoefficientType,
          class Real>
const DenseSubVector<CoefficientType> &
GalerkinSystem<dim, spacedim, CoefficientType, Real>::const_subexcitation_at(
    const unsigned int &i) const {
  assert(i < cell_subexcitation.size() && "Index i out of range!");

  return *cell_subexcitation[i];
}
template <unsigned int dim, unsigned int spacedim, class CoefficientType,
          class Real>
void GalerkinSystem<dim, spacedim, CoefficientType, Real>::replace_dof_handler(
    DoFHandler<dim, spacedim, CoefficientType, Real> *dof_handler) {
  this->dof_handler = dof_handler;
}
template <unsigned int dim, unsigned int spacedim, class CoefficientType,
          class Real>
template <class Integrator>
void GalerkinSystem<dim, spacedim, CoefficientType, Real>::
    fill_forward_excitation(
        Integrator &integrator,
        const Excitations::Excitation<spacedim, Real> &excitation,
        const MaterialData<Real> &materialData, const unsigned int &ngl_order) {
#ifndef DEBUG
#pragma omp parallel for
#endif
  for (auto &cell_test : (this->dof_handler)->get_dof_parents()) {
    cell_subexcitation[cell_test.index]->zero_out();
    integrator.fill_excitation(cell_test, *cell_dofmask[cell_test.index],
                               excitation, ngl_order,
                               cell_subexcitation[cell_test.index].get());
  }
}

template <unsigned int dim, unsigned int spacedim, class CoefficientType,
          class Real>
template <class Integrator>
void GalerkinSystem<dim, spacedim, CoefficientType, Real>::
    fill_forward_excitation_HOPS(
        Integrator &integrator,
        const Excitations::Excitation<spacedim, Real> &excitation,
        const MaterialData<Real> &materialData, const unsigned int &ngl_order) {
#ifndef DEBUG
#pragma omp parallel for
#endif
  for (auto &cell_test : (this->dof_handler)->get_dof_parents()) {
    cell_subexcitation[cell_test.index]->zero_out();
    integrator.fill_excitation_HOPS(cell_test, *cell_dofmask[cell_test.index],
                               excitation, ngl_order, materialData,
                               cell_subexcitation[cell_test.index].get());
  }
}

template <unsigned int dim, unsigned int spacedim, class CoefficientType,
          class Real>
template <class Integrator>
void GalerkinSystem<dim, spacedim, CoefficientType, Real>::
    fill_forward_excitation_perturb(
        Integrator &integrator,
        const Excitations::Excitation<spacedim, Real> &excitation,
        const MaterialData<Real> &materialData, int perturbVar, double perturbSize, const unsigned int &ngl_order) {
#ifndef DEBUG
#pragma omp parallel for
#endif
  for (auto &cell_test : (this->dof_handler)->get_dof_parents()) {
    cell_subexcitation[cell_test.index]->zero_out();
    integrator.fill_excitation_perturb(cell_test, *cell_dofmask[cell_test.index],
                               excitation, ngl_order, materialData, perturbVar, perturbSize,
                               cell_subexcitation[cell_test.index].get());
  }
}

template <unsigned int dim, unsigned int spacedim, class CoefficientType,
          class Real>
template <class AdjointExcitation>
void GalerkinSystem<dim, spacedim, CoefficientType, Real>::
    fill_adjoint_excitation(AdjointExcitation &adjointExcitation,
                            const unsigned int &ngl_order) {
#ifndef DEBUG
#pragma omp parallel for
#endif
  for (auto &cell_test : (this->dof_handler)->get_dof_parents()) {
    cell_adjoint_subexcitation[cell_test.index]->zero_out();
    adjointExcitation.fill_excitation(
        cell_test, *cell_dofmask[cell_test.index], ngl_order,
        cell_adjoint_subexcitation[cell_test.index].get());
  }
}
template <unsigned int dim, unsigned int spacedim, class CoefficientType,
          class Real>
DenseSubVector<CoefficientType> &
GalerkinSystem<dim, spacedim, CoefficientType, Real>::adjoint_subexcitation_at(
    const unsigned int &i) {
  assert(i < cell_adjoint_subexcitation.size() && "Index i out of range!");

  return *cell_adjoint_subexcitation[i];
}
template <unsigned int dim, unsigned int spacedim, class CoefficientType,
          class Real>
void GalerkinSystem<dim, spacedim, CoefficientType,
                    Real>::set_mask_from_dof_activity() {
  // Iterate through each cell and disable those DoFs on the DoFMask
  // that are not active
  //unsigned int cell_index = 0;
  for (const auto &cell : (this->dof_handler)->get_dof_parents()) {
    for (unsigned int dof_index = 0; dof_index < cell.n_dofs(); ++dof_index) {
      if (!cell.get_dof(dof_index).is_active)
        (this->cell_dofmask[cell.index])->mask_on(dof_index);
      else
        (this->cell_dofmask[cell.index])->mask_off(dof_index);
    }
  //  ++cell_index;
  }
}

DROMON_NAMESPACE_CLOSE

#endif // DROMON_GALERKINSYSTEM_H