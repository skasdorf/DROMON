//
// Created by Jake J. Harmon (jake.harmon@ieee.org) on 5/15/22.
//

#ifndef DROMON_ERRORESTIMATION_H
#define DROMON_ERRORESTIMATION_H

#include "DoFHandler.h"
#include "GalerkinSystem.h"
#include "config.h"

DROMON_NAMESPACE_OPEN
namespace ErrorEstimation {
/**
 * Computes the Dual Weighted Residual Error Contribution estimates given
 * the DoFHandlers to the forward discretization and the adjoint discretization
 * Note that these discretizations share the same GalerkinSystem, which is where
 * the integrals required are pulled from.
 * @tparam dim
 * @tparam spacedim
 * @tparam CoefficientType
 * @tparam Real
 * @param dof_handler_forward
 * @param dof_handler_adjoint
 * @param forward_solution
 * @param adjoint_solution
 * @param cg_sys
 * @param error_contribution_estimates
 */
template <unsigned int dim, unsigned int spacedim, class CoefficientType,
          class Real, class VectorType>
void EstimateDWR(
    const DoFHandler<dim, spacedim, CoefficientType, Real> &dof_handler_forward,
    const DoFHandler<dim, spacedim, CoefficientType, Real> &dof_handler_adjoint,
    const std::vector<CoefficientType> &forward_solution,
    const std::vector<CoefficientType> &adjoint_solution,
    const GalerkinSystem<dim, spacedim, CoefficientType, Real> &cg_sys,
    VectorType &error_contribution_estimates) {
  const bool is_symmetric = cg_sys.is_symmetric_system();
  // Loop through all the cells in the adjoint DoFHandler
  // then loop through all the cells in the forward DoFhandler
  for (auto &cell_adjoint : dof_handler_adjoint.get_dof_parents()) {
    const unsigned int n_dofs_adjoint = cell_adjoint.n_dofs();
    const unsigned int cell_adjoint_index = cell_adjoint.get_index();
    const auto &cell_adjoint_on_forward_disc =
        dof_handler_forward.get_dof_parent(cell_adjoint_index);

    const auto &subexcitation =
        cg_sys.const_subexcitation_at(cell_adjoint_index);
    const auto &subsystem = cg_sys.const_subsystem_at(cell_adjoint_index);
    // Loop through all the DoFs
    for (unsigned adjoint_dof_index = 0; adjoint_dof_index < n_dofs_adjoint;
         ++adjoint_dof_index) {
      unsigned int contribution_assigment_location = cell_adjoint_index;
      const auto &dof_adjoint = cell_adjoint.get_dof(adjoint_dof_index);
      if (!dof_adjoint.is_active)
        continue;
      const auto &dof_adjoint_on_forward =
          cell_adjoint_on_forward_disc.get_dof(adjoint_dof_index);
      // Now we must check whether or not this DoF is a higher order DoF, i.e.,
      // if it is active on the forward DoFHandler, we skip
            if (dof_adjoint_on_forward.is_active)
              continue;
      // Now we have selected an appropriate DoF to start accumulating
      // contributions for

      // If the DoF is an Edge/Face DoF, then we need to check if its neighbor
      // has a higher or lower expansion order whichever is lower receives the
      // contributions from that DoF as it is limiting further refinements
      if (dof_adjoint.dof_type == FaceDoF) {
        const auto &neighbor_face = dof_handler_adjoint.get_dof_face(
            dof_adjoint.host->get_neighbor_index());
        const auto &parent_cell_index = neighbor_face.parent_cell_index();
        const auto &neighbor_cell =
            dof_handler_adjoint.get_dof_parent(parent_cell_index);
        const auto &neighbor_degree =
            neighbor_cell.active_degree(dof_adjoint.cur_type);
        if (neighbor_degree < cell_adjoint.active_degree(dof_adjoint.cur_type))
          contribution_assigment_location = parent_cell_index;
      }

      // Now we add the contribution from the inner-product of the forward
      // excitation with the adjoint solution
      error_contribution_estimates.accumulate(
          contribution_assigment_location,
          std::conj(adjoint_solution[dof_adjoint.active_index]) *
              subexcitation.const_at(dof_adjoint.global_index));

      for (auto &cell_forward : dof_handler_forward.get_dof_parents()) {
        const unsigned int cell_forward_index = cell_forward.get_index();
        const unsigned int n_dofs_forward = cell_forward.n_dofs();
        const auto &subsystem_on_forward_cell =
            cg_sys.const_subsystem_at(cell_forward_index);
        const auto &cell_forward_on_adjoint_disc =
            dof_handler_adjoint.get_dof_parent(cell_forward_index);
        // Now we loop through the DoFs of the forward solution, taking only
        // those that are active
        for (unsigned int forward_dof_index = 0;
             forward_dof_index < n_dofs_forward; ++forward_dof_index) {
          const auto &dof_forward = cell_forward.get_dof(forward_dof_index);
          const auto &dof_forward_on_adjoint_disc =
              cell_forward_on_adjoint_disc.get_dof(forward_dof_index);
          if (!dof_forward.is_active)
            continue;

          if (is_symmetric && cell_forward_index < cell_adjoint_index) {
            error_contribution_estimates.accumulate(
                contribution_assigment_location,
                -forward_solution[dof_forward.active_index] *
                    std::conj(adjoint_solution[dof_adjoint.active_index]) *
                    subsystem_on_forward_cell.const_at(
                        dof_forward.global_index, dof_adjoint.global_index));
          } else if (is_symmetric && cell_forward_index == cell_adjoint_index &&
                     dof_forward_on_adjoint_disc.active_index >
                         dof_adjoint.active_index) {
            error_contribution_estimates.accumulate(
                contribution_assigment_location,
                -forward_solution[dof_forward.active_index] *
                    std::conj(adjoint_solution[dof_adjoint.active_index]) *
                    subsystem.const_at(dof_forward.global_index,
                                       dof_adjoint.global_index));
          } else {
            error_contribution_estimates.accumulate(
                contribution_assigment_location,
                -forward_solution[dof_forward.active_index] *
                    std::conj(adjoint_solution[dof_adjoint.active_index]) *
                    subsystem.const_at(dof_adjoint.global_index,
                                       dof_forward.global_index));
          }
        }
      }
    }
  }
  error_contribution_estimates.apply_carry();
  // In the Adjoint DoFHandler, first accumulate those contributions resulting
  // from the inner-product of the forward excitation with the higher order
  // components of the adjoint solution

  // In order to properly decide where contributions should be allocated,
  // in the event of edge DoFs, we assign it to the neighbor that owns the DoF
  // with the lower degree
}

template <unsigned int dim, unsigned int spacedim, class CoefficientType,
          class Real, class VectorType>
void EstimateDWR(
    const DoFHandler<dim, spacedim, CoefficientType, Real> &dof_handler,
    const std::vector<CoefficientType> &forward_solution,
    const std::vector<CoefficientType> &adjoint_solution,
    const GalerkinSystem<dim, spacedim, CoefficientType, Real> &cg_sys,
    VectorType &error_contribution_estimates) {
  const bool is_symmetric = cg_sys.is_symmetric_system();
  // Loop through all the cells in the adjoint DoFHandler
  // then loop through all the cells in the forward DoFhandler
  for (auto &cell_adjoint : dof_handler.get_dof_parents()) {
    const unsigned int n_dofs_adjoint = cell_adjoint.n_dofs();
    const unsigned int cell_adjoint_index = cell_adjoint.get_index();

    const auto &subexcitation =
        cg_sys.const_subexcitation_at(cell_adjoint_index);
    const auto &subsystem = cg_sys.const_subsystem_at(cell_adjoint_index);
    // Loop through all the DoFs
    for (unsigned adjoint_dof_index = 0; adjoint_dof_index < n_dofs_adjoint;
         ++adjoint_dof_index) {
      unsigned int contribution_assigment_location = cell_adjoint_index;
      const auto &dof_adjoint = cell_adjoint.get_dof(adjoint_dof_index);
      if (!dof_adjoint.is_active)
        continue;
      // Now we must check whether or not this DoF is a higher order DoF, i.e.,
      // if it is active on the forward DoFHandler, we skip
      //      if (dof_adjoint_on_forward.is_active)
      //        continue;
      // Now we have selected an appropriate DoF to start accumulating
      // contributions for

      // If the DoF is an Edge/Face DoF, then we need to check if its neighbor
      // has a higher or lower expansion order whichever is lower receives the
      // contributions from that DoF as it is limiting further refinements
      if (dof_adjoint.dof_type == FaceDoF) {
        const auto &neighbor_face =
            dof_handler.get_dof_face(dof_adjoint.host->get_neighbor_index());
        const auto &parent_cell_index = neighbor_face.parent_cell_index();
        const auto &neighbor_cell =
            dof_handler.get_dof_parent(parent_cell_index);
        const auto &neighbor_degree =
            neighbor_cell.active_degree(dof_adjoint.cur_type);
        if (neighbor_degree < cell_adjoint.active_degree(dof_adjoint.cur_type))
          contribution_assigment_location = parent_cell_index;
      }

      // Now we add the contribution from the inner-product of the forward
      // excitation with the adjoint solution
      error_contribution_estimates.accumulate(
          contribution_assigment_location,
          std::conj(adjoint_solution[dof_adjoint.active_index]) *
              subexcitation.const_at(dof_adjoint.global_index));

      for (auto &cell_forward : dof_handler.get_dof_parents()) {
        const unsigned int cell_forward_index = cell_forward.get_index();
        const unsigned int n_dofs_forward = cell_forward.n_dofs();
        const auto &subsystem_on_forward_cell =
            cg_sys.const_subsystem_at(cell_forward_index);

        // Now we loop through the DoFs of the forward solution, taking only
        // those that are active
        for (unsigned int forward_dof_index = 0;
             forward_dof_index < n_dofs_forward; ++forward_dof_index) {
          const auto &dof_forward = cell_forward.get_dof(forward_dof_index);

          if (!dof_forward.is_active)
            continue;

          if (is_symmetric && cell_forward_index < cell_adjoint_index) {
            error_contribution_estimates.accumulate(
                contribution_assigment_location,
                -forward_solution[dof_forward.active_index] *
                    std::conj(adjoint_solution[dof_adjoint.active_index]) *
                    subsystem_on_forward_cell.const_at(
                        dof_forward.global_index, dof_adjoint.global_index));
          } else if (is_symmetric && cell_forward_index == cell_adjoint_index &&
                     dof_forward.active_index > dof_adjoint.active_index) {
            error_contribution_estimates.accumulate(
                contribution_assigment_location,
                -forward_solution[dof_forward.active_index] *
                    std::conj(adjoint_solution[dof_adjoint.active_index]) *
                    subsystem.const_at(dof_forward.global_index,
                                       dof_adjoint.global_index));
          } else {
            error_contribution_estimates.accumulate(
                contribution_assigment_location,
                -forward_solution[dof_forward.active_index] *
                    std::conj(adjoint_solution[dof_adjoint.active_index]) *
                    subsystem.const_at(dof_adjoint.global_index,
                                       dof_forward.global_index));
          }
        }
      }
    }
  }
  error_contribution_estimates.apply_carry();
}

template <unsigned int dim, unsigned int spacedim, class CoefficientType,
          class Real, class MatrixType>
void EstimateDWRGlobalOnly(
    const DoFHandler<dim, spacedim, CoefficientType, Real> &dof_handler_forward,
    const DoFHandler<dim, spacedim, CoefficientType, Real> &dof_handler_adjoint,
    const std::vector<CoefficientType> &forward_solution,
    const std::vector<CoefficientType> &adjoint_solution, MatrixType &matrix,
    CoefficientType &total_error) {
  total_error = 0.0;
  // const bool is_symmetric = cg_sys.is_symmetric_system();
  // Loop through all the cells in the adjoint DoFHandler
  // then loop through all the cells in the forward DoFhandler
  for (auto &cell_adjoint : dof_handler_adjoint.get_dof_parents()) {
    const unsigned int n_dofs_adjoint = cell_adjoint.n_dofs();
    const unsigned int cell_adjoint_index = cell_adjoint.get_index();
    const auto &cell_adjoint_on_forward_disc =
        dof_handler_forward.get_dof_parent(cell_adjoint_index);

    // Loop through all the DoFs
    for (unsigned adjoint_dof_index = 0; adjoint_dof_index < n_dofs_adjoint;
         ++adjoint_dof_index) {
      unsigned int contribution_assigment_location = cell_adjoint_index;
      const auto &dof_adjoint = cell_adjoint.get_dof(adjoint_dof_index);
      if (!dof_adjoint.is_active)
        continue;
      const auto &dof_adjoint_on_forward =
          cell_adjoint_on_forward_disc.get_dof(adjoint_dof_index);
      // Now we must check whether or not this DoF is a higher order DoF, i.e.,
      // if it is active on the forward DoFHandler, we skip
      //      if (dof_adjoint_on_forward.is_active)
      //        continue;

      // Now we add the contribution from the inner-product of the forward
      // excitation with the adjoint solution
      total_error += std::conj(adjoint_solution[dof_adjoint.active_index]) *
                     matrix(dof_adjoint.active_index);
      matrix(dof_adjoint.active_index) = 0;

      for (auto &cell_forward : dof_handler_forward.get_dof_parents()) {
        const unsigned int cell_forward_index = cell_forward.get_index();
        const unsigned int n_dofs_forward = cell_forward.n_dofs();

        const auto &cell_forward_on_adjoint_disc =
            dof_handler_adjoint.get_dof_parent(cell_forward_index);
        // Now we loop through the DoFs of the forward solution, taking only
        // those that are active
        for (unsigned int forward_dof_index = 0;
             forward_dof_index < n_dofs_forward; ++forward_dof_index) {
          const auto &dof_forward = cell_forward.get_dof(forward_dof_index);
          const auto &dof_forward_on_adjoint_disc =
              cell_forward_on_adjoint_disc.get_dof(forward_dof_index);
          if (!dof_forward.is_active)
            continue;

          total_error -= forward_solution[dof_forward.active_index] *
                         std::conj(adjoint_solution[dof_adjoint.active_index]) *
                         matrix(dof_adjoint.active_index,
                                dof_forward_on_adjoint_disc.active_index);
          matrix(dof_adjoint.active_index,
                 dof_forward_on_adjoint_disc.active_index) = 0;
        }
      }
    }
  }
}

} // namespace ErrorEstimation
DROMON_NAMESPACE_CLOSE

#endif // DROMON_ERRORESTIMATION_H
