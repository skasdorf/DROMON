//
// Created by Jake J. Harmon (jake.harmon@ieee.org) on 5/7/22.
//

#ifndef DROMON_REFINEMENT_H
#define DROMON_REFINEMENT_H

//#include "DoFHandler.h"
#include "config.h"
//#include "utility.h"

DROMON_NAMESPACE_OPEN
namespace Refinement {
template <class DoFHandlerType>
void set_degrees(DoFHandlerType *dof_handler, const std::vector<unsigned int>& n_refinements)
{
  auto fe_collector_collector = dof_handler->get_fe_collection_collector();
  for (const auto &cur_type : fe_collector_collector->get_current_types()) {
    for (auto &cell : dof_handler->get_dof_parents())
    {
      const unsigned int& n_refine = n_refinements[cell.index];
      cell.set_fe_index(
          fe_collector_collector->get_fe_collection(cur_type)->next_fe(cell.get_current_fe_index(cur_type), n_refine), cur_type);
    }
    }
}

template <class DoFHandlerType>
void uniform_p_refinement(DoFHandlerType *dof_handler) {
  auto fe_collector_collector = dof_handler->get_fe_collection_collector();
  for (const auto &cur_type : fe_collector_collector->get_current_types()) {
    for (auto &cell : dof_handler->get_dof_parents()) {
      cell.set_fe_index(
          fe_collector_collector->get_fe_collection(cur_type)->next_fe(
              cell.get_current_fe_index(cur_type)),
          cur_type);
    }
  }
}

template <class DoFHandlerType>
void prediction_p_refinement(
    DoFHandlerType *dof_handler,
    const std::vector<typename DoFHandlerType::coefficient_type> &cell_errors,
    const typename DoFHandlerType::real_type &abs_local_tol,
    const int p_adjust = 0) {
  using Real = typename DoFHandlerType::real_type;

  auto fe_collector_collector = dof_handler->get_fe_collection_collector();
  for (auto &cell : dof_handler->get_dof_parents()) {
    // Get the cell on the Geom
    auto cell_geom = cell.get_cell_on_mesh();
    const Real diameter = cell_geom->diameter();
    const auto &current_error = cell_errors[cell.index];
    for (const auto &cur_type : fe_collector_collector->get_current_types()) {
      // At the moment, we conduct the same refinement regardless of the
      // numbers of current types. However, ideally the function should be
      // flexible to control refinements for Electric and Magnetic currents
      const auto active_degree = cell.active_degree(cur_type) + p_adjust;
      const auto abs_current_error = std::abs(current_error);

      Real hi_scal = (abs_local_tol / abs_current_error);
      if (hi_scal < 1.0) // Needs Refinement
      {
        Real newton_p = 0;
        auto fval = [&active_degree, &abs_current_error,
                     &abs_local_tol](const Real &x) {
          return std::pow(x / Real(active_degree), -x - Real(1.5)) -
                 abs_local_tol / abs_current_error;
        };
        auto fderiv = [&active_degree, &abs_current_error,
                       &abs_local_tol](const Real &x) {
          return -log(x / Real(active_degree)) /
                     std::pow(x / Real(active_degree), x + Real(1.5)) -
                 (x + Real(1.5)) /
                     (Real(active_degree) *
                      (std::pow(x / Real(active_degree), x + Real(2.5))));
        };
        newton_p = utility::newton_raphson(fval, fderiv, Real(active_degree),
                                           Real(0.1));
        int p_refine =
            std::ceil(Real(0.75) * newton_p) - active_degree;
        if (p_refine > 3)
          p_refine = 3;
        else if (p_refine < 1)
          p_refine = 1;
        cell.set_fe_index(
            fe_collector_collector->get_fe_collection(cur_type)->next_fe(
                cell.get_current_fe_index(cur_type), p_refine),
            cur_type);
      } else { // Needs coarsening
        const Real amplifier =
            std::pow(abs_local_tol / abs_current_error / Real(100.0),
                     Real(-1.0) / (active_degree + Real(1.5))) *
            Real(active_degree);
        int p_coarse = active_degree - std::ceil(amplifier) - 1;
        if (p_coarse >= 1)
          cell.set_fe_index(
              fe_collector_collector->get_fe_collection(cur_type)->prev_fe(
                  cell.get_current_fe_index(cur_type), p_coarse),
              cur_type);
      }
    }
  }
}
} // namespace Refinement

DROMON_NAMESPACE_CLOSE

#endif // DROMON_REFINEMENT_H
