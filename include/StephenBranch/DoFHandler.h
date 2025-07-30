//
// Created by Jake J. Harmon (jake.harmon@ieee.org) on 8/8/21.
//

#ifndef DROMON_DOFHANDLER_H
#define DROMON_DOFHANDLER_H

#include "DoFBase.h"
#include "DoFGeom.h"
#include "FEBase.h"
#include "FECollection.h"
#include "config.h"
#include "mesh.h"
#include <complex>
#include <vector>

DROMON_NAMESPACE_OPEN
/*
 * Note that for multiple scatterers, wires with plates, etc. a
 * DoFHandlerCollection (in the same vein as the FECollection) would be
 * required, which accepts DoFHandlers of different dim, but same spacedim...
 *
 * */
template <unsigned int dim, unsigned int spacedim,
          class CoefficientType = std::complex<double>, class Real = double>
class DoFHandler {
public:
  using coefficient_type = CoefficientType;
  using real_type = Real;

  DoFHandler() = default;
  explicit DoFHandler(MeshBase<dim, spacedim> &mesh);

  /**
   * distribute_dofs
   * Takes a finite element class instance and applies this template (Note: not
   * a template in C++ sense) to every cell in the discretization attached to
   * this DoFHandler
   */
  void distribute_dofs(const FiniteElement<dim, Real> &fe);

  /**
   * distribute_dofs
   * Takes an FECollection (e.g., different expansions in p, but could also be
   * different types more generally) and applies each to every cell in the
   * discretization attached to this DoFHandler. Unlike the previous function
   * call, which ignores the set index for each cell, this function call
   * requires that each cell in the discretization has been assigned an index
   * (by executing cell->set_fe_index) on each cell. Otherwise, the first entry
   * in the fe_collection is chosen
   *
   */
  void distribute_dofs(
      FECollectionCollector<dim, spacedim, Real> *fe_collection_collector,
      const unsigned int starting_index = 0);

  void distribute_dofs();

  //  void distribute_dofs(
  //      FECollectionCollector<dim, spacedim, Real> *fe_collection_collector,
  //      const std::vector<unsigned int>& starting_indices);

  void spawn_p_reduced_dof_handler(
      DoFHandler<dim, spacedim, CoefficientType, Real> *target);
  std::vector<CoefficientType> map_DoF_values_to_this_handler(
      DoFHandler<dim, spacedim, CoefficientType> *source_dof_handler,
      const std::vector<CoefficientType> &source_values) const;

  std::vector<DoFCell<dim, spacedim>> &get_dof_cells() { return dof_cells; }

  std::vector<DoFParent<dim, spacedim, Real>> &get_dof_parents() {
    return dof_parents;
  }

  const std::vector<DoFParent<dim, spacedim, Real>> &get_dof_parents() const {
    return dof_parents;
  }

  DoFParent<dim, spacedim, Real> &get_dof_parent(const unsigned int &i)

  {
#ifdef DEBUG
    assert(i < dof_parents.size() && "Index i out of range!");
#endif
    return dof_parents[i];
  }

  const DoFParent<dim, spacedim, Real> &
  get_dof_parent(const unsigned int &i) const {
#ifdef DEBUG
    assert(i < dof_parents.size() && "Index i out of range!");
#endif
    return dof_parents[i];
  }

  std::vector<DoFFace<dim, spacedim, Real>> &get_dof_faces() {
    return dof_faces;
  }

  DoFFace<dim, spacedim, Real> &get_dof_face(const unsigned int &i) {
#ifdef DEBUG
    assert(i < dof_faces.size() && "Index i out of range!");
#endif
    return dof_faces[i];
  }
  const DoFFace<dim, spacedim, Real> &
  get_dof_face(const unsigned int &i) const {
#ifdef DEBUG
    assert(i < dof_faces.size() && "Index i out of range!");
#endif
    return dof_faces[i];
  }

  std::vector<DoFFace<dim, spacedim, Real>> get_dof_faces_copy() const {
    return dof_faces;
  }
  std::vector<DoFCell<dim, spacedim>> get_dof_cells_copy() const {
    return dof_cells;
  }

  // Returns the total number of active DoFs associated with a cell, which
  // includes the number of DoFs attached to each of its halfedges...
  const unsigned int get_n_active_dofs_on_cell(unsigned int index) const;

  // Returns the total number of DoFs (both active and inactive) associated with
  // a cell, which includes the number of DoFs attached to each of its
  // halfedges...
  const unsigned int get_n_dofs_on_cell(unsigned int index) const;

  const unsigned int get_n_active_dofs() const;
  const unsigned int get_n_dofs() const;
  const unsigned int get_n_current_types() const;
  FECollectionCollector<dim, spacedim, Real> *get_fe_collection_collector();
  void
  replace_fe_collection_collector(FECollectionCollector<dim, spacedim, Real>
                                      *replacement_fe_collection_collector);
  void replace_mesh(MeshBase<dim, spacedim> *mesh);
  // void refine(const unsigned int& dof_cell_index, const unsigned int&
  // refinement_amount = 1);
private:
  void initialize_dofs();

  void assimilate_new_dofs();
  void finalize_dofs();
  MeshBase<dim, spacedim> *mesh;

  FECollectionCollector<dim, spacedim, Real> *fe_collection_collector;

  /**
   * dof_cells
   * Contains a list of cells (associated with a cell on the attached mesh, but
   * separate) that contain DoFs (both active and inactive). The order in which
   * the vector is accessed matches that of the attached mesh exactly
   */
  std::vector<DoFCell<dim, spacedim>> dof_cells;

  /**
   * dof_faces
   * Note that this handles Faces in the dimension dependent sense (i.e., for a
   * 2-D element, the face is a line (2 points), whereas in 1-D it is a node,
   * etc.
   */
  std::vector<DoFFace<dim, spacedim, Real>> dof_faces;

  /**
   * dof_parents
   * Contains a list of DoFParents associated with each cell in the
   * discretization. The order matches that of the cells in the mesh, but
   * includes the unification of the DoF information contained in dof_cells and
   * dof_faces...
   */
  std::vector<DoFParent<dim, spacedim, Real>> dof_parents;

  unsigned int n_dofs = 0;
  unsigned int n_active_dofs = 0;
};

template <unsigned int dim, unsigned int spacedim, class CoefficientType,
          class Real>
DoFHandler<dim, spacedim, CoefficientType, Real>::DoFHandler(
    MeshBase<dim, spacedim> &mesh) {
  this->mesh = &mesh;

  // Given the attached mesh, generate the DoFCell and DoFFace vectors...
  for (unsigned int i = 0; i < this->mesh->n_cells(); ++i) {
    auto cell_pt = this->mesh->get_pointer_to_cell(i);
    dof_cells.push_back(DoFCell<dim, spacedim>(i, cell_pt->is_active(),
                                               cell_pt->level(), cell_pt));

    dof_parents.push_back(DoFParent<dim, spacedim, Real>());
  }

  for (unsigned int i = 0; i < this->mesh->n_faces(); ++i) {
    auto face_pt = this->mesh->get_pointer_to_face(i);
    dof_faces.push_back(DoFFace<dim, spacedim, Real>(
        i, face_pt->is_active(), face_pt->level(), face_pt));
  }
}

// TODO: Support Anisotropy in $p$
template <unsigned int dim, unsigned int spacedim, class CoefficientType,
          class Real>
void DoFHandler<dim, spacedim, CoefficientType, Real>::distribute_dofs(
    FECollectionCollector<dim, spacedim, Real> *fe_collection_collector,
    const unsigned int starting_index) {
  if (this->fe_collection_collector != fe_collection_collector) {
    this->fe_collection_collector = fe_collection_collector;
    for (auto &cell : dof_cells)
      cell.initialize_fe_indices(
          starting_index, this->fe_collection_collector->get_current_types());
    for (auto &face : dof_faces)
      face.initialize_fe_indices(
          starting_index, this->fe_collection_collector->get_current_types());
    //            //Clear DoF information (if available)
    //            dof_cells.clear();
    //            dof_faces.clear();
    // Now run initialization
    this->initialize_dofs();
  } else {
    this->assimilate_new_dofs();
  }

  // Now that all of the DoFs have been spawned, we need to finalize their
  // insertion by matching them with their shape functions on one cell with its
  // neighbor (if applicable)

  // Secondly, we activate or deactivate DoFs based on the current choice of
  // parameters

  this->finalize_dofs();
}

template <unsigned int dim, unsigned int spacedim, class CoefficientType,
          class Real>
void DoFHandler<dim, spacedim, CoefficientType, Real>::distribute_dofs() {

  #ifdef DEBUG
  assert(this->fe_collection_collector != nullptr &&
         "The FE_collection_collector must be valid!");
  #endif
  this->assimilate_new_dofs();
  this->finalize_dofs();
}

// TODO: Support Anisotropy in $p$
template <unsigned int dim, unsigned int spacedim, class CoefficientType,
          class Real>
void DoFHandler<dim, spacedim, CoefficientType, Real>::initialize_dofs() {
  this->n_dofs = 0;
  // For 2-D:
  //  Expansion order for i, j
  //  Direction (u or v)
  //  Type (electric or magnetic)
  //  Cell index

  switch (dim) {
  case 1:
    //                for (const auto &type: fe_collection.types)
    //                    for (unsigned int i = 0; i <= max_i; ++i)
    break;
  case 2:
    // Loop through the expansion orders first (i,j) (and then, if necessary,
    // the electric vs magnetic) Loop through the cells in the discretization

    for (unsigned int type = 0;
         type < fe_collection_collector->get_number_of_types(); ++type) {
      auto fe_collection_ptr = fe_collection_collector->get_fe_collection(type);
      for (auto &cell : this->dof_cells) {
        cell.attach_active_fe_pointers(fe_collection_collector);
        auto fe_ptr = fe_collection_ptr->get_fe(cell.get_fe_index(type));
        assert(fe_ptr->get_n_components() == dim &&
               "The number of components does not match the dim!");
        for (unsigned int c = 0; c < fe_ptr->get_n_components(); ++c) {
          unsigned int max_i = 0, max_j = 0, min_i = 0, min_j = 0;
          switch (fe_ptr->get_conformity()) {
          case Hdiv:
            if (c == 0) {
              max_i = fe_ptr->get_degree();
              max_j = fe_ptr->get_degree() - 1;
              min_i = 2;
              min_j = 0;
            } else {
              max_i = fe_ptr->get_degree() - 1;
              max_j = fe_ptr->get_degree();
              min_i = 0;
              min_j = 2;
            }
            break;
          case Hcurl:
            if (c == 0) {
              max_i = fe_ptr->get_degree() - 1;
              max_j = fe_ptr->get_degree();
            } else {
              max_i = fe_ptr->get_degree();
              max_j = fe_ptr->get_degree() - 1;
            }
            break;
          default:
            assert(false &&
                   "The conforming space specified is not implemented!");
          }
          for (unsigned int i = min_i; i <= max_i; ++i)
            for (unsigned int j = min_j; j <= max_j; ++j) {
              // add this basis function to the DoFHandler for this cell
              cell.push_back(DoFBase<dim>(
                  &cell, cell.n_dofs(), this->n_dofs++, c == 0 ? u_dir : v_dir,
                  DoFType::CellDoF, fe_ptr->get_type(), MultiIndex<dim>(i, j)));
            }
        }
        cell.mark_as_populated(fe_ptr->get_degree());
      }
      for (auto &face : this->dof_faces) {
        auto fe_ptr = fe_collection_ptr->get_fe(face.get_fe_index(type));
        assert(fe_ptr->get_n_components() == dim &&
               "The number of components does not match the dim!");
        for (unsigned int c = 0; c < fe_ptr->get_n_components(); ++c) {
          unsigned int max_i = 0, max_j = 0, min_i = 0, min_j = 0;
          switch (fe_ptr->get_conformity()) {
          case Hdiv:
            if (c == 0) {
              max_i = 1;
              max_j = fe_ptr->get_degree() - 1;

              min_i = 0;
              min_j = 0;
            } else {
              max_i = fe_ptr->get_degree() - 1;
              max_j = 1;

              min_i = 0;
              min_j = 0;
            }
            break;
          case Hcurl:
            if (c == 0) {
              max_i = fe_ptr->get_degree() - 1;
              max_j = 1;
            } else {
              max_i = 1;
              max_j = fe_ptr->get_degree() - 1;
            }
            break;
          default:
            assert(false &&
                   "The conforming space specified is not implemented!");
          }
          for (unsigned int i = min_i; i <= max_i; ++i)
            for (unsigned int j = min_j; j <= max_j; ++j)
              // add this basis function to the DoFHandler for this cell if
              // compatible
              if (face.is_compatible(c == 0 ? u_dir : v_dir, i, j))
                face.push_back(
                    DoFBase<dim>(&face, face.n_dofs(), this->n_dofs++,
                                 c == 0 ? u_dir : v_dir, DoFType::FaceDoF,
                                 fe_ptr->get_type(), MultiIndex<dim>(i, j)));
        }
        face.mark_as_populated(fe_ptr->get_degree());
      }
    }
    break;
  }
}

template <unsigned int dim, unsigned int spacedim, class CoefficientType,
          class Real>
void DoFHandler<dim, spacedim, CoefficientType, Real>::assimilate_new_dofs() {
  switch (dim) {
  case 1:
    //                for (const auto &type: fe_collection.types)
    //                    for (unsigned int i = 0; i <= max_i; ++i)
    break;
  case 2:
    // Loop through the expansion orders first (i,j) (and then, if necessary,
    // the electric vs magnetic) Loop through the cells in the discretization

    for (unsigned int type = 0;
         type < fe_collection_collector->get_number_of_types(); ++type) {
      auto fe_collection_ptr = fe_collection_collector->get_fe_collection(type);
      for (auto &cell : this->dof_cells) {
        cell.attach_active_fe_pointers(fe_collection_collector);
        if (cell.fully_populated)
          continue;
        auto fe_ptr = fe_collection_ptr->get_fe(cell.get_fe_index(type));
        assert(fe_ptr->get_n_components() == dim &&
               "The number of components does not match the dim!");
        for (unsigned int c = 0; c < fe_ptr->get_n_components(); ++c) {
          unsigned int max_i = 0, max_j = 0, min_i = 0, min_j = 0;
          switch (fe_ptr->get_conformity()) {
          case Hdiv:
            if (c == 0) {
              max_i = fe_ptr->get_degree();
              max_j = fe_ptr->get_degree() - 1;
              min_i = 2;
              min_j = 0;
            } else {
              max_i = fe_ptr->get_degree() - 1;
              max_j = fe_ptr->get_degree();
              min_i = 0;
              min_j = 2;
            }
            break;
          case Hcurl:
            if (c == 0) {
              max_i = fe_ptr->get_degree() - 1;
              max_j = fe_ptr->get_degree();
            } else {
              max_i = fe_ptr->get_degree();
              max_j = fe_ptr->get_degree() - 1;
            }
            break;
          default:
            assert(false &&
                   "The conforming space specified is not implemented!");
          }
          for (unsigned int i = min_i; i <= max_i; ++i)
            for (unsigned int j = min_j; j <= max_j; ++j) {
              MultiIndex<dim> shape_fcn(i, j);
              if (!fe_ptr->is_subset_of_max(shape_fcn, c == 0 ? u_dir : v_dir,
                                            cell.populated_expansion_order))
                // add this basis function to the DoFHandler for this cell
                cell.push_back(
                    DoFBase<dim>(&cell, cell.n_dofs(), this->n_dofs++,
                                 c == 0 ? u_dir : v_dir, DoFType::CellDoF,
                                 fe_ptr->get_type(), MultiIndex<dim>(i, j)));
            }
        }
        cell.mark_as_populated(fe_ptr->get_degree());
      }
      for (auto &face : this->dof_faces) {
        if (face.fully_populated)
          continue;
        auto fe_ptr = fe_collection_ptr->get_fe(face.get_fe_index(type));
        assert(fe_ptr->get_n_components() == dim &&
               "The number of components does not match the dim!");
        for (unsigned int c = 0; c < fe_ptr->get_n_components(); ++c) {
          unsigned int max_i = 0, max_j = 0, min_i = 0, min_j = 0;
          switch (fe_ptr->get_conformity()) {
          case Hdiv:
            if (c == 0) {
              max_i = 1;
              max_j = fe_ptr->get_degree() - 1;

              min_i = 0;
              min_j = 0;
            } else {
              max_i = fe_ptr->get_degree() - 1;
              max_j = 1;

              min_i = 0;
              min_j = 0;
            }
            break;
          case Hcurl:
            if (c == 0) {
              max_i = fe_ptr->get_degree() - 1;
              max_j = 1;
            } else {
              max_i = 1;
              max_j = fe_ptr->get_degree() - 1;
            }
            break;
          default:
            assert(false &&
                   "The conforming space specified is not implemented!");
          }
          for (unsigned int i = min_i; i <= max_i; ++i)
            for (unsigned int j = min_j; j <= max_j; ++j) {
              MultiIndex<dim> shape_fcn(i, j);
              // add this basis function to the DoFHandler for this cell if
              // compatible
              if (face.is_compatible(c == 0 ? u_dir : v_dir, i, j))
                if (!fe_ptr->is_subset_of_max(shape_fcn, c == 0 ? u_dir : v_dir,
                                              face.populated_expansion_order))
                  face.push_back(
                      DoFBase<dim>(&face, face.n_dofs(), this->n_dofs++,
                                   c == 0 ? u_dir : v_dir, DoFType::FaceDoF,
                                   fe_ptr->get_type(), MultiIndex<dim>(i, j)));
            }
        }
        face.mark_as_populated(fe_ptr->get_degree());
      }
    }
    break;
  }
}

/**
 * This function matches all the DoFs and activates and deactivates the
 * other DoFs
 * @tparam dim
 * @tparam spacedim
 * @tparam CoefficientType
 * @tparam Real
 */
template <unsigned int dim, unsigned int spacedim, class CoefficientType,
          class Real>
void DoFHandler<dim, spacedim, CoefficientType, Real>::finalize_dofs() {
  this->n_active_dofs = 0;

  // Deactivate all the DoFs, then activate/reactive the ones that need to be
  for (auto &cell : this->dof_cells)
    cell.deactivate_dofs();

  for (auto &face : this->dof_faces)
    face.deactivate_dofs();
  // Note  that the DoFParent essentially just stores pointers, in which
  // case the overhead is not too significant (the real data is retained in
  // the GalerkinSystem)
  for (auto &parent : this->dof_parents)
    parent.reset();

  // Loop through all cell DoFs and activate only those that satisfy the
  // expansion requirements
  //  in the future, could be extended to the multi-level problem
  for (auto &cell : this->dof_cells) {
    for (auto &dof : (cell.get_dofs())) {
      auto fe_collection_ptr =
          fe_collection_collector->get_fe_collection(dof.cur_type);
      auto fe_ptr = fe_collection_ptr->get_fe(cell.get_fe_index(dof.cur_type));
      if (!fe_ptr->contains_this_dof(dof))
        dof.is_active = false;
      else {
        dof.is_active = true;
        dof.active_index = this->n_active_dofs;
        ++this->n_active_dofs;
      }
    }
    // Now, incorporate these DoFs into the DoFParent associated with this cell
    this->dof_parents[cell.index].attach_dof_cell(&cell);
  }

  // Then Loop through the face dofs
  // Here we have two requirements: that the DoF satisfies the expansion
  // requirements, and secondly, has its twin on a neighboring cell/face

  // In this case, if the neighbor cell does not support the expansion
  // order, simply deactivate and move on. Otherwise, we need to find
  // the partner DoF

  // Note that for matching the DoFs, it matters both the orientation
  // of $u$ and the orientation in $v$

  for (auto &face : this->dof_faces) {
    bool has_sibling = false, orientation = true;
    unsigned int sibling;

    face.get_object_context(&has_sibling, &sibling, &orientation);

    DoFFace<dim, spacedim, Real> *sibling_face;

    if (!has_sibling) {
      face.deactivate_dofs();
      this->dof_parents[face.parent_cell_index()].attach_dof_face(&face);
      continue;
    } else
      sibling_face = &this->dof_faces[sibling];

    for (auto &dof : (face.get_dofs())) {
      // Check if the DoF matches the current settings for the cell
      auto fe_collection_ptr =
          fe_collection_collector->get_fe_collection(dof.cur_type);
      auto fe_ptr = fe_collection_ptr->get_fe(face.get_fe_index(dof.cur_type));

      if (!fe_ptr->contains_this_dof(dof))
        dof.is_active = false;
      else {
        // We need to find if this has its sibling on another cell...
        const auto &dof_dir = dof.direction;
        const auto &dof_cur = dof.cur_type;

        // We only need to match one component, since the other
        // is simply fixed over the edge shape function

        unsigned int dof_match_order =
            dof_dir == u_dir ? dof.orders(1) : dof.orders(0);
        for (auto &dof_sibling : (sibling_face->get_dofs())) {
          if (dof_cur != dof_sibling.cur_type)
            continue;
          unsigned int dof_sibling_match_order = dof_sibling.direction == u_dir
                                                     ? dof_sibling.orders(1)
                                                     : dof_sibling.orders(0);

          if (dof_match_order == dof_sibling_match_order) {
            // We also need dof_sibling_match to be active...
            auto fe_ptr_sibling = fe_collection_ptr->get_fe(
                sibling_face->get_fe_index(dof.cur_type));
            if (!fe_ptr_sibling->contains_this_dof(dof_sibling)) {
              dof.is_active = false;
              dof_sibling.is_active = false;
            } else if (dof.is_active) {
              // If this DoF is already active, because a cell was assigned
              // that DoF already, go through and determine the pairing
              fe_ptr->determine_pairing(
                  face.get_local_index(), sibling_face->get_local_index(),
                  face.get_master_orientation(), &(dof.coordinate_orientation),
                  &(dof.invert_value));
            } else {
              dof.is_active = true;
              dof_sibling.is_active = true;

              dof.sibling_global_index = dof_sibling.global_index;
              dof_sibling.sibling_global_index = dof.global_index;

              dof.active_index = this->n_active_dofs;
              dof_sibling.active_index = this->n_active_dofs;
              ++this->n_active_dofs;
            }
          }
        }
      }
    }
    // Get the parent index of this half edge and append these DoFs to it
    this->dof_parents[face.parent_cell_index()].attach_dof_face(&face);
  }

  // Now we must finalize all of the DoFParents...
  for (auto &dof_parent : this->dof_parents)
    dof_parent.finalize();
}

template <unsigned int dim, unsigned int spacedim, class CoefficientType,
          class Real>
const unsigned int
DoFHandler<dim, spacedim, CoefficientType, Real>::get_n_active_dofs_on_cell(
    unsigned int index) const {
  //  assert(index < dof_cells.size() && "Index out of range!");
  //  unsigned total_dofs = dof_cells[index].n_active_dofs();
  //  for (const auto& boundary_id : dof_cells[index].get_boundary_ids())
  //  {
  //    total_dofs += dof_faces[boundary_id].n_active_dofs();
  //  }
  //  return total_dofs;
  assert(index < dof_parents.size() && "Index out of range!");
  assert(dof_parents[index].is_finalized &&
         "The DofParent must be finalized before you can get the n_dofs!");
  return this->dof_parents[index].n_active_dofs();
}

template <unsigned int dim, unsigned int spacedim, class CoefficientType,
          class Real>
const unsigned int
DoFHandler<dim, spacedim, CoefficientType, Real>::get_n_dofs_on_cell(
    unsigned int index) const {
  //  assert(index < dof_cells.size() && "Index out of range!");
  //  unsigned total_dofs = dof_cells[index].n_dofs();
  //  for (const auto& boundary_id : dof_cells[index].get_boundary_ids())
  //  {
  //    total_dofs += dof_faces[boundary_id].n_dofs();
  //  }
  //  return total_dofs;
  assert(index < dof_parents.size() && "Index out of range!");
  // assert(dof_parents[index].is_finalized && "The DofParent must be finalized
  // before you can get the n_dofs!");

  return dof_parents[index].is_finalized ? this->dof_parents[index].n_dofs()
                                         : 0;
}
template <unsigned int dim, unsigned int spacedim, class CoefficientType,
          class Real>
const unsigned int
DoFHandler<dim, spacedim, CoefficientType, Real>::get_n_dofs() const {
  return n_dofs;
}
template <unsigned int dim, unsigned int spacedim, class CoefficientType,
          class Real>
const unsigned int
DoFHandler<dim, spacedim, CoefficientType, Real>::get_n_active_dofs() const {
  return n_active_dofs;
}
template <unsigned int dim, unsigned int spacedim, class CoefficientType,
          class Real>
const unsigned int
DoFHandler<dim, spacedim, CoefficientType, Real>::get_n_current_types() const {
  return this->fe_collection_collector->get_number_of_types();
}
template <unsigned int dim, unsigned int spacedim, class CoefficientType,
          class Real>
FECollectionCollector<dim, spacedim, Real> *
DoFHandler<dim, spacedim, CoefficientType,
           Real>::get_fe_collection_collector() {
  return fe_collection_collector;
}
template <unsigned int dim, unsigned int spacedim, class CoefficientType,
          class Real>
void DoFHandler<dim, spacedim, CoefficientType, Real>::
    spawn_p_reduced_dof_handler(
        DoFHandler<dim, spacedim, CoefficientType, Real> *target) {
  // In this function, we take from this existing DoFHandler, we copy over its
  // DoFCells and DoFFaces Then, after reducing the expansion order on the
  // cloned DoFs, we repopulate the DoF parents Since the the GalerkinSystem is
  // separate, we can then use this DoFHandler to immediately generate the
  // solution and run post-processing for a p-reduced model

  // First, copy over the correct pointer to the FE_Collection_Collector (which
  // may contain certain memoizations) and the mesh pointer
  target->replace_fe_collection_collector(this->fe_collection_collector);
  target->replace_mesh(this->mesh);
  // Next copy over the cells and the faces
  target->get_dof_cells() = this->get_dof_cells_copy();
  target->get_dof_faces() = this->get_dof_faces_copy();
  target->get_dof_parents() =
      std::vector(this->dof_parents.size(), DoFParent<dim, spacedim, Real>());
  target->n_dofs = this->n_dofs;
  // We are now ready to call the coarsen routine for the cells and faces
  for (const auto &cur_type :
       this->fe_collection_collector->get_current_types()) {
    for (auto &cell : target->get_dof_cells()) {
      cell.set_fe_index(
          this->fe_collection_collector->get_fe_collection(cur_type)->prev_fe(
              cell.get_fe_index(cur_type)),
          cur_type);
      //  cell.attach_active_fe_pointers(target->fe_collection_collector);
    }
    for (auto &face : target->get_dof_faces()) {
      face.set_fe_index(
          this->fe_collection_collector->get_fe_collection(cur_type)->prev_fe(
              face.get_fe_index(cur_type)),
          cur_type);
    }
  }
  // The new target DoFHandler is now ready for "distribute dofs"
}
template <unsigned int dim, unsigned int spacedim, class CoefficientType,
          class Real>
void DoFHandler<dim, spacedim, CoefficientType, Real>::
    replace_fe_collection_collector(FECollectionCollector<dim, spacedim, Real>
                                        *replacement_fe_collection_collector) {
  this->fe_collection_collector = replacement_fe_collection_collector;
}
template <unsigned int dim, unsigned int spacedim, class CoefficientType,
          class Real>
void DoFHandler<dim, spacedim, CoefficientType, Real>::replace_mesh(
    MeshBase<dim, spacedim> *mesh) {
  this->mesh = mesh;
}
template <unsigned int dim, unsigned int spacedim, class CoefficientType,
          class Real>
std::vector<CoefficientType> DoFHandler<dim, spacedim, CoefficientType, Real>::
    map_DoF_values_to_this_handler(
        DoFHandler<dim, spacedim, CoefficientType> *source_dof_handler,
        const std::vector<CoefficientType> &source_values) const {
  // First spawn a vector that matches the number of active unknowns on this
  // dofhandler
  std::vector<CoefficientType> target_values(this->n_active_dofs,
                                             CoefficientType(0));
  // Now, we loop through all the DoFs on the source_dof_handler, finding the
  // matched global dof indices With the matched Global DoF indices, we then
  // assign based on the active DoF indices Loop through all DoF parents on the
  // source_dof_handler
  for (auto &source_dof_parent : source_dof_handler->get_dof_parents()) {
    // Get the index of this DoF parent
    const unsigned int &dof_parent_index = source_dof_parent.index;
    const auto &target_dof_parent = this->get_dof_parent(dof_parent_index);
    // Loop through DoFs
    for (unsigned int source_dof_index = 0;
         source_dof_index < source_dof_parent.n_dofs(); ++source_dof_index) {
      const auto &source_dof = source_dof_parent.get_dof(source_dof_index);
      if (!source_dof.is_active)
        continue;
      else {
        // Get the current value
        const auto value_to_copy = source_values[source_dof.active_index];
        // Get global index
        const auto global_index = source_dof.global_index;
        const auto &target_dof = target_dof_parent.get_dof(source_dof_index);
        // Finally, map the value
        target_values[target_dof.active_index] = value_to_copy;
      }
    }
  }
  return target_values;
}

// template <unsigned int dim, unsigned int spacedim, class CoefficientType,
//           class Real>
// const unsigned int
// DoFHandler<dim, spacedim, CoefficientType, Real>::n_active_dofs() const {
//   return 0;
// }

DROMON_NAMESPACE_CLOSE
#endif // DROMON_DOFHANDLER_H
