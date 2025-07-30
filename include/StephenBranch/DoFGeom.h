//
// Created by Jake J. Harmon (jake.harmon@ieee.org) on 10/5/21.
//

#ifndef DROMON_DOFGEOM_H
#define DROMON_DOFGEOM_H

#include "DoFBase.h"
#include "DoFGeomBase.h"
#include "FEBase.h"
#include "FECollection.h"
#include "GeomBase.h"
#include "config.h"

DROMON_NAMESPACE_OPEN

template <unsigned int dim, unsigned int spacedim, class Real = double>
class DoFCell : public DoFGeomBase<dim> {
public:
  using CellBase = GeomBase<dim, dim, spacedim>;

  static constexpr unsigned int cell_dim = dim;
  DoFCell(unsigned int index, bool active, unsigned int geom_level,
          CellBase *cb_pointer);
  virtual void clear_dofs() override;
  virtual void activate_dofs() override;
  virtual void deactivate_dofs() override;
  virtual unsigned int n_dofs() const override;
  virtual unsigned int n_active_dofs() const override;
  void push_back(const DoFBase<dim> &dof);
  // TODO: Generalize for arbitrary dim...
  virtual const bool is_compatible(const DoFDirection &dir,
                                   const unsigned int &i,
                                   const unsigned int &j) const override {
    return true;
  }

  CellBase *get_cell_on_mesh() const;
  const std::vector<unsigned int> get_boundary_ids() const;
  std::vector<DoFBase<dim>> &get_dofs();
  const DoFBase<dim> &get_dof(const unsigned int &i) const;

  friend bool operator==(const DoFCell<dim, spacedim> &m1,
                         const DoFCell<dim, spacedim> &m2) {
    return (m1.cell == m2.cell && m1.n_dofs() == m2.n_dofs() &&
            m1.n_active_dofs() == m2.n_active_dofs());
  }
  friend bool operator!=(const DoFCell<dim, spacedim> &m1,
                         const DoFCell<dim, spacedim> &m2) {
    return !(m1 == m2);
  }

  void attach_active_fe_pointers(
      FECollectionCollector<dim, spacedim, Real> *feCollectionCollector);
  FiniteElement<dim, Real> *
  get_active_fe_pointer(const unsigned int &dof_index) const;
  FiniteElement<dim, Real> *
  get_active_fe_pointer(const CurrentType &cur_type) const;

private:
  std::vector<DoFBase<dim>> dofs;

  std::vector<FiniteElement<dim, Real> *> active_fe_pointers;

  // A pointer to the geometric object this DoFCell is associated with
  CellBase *cell;
};

template <unsigned int dim, unsigned int spacedim, class Real>
DoFCell<dim, spacedim, Real>::DoFCell(unsigned int index, bool active,
                                      unsigned int geom_level,
                                      CellBase *cb_pointer)
    : DoFGeomBase<dim>(index, active, geom_level, DoFGeom::DoFGeomType::Cell) {
  this->cell = cb_pointer;
}

template <unsigned int dim, unsigned int spacedim, class Real>
void DoFCell<dim, spacedim, Real>::clear_dofs() {
  dofs.clear();
}

template <unsigned int dim, unsigned int spacedim, class Real>
void DoFCell<dim, spacedim, Real>::activate_dofs() {
  this->active = true;
  for (auto &dof : dofs) {
    dof.is_active = true;
  }
}

template <unsigned int dim, unsigned int spacedim, class Real>
void DoFCell<dim, spacedim, Real>::deactivate_dofs() {
  this->active = false;
  for (auto &dof : dofs) {
    dof.is_active = false;
    dof.active_index = std::numeric_limits<int>::max();
  }
}

template <unsigned int dim, unsigned int spacedim, class Real>
void DoFCell<dim, spacedim, Real>::push_back(const DoFBase<dim> &dof) {
  dofs.push_back(dof);
}

template <unsigned int dim, unsigned int spacedim, class Real>
unsigned int DoFCell<dim, spacedim, Real>::n_dofs() const {
  return dofs.size();
}
template <unsigned int dim, unsigned int spacedim, class Real>
std::vector<DoFBase<dim>> &DoFCell<dim, spacedim, Real>::get_dofs() {
  return dofs;
}
template <unsigned int dim, unsigned int spacedim, class Real>
const std::vector<unsigned int>
DoFCell<dim, spacedim, Real>::get_boundary_ids() const {
  return cell->boundary_ids;
}
template <unsigned int dim, unsigned int spacedim, class Real>
unsigned int DoFCell<dim, spacedim, Real>::n_active_dofs() const {
  unsigned int active_dofs = 0;
  for (const auto &dof : dofs) {
    if (dof.is_active)
      ++active_dofs;
  }
  return active_dofs;
}
template <unsigned int dim, unsigned int spacedim, class Real>
typename DoFCell<dim, spacedim, Real>::CellBase *
DoFCell<dim, spacedim, Real>::get_cell_on_mesh() const {
  return cell;
}
template <unsigned int dim, unsigned int spacedim, class Real>
FiniteElement<dim, Real> *DoFCell<dim, spacedim, Real>::get_active_fe_pointer(
    const unsigned int &dof_index) const {
  return this->active_fe_pointers[this->get_relevant_fe_index(
      dofs[dof_index].cur_type)];
}
template <unsigned int dim, unsigned int spacedim, class Real>
const DoFBase<dim> &
DoFCell<dim, spacedim, Real>::get_dof(const unsigned int &i) const {
  #ifdef DEBUG
  assert(i < dofs.size() && "Index out of range!");
  #endif
  return dofs[i];
}
template <unsigned int dim, unsigned int spacedim, class Real>
FiniteElement<dim, Real> *DoFCell<dim, spacedim, Real>::get_active_fe_pointer(
    const CurrentType &cur_type) const {
  return this->active_fe_pointers[this->get_relevant_fe_index(cur_type)];
}
template <unsigned int dim, unsigned int spacedim, class Real>
void DoFCell<dim, spacedim, Real>::attach_active_fe_pointers(
    FECollectionCollector<dim, spacedim, Real> *feCollectionCollector) {
  active_fe_pointers.resize(this->fe_indices.size());
  this->active_fe_degrees.resize(this->fe_indices.size());

  // Add the necessary pointers and then set the active degrees...
  // For the current cell, get the indices in the FECollectionCollector
  for (unsigned int i = 0; i < this->fe_indices.size(); ++i) {
    auto feCollection = feCollectionCollector->get_fe_collection(
        i); // pointer to the FECollection
     active_fe_pointers[i] = feCollection->get_fe(this->fe_indices[i]);
    this->active_fe_degrees[i] = active_fe_pointers[i]->get_degree();
  }
}

template <unsigned int dim, unsigned int spacedim, class Real = double>
class DoFFace : public DoFGeomBase<dim> {
public:
  using FaceBase = GeomBase<dim - 1, dim, spacedim>;

  DoFFace(unsigned int index, bool active, unsigned int geom_level,
          FaceBase *fb_pointer);
  virtual void clear_dofs() override;
  virtual void activate_dofs() override;
  virtual void deactivate_dofs() override;
  virtual unsigned int n_dofs() const override;
  virtual unsigned int n_active_dofs() const override;

  // TODO: Generalize for arbitrary dim...
  virtual const bool is_compatible(const DoFDirection &dir,
                                   const unsigned int &i,
                                   const unsigned int &j) const override;
  void push_back(const DoFBase<dim> &dof);
  std::vector<DoFBase<dim>> &get_dofs();
  virtual const unsigned int get_neighbor_index() const override;
  const unsigned int get_local_index() const;
  virtual const unsigned int parent_cell_index() const override;
  const bool get_master_orientation() const;
  void get_object_context(bool *has_sibling, unsigned int *sibling_index,
                          bool *orientation) const;

private:
  std::vector<DoFBase<dim>> dofs;
  // A pointer to the geometric object this DoFFace is associated with
  FaceBase *face;
};

template <unsigned int dim, unsigned int spacedim, class Real>
DoFFace<dim, spacedim, Real>::DoFFace(unsigned int index, bool active,
                                      unsigned int geom_level,
                                      FaceBase *fb_pointer)
    : DoFGeomBase<dim>(index, active, geom_level, DoFGeom::DoFGeomType::Face) {
  this->face = fb_pointer;
}

template <unsigned int dim, unsigned int spacedim, class Real>
void DoFFace<dim, spacedim, Real>::clear_dofs() {
  dofs.clear();
}

template <unsigned int dim, unsigned int spacedim, class Real>
void DoFFace<dim, spacedim, Real>::activate_dofs() {
  this->active = true;
  for (auto &dof : dofs) {
    dof.is_active = true;
  }
}

template <unsigned int dim, unsigned int spacedim, class Real>
void DoFFace<dim, spacedim, Real>::deactivate_dofs() {
  this->active = false;
  for (auto &dof : dofs) {
    dof.is_active = false;
    dof.active_index = std::numeric_limits<int>::max();
  }
}

template <unsigned int dim, unsigned int spacedim, class Real>
void DoFFace<dim, spacedim, Real>::push_back(const DoFBase<dim> &dof) {
  dofs.push_back(dof);
}

template <unsigned int dim, unsigned int spacedim, class Real>
unsigned int DoFFace<dim, spacedim, Real>::n_dofs() const {
  return dofs.size();
}

template <unsigned int dim, unsigned int spacedim, class Real>
const bool
DoFFace<dim, spacedim, Real>::is_compatible(const DoFDirection &dir,
                                            const unsigned int &i,
                                            const unsigned int &j) const {
  // Check the Face object on the mesh to see if the expansion orders chosen are
  // acceptable
  auto local_half_edge_index = face->local_index;
  switch (dim) {
  case 2:
    // In 2-D, the local index of the half edge must either be 3 or 1 (for i = 0
    // and i = 1, respectively)
    switch (dir) {
    case u_dir:
      if ((i == 0 && local_half_edge_index == 2) ||
          (i == 1 && local_half_edge_index == 1))
        return true;
      else
        return false;
      break;
    case v_dir:
      if ((j == 0 && local_half_edge_index == 0) ||
          (j == 1 && local_half_edge_index == 3))
        return true;
      else
        return false;
      break;
    default:
      assert(false && "Direction not supported for this dim!");
    }
    break;
  default:
    assert(false && "Case not supported!");
  }
  return false;
}
template <unsigned int dim, unsigned int spacedim, class Real>
std::vector<DoFBase<dim>> &DoFFace<dim, spacedim, Real>::get_dofs() {
  return dofs;
}

template <unsigned int dim, unsigned int spacedim, class Real>
const unsigned int DoFFace<dim, spacedim, Real>::get_neighbor_index() const {
  return face->get_sibling_index();
}
template <unsigned int dim, unsigned int spacedim, class Real>
const unsigned int DoFFace<dim, spacedim, Real>::get_local_index() const {
  return face->local_index;
}
template <unsigned int dim, unsigned int spacedim, class Real>
const bool DoFFace<dim, spacedim, Real>::get_master_orientation() const {
  return face->get_orientation();
}
template <unsigned int dim, unsigned int spacedim, class Real>
void DoFFace<dim, spacedim, Real>::get_object_context(
    bool *has_sibling, unsigned int *sibling_index, bool *orientation) const {
  face->get_object_context(has_sibling, sibling_index, orientation);
}
template <unsigned int dim, unsigned int spacedim, class Real>
unsigned int DoFFace<dim, spacedim, Real>::n_active_dofs() const {
  unsigned int active_dofs = 0;
  for (const auto &dof : dofs) {
    if (dof.is_active)
      ++active_dofs;
  }
  return active_dofs;
}
template <unsigned int dim, unsigned int spacedim, class Real>
const unsigned int DoFFace<dim, spacedim, Real>::parent_cell_index() const {
  return face->parent_index;
}

/**
 * The DoFParent unifies the functionality of the DoFCell and DoFFace
 * It contains a pointer to the DoFCell and a vector of pointers to the DoFFace
 * along with a single vector of pointers to the DoFs along with some
 * functionality borrowed from the DoFCell...
 * @tparam dim
 * @tparam spacedim
 * @tparam Real
 */
template <unsigned int dim, unsigned int spacedim, class Real = double>
class DoFParent {
public:
  static constexpr unsigned int cell_dim = dim;
  static constexpr unsigned int space_dim = spacedim;

  using CellBase = GeomBase<dim, dim, spacedim>;
  using FaceBase = GeomBase<dim - 1, dim, spacedim>;

  DoFParent() = default;

  // Attaches the pointer to the DoFCell and appends all its dofs
  void attach_dof_cell(DoFCell<dim, spacedim, Real> *dof_cell);
  const unsigned int get_current_fe_index(const CurrentType& cur_type) const;
  void set_fe_index(const unsigned int &new_fe_index, const CurrentType& cur_type);
  // Attaches the pointer to the DoFFace and appends all its dofs
  void attach_dof_face(DoFFace<dim, spacedim, Real> *dof_face);

  void finalize();

  void reset();

  const DoFBase<dim> &get_dof(const unsigned int &i) const;
  const unsigned int &get_index() const;

  // The following functions access the same functionality on the DoFCell
  unsigned int n_dofs() const;
  unsigned int n_active_dofs() const;
  const unsigned int active_degree(const unsigned int &cur_type) const;
  CellBase *get_cell_on_mesh() const;
  FiniteElement<dim, Real> *
  get_active_fe_pointer(const unsigned int &dof_index) const;

  bool is_finalized = false;
  unsigned int index;

  friend bool operator==(const DoFParent<dim, spacedim, Real> &m1,
                         const DoFParent<dim, spacedim, Real> &m2) {
    return *(m1.cell) == *(m2.cell);
  }
  friend bool operator!=(const DoFParent<dim, spacedim, Real> &m1,
                         const DoFParent<dim, spacedim, Real> &m2) {
    return !(m1 == m2);
  }

private:
  DoFCell<dim, spacedim, Real> *cell = nullptr;
  std::vector<DoFFace<dim, spacedim, Real> *> faces;
  std::vector<DoFBase<dim> *> dofs;

  // If none of the DoFs are active, then is_active is false; otherwise,
  // it is true
  bool is_active = false;

  unsigned int active_dofs = 0;
  void assign_active_fe_pointers_to_dofs();
};

template <unsigned int dim, unsigned int spacedim, class Real>
void DoFParent<dim, spacedim, Real>::attach_dof_cell(
    DoFCell<dim, spacedim, Real> *dof_cell) {
  if (cell == nullptr)
    this->cell = dof_cell;

  // Go through each DoF on the cell and pushback their pointers...
  for (auto &dof : this->cell->get_dofs()) {
    //if (dof.is_new)
      dofs.push_back(&dof);
   // dof.is_new = false;
  }
  if (this->cell->active)
    this->is_active = true;

  this->index = this->cell->index;
}
template <unsigned int dim, unsigned int spacedim, class Real>
void DoFParent<dim, spacedim, Real>::attach_dof_face(
    DoFFace<dim, spacedim, Real> *dof_face) {
  if (std::find(faces.begin(), faces.end(), dof_face) == faces.end())
    this->faces.push_back(dof_face);

  // Go through each DoF on the face and pushback their pointers
  for (auto &dof : dof_face->get_dofs()) {
   // if (dof.is_new)
      dofs.push_back(&dof);
  //  dof.is_new = false;
  }
  // dofs.insert(std::end(dofs), std::begin(dof_face->get_dofs()),
  // std::end(dof_face->get_dofs()));

  if (dof_face->active)
    this->is_active = true;
}

template <unsigned int dim, unsigned int spacedim, class Real>
void DoFParent<dim, spacedim, Real>::reset() {
  this->cell = nullptr;
  this->faces.clear();
  this->dofs.clear();
  this->is_active = false;
  active_dofs = 0;
}
template <unsigned int dim, unsigned int spacedim, class Real>
unsigned int DoFParent<dim, spacedim, Real>::n_dofs() const {
  return dofs.size();
}
template <unsigned int dim, unsigned int spacedim, class Real>
unsigned int DoFParent<dim, spacedim, Real>::n_active_dofs() const {
  return active_dofs;
}
template <unsigned int dim, unsigned int spacedim, class Real>
void DoFParent<dim, spacedim, Real>::finalize() {
  active_dofs = 0;
  for (const auto &dof : dofs)
    if (dof->is_active)
      ++active_dofs;
  is_finalized = true;
  this->assign_active_fe_pointers_to_dofs();
}
template <unsigned int dim, unsigned int spacedim, class Real>
const unsigned int DoFParent<dim, spacedim, Real>::active_degree(
    const unsigned int &cur_type) const {
  return cell->active_degree(cur_type);
}
template <unsigned int dim, unsigned int spacedim, class Real>
typename DoFParent<dim, spacedim, Real>::CellBase *
DoFParent<dim, spacedim, Real>::get_cell_on_mesh() const {
  return cell->get_cell_on_mesh();
}
template <unsigned int dim, unsigned int spacedim, class Real>
const DoFBase<dim> &
DoFParent<dim, spacedim, Real>::get_dof(const unsigned int &i) const {
  assert(i < dofs.size() && "Index i out of range!");

  return *(this->dofs[i]);
}
template <unsigned int dim, unsigned int spacedim, class Real>
FiniteElement<dim, Real> *DoFParent<dim, spacedim, Real>::get_active_fe_pointer(
    const unsigned int &dof_index) const {
  return this->cell->get_active_fe_pointer(dofs[dof_index]->cur_type);
}
template <unsigned int dim, unsigned int spacedim, class Real>
void DoFParent<dim, spacedim, Real>::assign_active_fe_pointers_to_dofs() {
  for (auto &dof : dofs)
    dof->set_active_fe_pointer(
        this->cell->get_active_fe_pointer(dof->cur_type));
}
template <unsigned int dim, unsigned int spacedim, class Real>
void DoFParent<dim, spacedim, Real>::set_fe_index(const unsigned int &new_fe_index, const CurrentType& cur_type) {
  // Set the fe_indices for the cell and faces
  cell->set_fe_index(new_fe_index, cur_type);

  for (auto face : faces)
    face->set_fe_index(new_fe_index, cur_type);
}
template <unsigned int dim, unsigned int spacedim, class Real>
const unsigned int DoFParent<dim, spacedim, Real>::get_current_fe_index(const CurrentType& cur_type) const {
  return this->cell->get_fe_index(cur_type);
}
template <unsigned int dim, unsigned int spacedim, class Real>
const unsigned int &DoFParent<dim, spacedim, Real>::get_index() const {
  return this->cell->get_index();
}

DROMON_NAMESPACE_CLOSE
#endif // DROMON_DOFGEOM_H
