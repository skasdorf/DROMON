//
// Created by Jake J. Harmon (jake.harmon@ieee.org) on 7/17/21.

#include "GeomBase.h"
#include "MeshBase.h"

template<unsigned int structdim, unsigned int dim, unsigned int spacedim, class Real>
dromon::GeomBase<structdim, dim, spacedim, Real>::GeomBase(MeshBase<dim, spacedim> *mesh, const unsigned int &level, bool is_active) : geom_level(level), active(is_active)//, boundary_ids(dim == 3 ? 4 : dim == 2 ? 2 : 1)
{
  this->mesh = mesh;
}

template<unsigned int structdim, unsigned int dim, unsigned int spacedim, class Real>
unsigned int dromon::GeomBase<structdim, dim, spacedim, Real>::type() const {
  return 0;
}

template<unsigned int structdim, unsigned int dim, unsigned int spacedim, class Real>
unsigned int dromon::GeomBase<structdim, dim, spacedim, Real>::n_nodes() const {
  return 1;
}

template<unsigned int structdim, unsigned int dim, unsigned int spacedim, class Real>
unsigned int dromon::GeomBase<structdim, dim, spacedim, Real>::level() const {
  return this->geom_level;
}

template<unsigned int structdim, unsigned int dim, unsigned int spacedim, class Real>
bool dromon::GeomBase<structdim, dim, spacedim, Real>::is_active() const {
  return this->active;
}

template<unsigned int structdim, unsigned int dim, unsigned int spacedim, class Real>
void dromon::GeomBase<structdim, dim, spacedim, Real>::deactivate() {
  this->active = false;
}

template<unsigned int structdim, unsigned int dim, unsigned int spacedim, class Real>
void dromon::GeomBase<structdim, dim, spacedim, Real>::activate() {
  this->active = true;
}

template<unsigned int structdim, unsigned int dim, unsigned int spacedim, class Real>
Point<spacedim, Real> dromon::GeomBase<structdim, dim, spacedim, Real>::r(const Point<structdim, Real> &uv) const {
  return Point<spacedim, Real>();
}

template<unsigned int structdim, unsigned int dim, unsigned int spacedim, class Real>
const unsigned int dromon::GeomBase<structdim, dim, spacedim, Real>::get_sibling_index() const {
  assert(false && "Exception: function not implemented for this class!");
  return 0;
}
template <unsigned int structdim, unsigned int dim, unsigned int spacedim,
          class Real>
void dromon::GeomBase<structdim, dim, spacedim, Real>::get_object_context(
    bool* has_sibling, unsigned int *sibling_index, bool *orientation) const
{
  assert(false && "Exception: function not implemented for this class!");
}
template<unsigned int structdim, unsigned int dim, unsigned int spacedim, class Real>
const bool dromon::GeomBase<structdim, dim, spacedim, Real>::get_orientation() const {
  assert(false && "Exception: function not implemented for this class!");
  return 0;
}
template <unsigned int structdim, unsigned int dim, unsigned int spacedim,
          class Real>
Point<spacedim, Real> dromon::GeomBase<structdim, dim, spacedim, Real>::drdu(
    const Point<structdim, Real> &uv) const {
  return Point<spacedim, Real>();
}
template <unsigned int structdim, unsigned int dim, unsigned int spacedim,
          class Real>
Point<spacedim, Real> dromon::GeomBase<structdim, dim, spacedim, Real>::drdv(
    const Point<structdim, Real> &uv) const {
  return Point<spacedim, Real>();
}
template <unsigned int structdim, unsigned int dim, unsigned int spacedim,
          class Real>
Point<spacedim, Real> dromon::GeomBase<structdim, dim, spacedim, Real>::drdw(
    const Point<structdim, Real> &uv) const {
  assert(dim == 3 && "Not valid for dim != 3");
  return Point<spacedim, Real>();
}
template <unsigned int structdim, unsigned int dim, unsigned int spacedim,
          class Real>
Real dromon::GeomBase<structdim, dim, spacedim, Real>::jacobian(
    const Point<dim, Real> &uv) const {
  return 0;
}
template <unsigned int structdim, unsigned int dim, unsigned int spacedim,
          class Real>
const unsigned int
dromon::GeomBase<structdim, dim, spacedim, Real>::get_vertex_index(const unsigned int) const {
  assert(false && "Function not implemented by this base class!");
  return 0;
}
template <unsigned int structdim, unsigned int dim, unsigned int spacedim,
          class Real>
unsigned int
dromon::GeomBase<structdim, dim, spacedim, Real>::n_vertices() const {
  return 0;
}
template <unsigned int structdim, unsigned int dim, unsigned int spacedim,
          class Real>
Point<spacedim, Real>
dromon::GeomBase<structdim, dim, spacedim, Real>::unitary_vector(
    const Point<structdim, Real> &uv, const DoFDirection& dof_direction) const {
  return Point<spacedim, Real>();
}
template <unsigned int structdim, unsigned int dim, unsigned int spacedim,
          class Real>
unsigned int
dromon::GeomBase<structdim, dim, spacedim, Real>::material_domain_id() const {
  return mesh->material_domain_id;
}
template <unsigned int structdim, unsigned int dim, unsigned int spacedim,
          class Real>
const unsigned int
dromon::GeomBase<structdim, dim, spacedim, Real>::get_face_index_on_neighbor(
    const unsigned int&) const {
  assert (false && "Not implemented for the base class!");
  return 0;
}
template <unsigned int structdim, unsigned int dim, unsigned int spacedim,
          class Real>
void
dromon::GeomBase<structdim, dim, spacedim, Real>::check_if_adjacent(
    dromon::GeomBase<structdim, dim, spacedim, Real> *,
    dromon::AdjacencyType *, dromon::MultiIndex<2, unsigned int> *) const {
  return;
}
template <unsigned int structdim, unsigned int dim, unsigned int spacedim,
          class Real>
Real dromon::GeomBase<structdim, dim, spacedim, Real>::diameter() const {
  return 0;
}

// -------------------------------------------------------------------
// -- Explicit specializations

//In the event of needing to use other floating-point representations, e.g., float or quad,
//add the specializations as below...
template class dromon::GeomBase<0,1,1, double>;
template class dromon::GeomBase<0,1,2, double>;
template class dromon::GeomBase<0,1,3, double>;
template class dromon::GeomBase<0,2,1, double>;
template class dromon::GeomBase<0,2,2, double>;
template class dromon::GeomBase<0,2,3, double>;
template class dromon::GeomBase<0,3,1, double>;
template class dromon::GeomBase<0,3,2, double>;
template class dromon::GeomBase<0,3,3, double>;

template class dromon::GeomBase<1,1,1, double>;
template class dromon::GeomBase<1,1,2, double>;
template class dromon::GeomBase<1,1,3, double>;
template class dromon::GeomBase<1,2,1, double>;
template class dromon::GeomBase<1,2,2, double>;
template class dromon::GeomBase<1,2,3, double>;
template class dromon::GeomBase<1,3,1, double>;
template class dromon::GeomBase<1,3,2, double>;
template class dromon::GeomBase<1,3,3, double>;

template class dromon::GeomBase<2,2,2, double>;
template class dromon::GeomBase<2,2,3, double>;
template class dromon::GeomBase<2,3,3, double>;

template class dromon::GeomBase<3,3,3, double>;

