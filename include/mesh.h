//
// Created by Jake J. Harmon (jake.harmon@ieee.org) on 7/16/21.
//

#ifndef DROMON_MESH_H
#define DROMON_MESH_H

#include "GeomBase.h"
#include "IteratorRanger.h"
#include "MeshBase.h"
#include "Point.h"
#include "config.h"
#include "utility.h"
#include <cmath>
#include <memory>
#include <vector>
#include <algorithm>

DROMON_NAMESPACE_OPEN

/**
 *
 * @tparam dim The dimension of the object
 * @tparam spacedim The dimension of the space that the object is embedded in
 * @tparam celltype The type of mapping, i.e. dim-linear, quadratic, etc.
 * The class which includes anything related to the geometry.
 * For a surface in 3-D, dim == 2, and spacedim == 3
 * For a wire in 3-D, then dim == 1, and spacedim == 3, and so on.
 * By default, spacedim is assumed to be 3.
 *
 * When a cell is created, it is of course granted a pointer
 * to the entire mesh that spawned it as well as the nodes
 * that form it.
 * The cell then inserts new edges to the overall mesh, as well as faces when
 * necessary. The cell is not "fed" the edge_indices and the face_indices is for
 * future support for multi-level hp refinement. When an edge/face exists
 * already, the index is simply stored.
 */

template <unsigned int dim, unsigned int spacedim, unsigned int celltype,
          class Real = double>
class Cell : public GeomBase<dim, dim, spacedim, double> {
private:
  // Each cell will have (2+celltype)^dim nodes
  // unsigned int node_indices [std::pow((int)(2+celltype), dim)];
  //    unsigned int node_indices [dim == 3 ?
  //    (2+celltype)*(2+celltype)*(2+celltype)
  //                                        : dim == 2 ?
  //                                        (2+celltype)*(2+celltype) :
  //                                        (2+celltype)];
  std::vector<unsigned int> node_indices;

  // Each cell will have 2*dim faces (nodes in 1-D, edges in 2-D, faces in 3-D),
  // noting, of course, that face is used here in the dimension-dependent sense,
  // namely that it is a $dim-1$ object that demarks the boundary of the cell.
  // unsigned int face_indices [2*dim];
  std::vector<int> face_indices;

public:
  Cell();
  Cell(Mesh<dim, spacedim, celltype> *mesh, unsigned int index,
       std::vector<unsigned int> node_indices);
  Cell(Mesh<dim, spacedim, celltype> *mesh, unsigned int index,
       std::array<unsigned int,
                  GeometryInfo<dim, spacedim, celltype>::nodes_per_cell>
           node_indices);

  std::vector<unsigned int> get_node_indices() const;
  unsigned int get_node_index(unsigned int i) const;
  virtual const unsigned int get_vertex_index(const unsigned int i) const override;

  virtual unsigned int type() const override;
  virtual unsigned int n_nodes() const override;
  virtual unsigned int n_vertices() const override;
  void push_back_face_index(unsigned int index) {
    face_indices.push_back(index);
    this->boundary_ids.push_back(index);
  }
  unsigned int n_tracked_faces() const { return face_indices.size(); }
  virtual Point<spacedim, Real> r(const Point<dim, Real> &uv) const override;
  virtual Point<spacedim, Real> drdu(const Point<dim, Real> &uv) const override;
  virtual Point<spacedim, Real> drdv(const Point<dim, Real> &uv) const override;
  virtual Point<spacedim, Real> drdw(const Point<dim, Real> &uv) const override;
  virtual Point<spacedim, Real> unitary_vector(const Point<dim, Real> &uv, const DoFDirection& dof_direction) const override;

  virtual Real jacobian(const Point<dim, Real> &uv) const override;
  virtual Real diameter() const override;
  //virtual const unsigned int get_face_index_on_neighbor(const unsigned int& local_index) const override;
  virtual void check_if_adjacent(GeomBase<dim, dim, spacedim, double>* neighbor, AdjacencyType* adjacencyType, MultiIndex<2, unsigned int>* local_indices_of_adjacency) const;
};

template <unsigned int dim, unsigned int spacedim, unsigned int celltype,
          class Real>
Cell<dim, spacedim, celltype, Real>::Cell(
    Mesh<dim, spacedim, celltype> *mesh, unsigned int index,
    std::vector<unsigned int> node_indices)
    : GeomBase<dim, dim, spacedim>(mesh),
      node_indices(dim == 3   ? (2 + celltype) * (2 + celltype) * (2 + celltype)
                   : dim == 2 ? (2 + celltype) * (2 + celltype)
                              : (2 + celltype))
//,face_indices(2*dim)
{
  assert(node_indices.size() == this->node_indices.size() &&
         "Mismatch on node_indices size!");
  this->index = index;
  for (unsigned int i = 0; i < this->node_indices.size(); ++i)
    this->node_indices[i] = node_indices[i];
}

template <unsigned int dim, unsigned int spacedim, unsigned int celltype,
          class Real>
Cell<dim, spacedim, celltype, Real>::Cell(
    Mesh<dim, spacedim, celltype> *mesh, unsigned int index,
    std::array<unsigned int,
               GeometryInfo<dim, spacedim, celltype>::nodes_per_cell>
        node_indices)
    : GeomBase<dim, dim, spacedim>(mesh),
      node_indices(GeometryInfo<dim, spacedim, celltype>::nodes_per_cell)
//,face_indices(2*dim)
{
  assert(node_indices.size() == this->node_indices.size() &&
         "Mismatch on node_indices size!");
  this->index = index;
  for (unsigned int i = 0; i < this->node_indices.size(); ++i)
    this->node_indices[i] = node_indices[i];
}

template <unsigned int dim, unsigned int spacedim, unsigned int celltype,
          class Real>
unsigned int Cell<dim, spacedim, celltype, Real>::type() const {
  return celltype;
}

template <unsigned int dim, unsigned int spacedim, unsigned int celltype,
          class Real>
std::vector<unsigned int>
Cell<dim, spacedim, celltype, Real>::get_node_indices() const {
  return std::vector<unsigned int>(node_indices.begin(), node_indices.end());
}

template <unsigned int dim, unsigned int spacedim, unsigned int celltype,
          class Real>
unsigned int Cell<dim, spacedim, celltype, Real>::n_nodes() const {
  return (this->node_indices).size();
}

template <unsigned int dim, unsigned int spacedim, unsigned int celltype,
          class Real>
Point<spacedim, Real>
Cell<dim, spacedim, celltype, Real>::r(const Point<dim, Real> &uv) const {
  Point<spacedim, Real> out;
  switch (dim) {
  case 2: {

    for (unsigned int m = 0; m <= celltype + 1; ++m)
    {
      const Real Lag_p_m = utility::lagrange_p(celltype+1, m, uv(0));
      for (unsigned int n = 0; n <= celltype+1; ++n)
      {
        // Jake: Note that we don't need to recompute the "m" terms (they could
        // be stored in a vector) In fact, both of these sets of lagrange terms
        // could be stored elsewhere and used across all cells.
        //out += this->mesh->get_node(node_indices[(celltype+2)*n + m])->get_point() * utility::lagrange_p(celltype+1, n, uv(1)) * Lag_p_m;
        const Real m_term = utility::lagrange_p(celltype+1, n, uv(1)) * Lag_p_m;
        out.add_multiple_of(this->mesh->get_node(node_indices[(celltype+2)*n + m])->get_point(), m_term);
      }
    }
    break;
  }
  default:
    assert(false && "Case not implemented!");
    break;
  }
  return out;
}

template <unsigned int dim, unsigned int spacedim, unsigned int celltype,
          class Real>
unsigned int
Cell<dim, spacedim, celltype, Real>::get_node_index(unsigned int i) const {
  assert(i < node_indices.size() && "Index out of range!");
  return node_indices[i];
}
template <unsigned int dim, unsigned int spacedim, unsigned int celltype,
          class Real>
Point<spacedim, Real>
Cell<dim, spacedim, celltype, Real>::drdu(const Point<dim, Real> &uv) const {
  Point<spacedim, Real> out;
  switch (dim) {
  case 2: {

    for (unsigned int m = 0; m <= celltype + 1; ++m)
    {
     const Real Lag_p_m = utility::lagrange_p_prime(celltype+1, m, uv(0));
      for (unsigned int n = 0; n <= celltype+1; ++n)
      {
        // Jake: Note that we don't need to recompute the "m" terms (they could
        // be stored in a vector) In fact, both of these sets of lagrange terms
        // could be stored elsewhere and used across all cells.
      //  out += this->mesh->get_node(node_indices[(celltype+2)*n + m])->get_point() * (utility::lagrange_p(celltype+1, n, uv(1)) * Lag_p_m);
        const Real m_term = utility::lagrange_p(celltype+1, n, uv(1)) * Lag_p_m;
        out.add_multiple_of(this->mesh->get_node(node_indices[(celltype+2)*n + m])->get_point(), m_term);
      }
    }
    break;
  }
  default:
    assert(false && "Case not implemented!");
    break;
  }
  return out;
}
template <unsigned int dim, unsigned int spacedim, unsigned int celltype,
          class Real>
Point<spacedim, Real>
Cell<dim, spacedim, celltype, Real>::drdv(const Point<dim, Real> &uv) const {
  Point<spacedim, Real> out;
  switch (dim) {
  case 2: {

    for (unsigned int m = 0; m <= celltype + 1; ++m)
    {
      const Real Lag_p_m = utility::lagrange_p(celltype+1, m, uv(0));
      for (unsigned int n = 0; n <= celltype+1; ++n)
      {
        // Jake: Note that we don't need to recompute the "m" terms (they could
        // be stored in a vector) In fact, both of these sets of lagrange terms
        // could be stored elsewhere and used across all cells.
       // out += this->mesh->get_node(node_indices[(celltype+2)*n + m])->get_point() * utility::lagrange_p_prime(celltype+1, n, uv(1)) * Lag_p_m;
        const Real m_term = utility::lagrange_p_prime(celltype+1, n, uv(1)) * Lag_p_m;
        out.add_multiple_of(this->mesh->get_node(node_indices[(celltype+2)*n + m])->get_point(), m_term);
      }
    }
    break;
  }
  default:
    assert(false && "Case not implemented!");
    break;
  }
  return out;
}
template <unsigned int dim, unsigned int spacedim, unsigned int celltype,
          class Real>
Real Cell<dim, spacedim, celltype, Real>::jacobian(
    const Point<dim, Real> &uv) const {
  if (dim == 2)
    return cross(drdu(uv), drdv(uv)).norm();
  else {
    assert(false && "Case not implemented!");
    return 0;
  }
}
template <unsigned int dim, unsigned int spacedim, unsigned int celltype,
          class Real>
Point<spacedim, Real>
Cell<dim, spacedim, celltype, Real>::drdw(const Point<dim, Real> &uv) const {
  assert(false && "Function not implemented!");
  return Point<spacedim, Real>();
}
template <unsigned int dim, unsigned int spacedim, unsigned int celltype,
          class Real>
const unsigned int
Cell<dim, spacedim, celltype, Real>::get_vertex_index(const unsigned int i) const {
  // assert(i < node_indices.size() && "Index out of range!");
  unsigned int node_index;
  switch (dim)
  {
  case 2:
    node_index = i > 1 ? (celltype+1)*(i + (celltype)) : i*(celltype+1);
    break;
  default:
    assert(false && "Not implemented for dim != 2!");
  }
  assert(node_index < this->n_nodes() && "Node_index exceeds maximum!");
  return node_indices[node_index];
}
template <unsigned int dim, unsigned int spacedim, unsigned int celltype,
          class Real>
unsigned int Cell<dim, spacedim, celltype, Real>::n_vertices() const {
  return GeometryInfo<dim, spacedim, celltype>::vertices_per_cell;
}
template <unsigned int dim, unsigned int spacedim, unsigned int celltype,
          class Real>
Point<spacedim, Real> Cell<dim, spacedim, celltype, Real>::unitary_vector(
    const Point<dim, Real> &uv, const DoFDirection &dof_direction) const {
  if (dof_direction == u_dir)
    return drdu(uv);
  else if (dof_direction == v_dir)
    return drdv(uv);
  else if (dof_direction == w_dir)
    return drdw(uv);
  else
  {
    assert(false && "Invalid DoFDirection!");
    return Point<spacedim, Real>();
  }
}
template <unsigned int dim, unsigned int spacedim, unsigned int celltype,
          class Real>
void Cell<dim, spacedim, celltype, Real>::check_if_adjacent(
    GeomBase<dim, dim, spacedim, double> *neighbor,
    AdjacencyType *adjacencyType,
    MultiIndex<2, unsigned int> *local_indices_of_adjacency) const {
  std::vector<unsigned int> matched_vertices_this, matched_vertices_neighbor;
 // Loop through the vertices of the two cells and track which ones match (if any)
  for (unsigned int i = 0; i < this->n_vertices(); ++i)
  {
    const unsigned int vertex_this = this->get_vertex_index(i);
    for (unsigned int j = 0; j < neighbor->n_vertices(); ++j)
    {
      const unsigned int vertex_neighbor = neighbor->get_vertex_index(j);
      if (vertex_this == vertex_neighbor)
      {
        matched_vertices_this.push_back(i);
        matched_vertices_neighbor.push_back(j);
      }
    }
  }
  if (matched_vertices_this.size() == 0)
  {
    *adjacencyType = Disjoint;
  }
  else if (matched_vertices_this.size() == 1)
  {
    *adjacencyType = Vertex;
    *local_indices_of_adjacency = {matched_vertices_this[0], matched_vertices_neighbor[0]};
  }else if (matched_vertices_this.size() == 2)
  {
    *adjacencyType = Edge;
    unsigned int local_index_this;
    unsigned int local_index_neighbor;
    // Need to figure out the local indices of the edges based on the local vertex indices...
    // To simplify conditionals, let's first sort these...
    std::sort(matched_vertices_this.begin(), matched_vertices_this.end());
    std::sort(matched_vertices_neighbor.begin(), matched_vertices_neighbor.end());
    if (matched_vertices_this[0] == 0)
    {
      if (matched_vertices_this[1] == 1)
        local_index_this = 0;
      else
        local_index_this = 2;
    }else if (matched_vertices_this[0] == 1)
      local_index_this = 1;
    else if (matched_vertices_this[0] == 2)
      local_index_this = 3;

    if (matched_vertices_neighbor[0] == 0)
    {
      if (matched_vertices_neighbor[1] == 1)
        local_index_neighbor = 0;
      else
        local_index_neighbor = 2;
    }else if (matched_vertices_neighbor[0] == 1)
      local_index_neighbor = 1;
    else if (matched_vertices_neighbor[0] == 2)
      local_index_neighbor = 3;

    *local_indices_of_adjacency = {local_index_this, local_index_neighbor};
  } else if (matched_vertices_this.size() == this->n_vertices())
    *adjacencyType = Self;
  else
    assert(false && "Not implemented yet, would be 3-D structures!");

}

/**
 * This function computes an estimate of the diameter.
 * Specifically, we take the maximum diagonal in the cell.
 * @tparam dim
 * @tparam spacedim
 * @tparam celltype
 * @tparam Real
 * @return
 */
template <unsigned int dim, unsigned int spacedim, unsigned int celltype,
          class Real>
Real Cell<dim, spacedim, celltype, Real>::diameter() const {
  // Depending on the dim, we must compute the estimate of the diameter differently
  switch (dim)
  {
  case 1: {
    const auto &p0 =
        this->mesh->get_node(this->get_vertex_index(0))->get_point();
    const auto &p1 =
        this->mesh->get_node(this->get_vertex_index(1))->get_point();
    return (p1-p0).norm();
  }
  case 2:
  {
    const auto& p0 = this->mesh->get_node(this->get_vertex_index(0))->get_point();
    const auto& p3 = this->mesh->get_node(this->get_vertex_index(3))->get_point();

    const auto& p1 = this->mesh->get_node(this->get_vertex_index(1))->get_point();
    const auto& p2 = this->mesh->get_node(this->get_vertex_index(2))->get_point();
    return std::max((p3-p0).norm(), (p2-p1).norm());
  }
  default:
    return 0;
  }
}


template <unsigned int dim, unsigned int spacedim>
/**
 *
 * @tparam dim
 * @tparam spacedim
 * Note that the way faces are defined on a cell are such that
 * In other words, the orientation matches the direction of the associated
 * unitary vector on that edge
 *           3
 *   2--------------->3
 *   ^                ^
 *   |                |
 * 2 |                | 1
 *   |                |
 *   |                |
 *   0--------------->1
 *           0
 */
class HalfFace : public GeomBase<dim - 1, dim, spacedim> {
public:
  std::vector<unsigned int> node_indices;

  unsigned int cell_index; // The index that the half face belongs to

  // The (signed) index of the sibling half face, where the sign keeps track of
  // different orientations
  unsigned int sibling_index;

  // If relative orientation is true, then this HalfFace has the master
  // orientation Otherwise, the orientation is not the master orientation and
  // does not match the master orientation
  bool master_orientation = true;

  // A half face has only one cell; however, with its sibling half face, it
  // forms a full face
  bool face_is_full = false;
  HalfFace();
  explicit HalfFace(MeshBase<dim, spacedim> *mesh,
                    const std::vector<unsigned int> &node_indices);
  void set_cell(unsigned int cell_index, unsigned int local_index);
  void set_sibling(unsigned int sibling_index, bool master_orientation = true);
  virtual unsigned int n_nodes() const override;
  virtual const unsigned int get_sibling_index() const override;
  virtual const bool get_orientation() const override;
  virtual void get_object_context(bool *has_sibling,
                                  unsigned int *sibling_index,
                                  bool *orientation) const override;
};

template <unsigned int dim, unsigned int spacedim>
HalfFace<dim, spacedim>::HalfFace(MeshBase<dim, spacedim> *mesh,
                                  const std::vector<unsigned int> &node_indices)
    : GeomBase<dim - 1, dim, spacedim>(mesh), node_indices(dim == 3   ? 4
                                                           : dim == 2 ? 2
                                                                      : 1) {
  assert(node_indices.size() == this->node_indices.size() &&
         "Input vector size of node_indices does not match!");
  for (unsigned int i = 0; i < this->node_indices.size(); ++i) {
    this->node_indices[i] = node_indices[i];
    this->boundary_ids.push_back(node_indices[i]);
  }
}

template <unsigned int dim, unsigned int spacedim>
HalfFace<dim, spacedim>::HalfFace() {}

template <unsigned int dim, unsigned int spacedim>
void HalfFace<dim, spacedim>::set_cell(unsigned int cell_index,
                                       unsigned int local_index) {
  this->cell_index = cell_index;
  this->local_index = local_index;

  // A bit of redundancy for the generalization to GeomBase...
  // Cell_index could be replaced by an alias, i.e., & cell_index =
  // parent_index...
  this->parent_index = cell_index;
}

template <unsigned int dim, unsigned int spacedim>
void HalfFace<dim, spacedim>::set_sibling(unsigned int sibling_index,
                                          bool master_orientation) {
  this->face_is_full = true;
  this->sibling_index = sibling_index;
  this->master_orientation = master_orientation;
}

template <unsigned int dim, unsigned int spacedim>
unsigned int HalfFace<dim, spacedim>::n_nodes() const {
  return (this->node_indices).size();
}
template <unsigned int dim, unsigned int spacedim>
const unsigned int HalfFace<dim, spacedim>::get_sibling_index() const {
  return sibling_index;
}
template <unsigned int dim, unsigned int spacedim>
void HalfFace<dim, spacedim>::get_object_context(bool *has_sibling,
                                                 unsigned int *sibling_index,
                                                 bool *orientation) const {
  *has_sibling = this->face_is_full;
  *sibling_index = this->sibling_index;
  *orientation = this->master_orientation;
}
template <unsigned int dim, unsigned int spacedim>
const bool HalfFace<dim, spacedim>::get_orientation() const {
  return this->master_orientation;
}

// template <unsigned int dim, unsigned int celltype = LINEARP, unsigned int
// spacedim = 3, class Real = double>
template <unsigned int dim, unsigned int spacedim, unsigned int celltype>
class Mesh : public MeshBase<dim, spacedim> {
public:
  static constexpr unsigned int cell_type = celltype;

  void push_back_face(HalfFace<dim, spacedim> *face_pb,
                      unsigned int host_cell_index);
  void push_back_cell(Cell<dim, spacedim, celltype> *cell_pb);

  unsigned int n_cells() const override;
  unsigned int n_faces() const override;
  // GeomBase<0, dim, spacedim>* get_node(unsigned int i);
  const unsigned int get_cell_type();

  void cells_from_grid(std::vector<unsigned int> &grid,
                       const unsigned int &grid_size);
  void insert_cell(
      const std::array<unsigned int,
                       GeometryInfo<dim, spacedim, celltype>::nodes_per_cell>
          &cell);

  void spawn_and_assign_faces();

  void activate_cell(const unsigned int cell_index);
  void deactivate_cell(const unsigned int cell_index);

  // Overide Iterator related functions...
  virtual typename MeshBase<dim, spacedim>::ActiveCellIterator
  end_cell() const override {
    return typename MeshBase<dim, spacedim>::ActiveCellIterator(
        const_cast<Mesh<dim, spacedim, celltype> *>(this),
        const_cast<Cell<dim, spacedim, celltype> *>(&cells[cells.size()]),
        n_cells());
  }

  virtual typename MeshBase<dim, spacedim>::ActiveCellIterator
  begin_active_cell() const override {
    for (unsigned int i = 0; i < this->cells.size(); ++i)
      if (cells[i].is_active())
        return typename MeshBase<dim, spacedim>::ActiveCellIterator(
            const_cast<Mesh<dim, spacedim, celltype> *>(this),
            const_cast<Cell<dim, spacedim, celltype> *>(&cells[i]), i);
    // return MeshBase<dim, spacedim>::ActiveCellIterator();
    return end_cell();
  }

  virtual IteratorRanger<typename MeshBase<dim, spacedim>::ActiveCellIterator>
  active_cell_iterators() const override {
    return IteratorRanger<typename MeshBase<dim, spacedim>::ActiveCellIterator>(
        begin_active_cell(), end_cell());
  }
  virtual GeomBase<dim, dim, spacedim> *
  get_pointer_to_cell(const unsigned int &index) const override;
  virtual GeomBase<dim - 1, dim, spacedim> *
  get_pointer_to_face(const unsigned int &index) const override;
  std::vector<Cell<dim, spacedim, celltype>> cells;
  std::vector<HalfFace<dim, spacedim>>
      half_faces; // The boundary of the element, nodes in 1-D, edges in 2-D,
                  // faces in 3-D
private:
  // Since we relying on the GeomBase here, these must be stored as pointers to
  // prevent object slicing...
  // std::vector<std::unique_ptr<Face<dim, spacedim>>> faces; //The boundary of
  // the element, nodes in 1-D, edges in 2-D, faces in 3-D

  // The entire element, a line in 1-D, a face in 2-D, a volume in 3-D
  //  std::vector<std::unique_ptr<Cell<dim, spacedim, celltype>>> cells;
};

DROMON_NAMESPACE_CLOSE
#endif // DROMON_MESH_H
