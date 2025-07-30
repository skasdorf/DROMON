//
// Created by Jake J. Harmon (jake.harmon@ieee.org) on 7/17/21.
//

#ifndef DROMON_GEOMBASE_H
#define DROMON_GEOMBASE_H
#include "config.h"
#include "Point.h"
#include "MultiIndex.h"

DROMON_NAMESPACE_OPEN

//Forward declarations
template <unsigned int dim, unsigned int spacedim>
class MeshBase;

template <unsigned int dim, unsigned int spacedim, unsigned int celltype>
class Mesh;


template <unsigned int structdim, unsigned int dim, unsigned int spacedim, class Real = double>
class GeomBase
{
public:
  static constexpr unsigned int geom_spacedim = spacedim;

  GeomBase(MeshBase<dim, spacedim>* mesh,const unsigned int& level = 0, bool is_active = true);

  virtual unsigned int level() const;
  virtual bool is_active() const;
  virtual unsigned int type() const;
  virtual unsigned int n_nodes() const;
  virtual unsigned int n_vertices() const;
  virtual void deactivate();
  virtual void activate();

  virtual const unsigned int get_sibling_index() const;
  virtual const unsigned int get_face_index_on_neighbor(const unsigned int&) const;
  virtual const unsigned int get_vertex_index(const unsigned int) const;
  virtual void get_object_context(bool* has_sibling, unsigned int* sibling_index, bool* orientation) const;
  virtual const bool get_orientation() const;
  virtual void check_if_adjacent(GeomBase<structdim, dim, spacedim, Real>*, AdjacencyType*, MultiIndex<2, unsigned int>*) const;






  virtual Point<spacedim, Real> r(const Point<structdim, Real>& uv) const;
  virtual Point<spacedim, Real> drdu(const Point<structdim, Real>& uv) const;
  virtual Point<spacedim, Real> drdv(const Point<structdim, Real>& uv) const;
  virtual Point<spacedim, Real> drdw(const Point<structdim, Real>& uv) const;
  virtual unsigned int material_domain_id() const;
  virtual Point<spacedim, Real> unitary_vector(const Point<structdim, Real>& uv, const DoFDirection& dof_direction) const;
  virtual Real jacobian(const Point<dim, Real>& uv) const;
  virtual Real diameter() const;

  unsigned int index;
  unsigned int local_index; //A local index for the geometrical object
  unsigned int parent_index = std::numeric_limits<unsigned int>::max();

  std::vector<unsigned int> boundary_ids;

protected:

  const unsigned int geom_level;
  bool active;
  MeshBase<dim, spacedim>* mesh;
};




DROMON_NAMESPACE_CLOSE
#endif //DROMON_GEOMBASE_H
