//
// Created by Jake J. Harmon (jake.harmon@ieee.org) on 7/16/21.
//
#include "mesh.h"

template<unsigned int dim, unsigned int spacedim, unsigned int celltype>
void dromon::Mesh<dim, spacedim, celltype>::push_back_cell(Cell<dim, spacedim, celltype> *cell_pb)
{
  this->cells.push_back(*cell_pb);
}

template<unsigned int dim, unsigned int spacedim, unsigned int celltype>
const unsigned int dromon::Mesh<dim, spacedim, celltype>::get_cell_type() {
  return celltype;
}

template<unsigned int dim, unsigned int spacedim, unsigned int celltype>
void dromon::Mesh<dim, spacedim, celltype>::push_back_face(dromon::HalfFace<dim, spacedim> *face_pb, unsigned int host_cell_index) {
  //check if this face already exists on the same level...
  int found_flag = 0;
  int sibling_half_edge_index;

  for (unsigned int i = 0; i < half_faces.size(); ++i)
  {
    if (half_faces[i].level() != face_pb->level() || half_faces[i].n_nodes() != face_pb->n_nodes())
      continue;
    else{
      //loop through the vertices
      unsigned int start = 0;

      unsigned int end = half_faces[i].n_nodes() - 1;
      int matched_vertices = 0;
      while (start < end)
      {
        if (half_faces[i].node_indices[start] == face_pb->node_indices[start] && half_faces[i].node_indices[end] == face_pb->node_indices[end])
        {
          matched_vertices +=2;
        }else if (half_faces[i].node_indices[end] == face_pb->node_indices[start] && half_faces[i].node_indices[start] == face_pb->node_indices[end])
        {
          matched_vertices -= 2;
        }
        else
          break;

        ++start;
        --end;
      }
      if (matched_vertices == face_pb->n_nodes())
      {
        found_flag = 1;
        sibling_half_edge_index = i;

        break;
      }else if (matched_vertices == -face_pb->n_nodes())
      {
        found_flag = -1;
        sibling_half_edge_index = i;

        break;
      }
    }
  }
  this->half_faces.push_back(*face_pb);
  half_faces[half_faces.size()-1].set_cell(host_cell_index, this->cells[host_cell_index].n_tracked_faces());
  this->cells[host_cell_index].push_back_face_index(half_faces.size()-1);

  if (found_flag == 1)
  {
    half_faces[half_faces.size()-1].set_sibling(sibling_half_edge_index, true);
    half_faces[sibling_half_edge_index].set_sibling(half_faces.size()-1, true);
  }else if (found_flag == -1){
    half_faces[half_faces.size()-1].set_sibling(sibling_half_edge_index, false);
    half_faces[sibling_half_edge_index].set_sibling(half_faces.size()-1, true);
  }
}

template<unsigned int dim, unsigned int spacedim, unsigned int celltype>
void dromon::Mesh<dim, spacedim, celltype>::cells_from_grid(std::vector<unsigned int> &grid, const unsigned int& grid_size) {
  unsigned cells_along_dim = (grid_size-1)/(celltype + 1);
  switch (dim)
  {
  case 1:
    break;
  case 2:
    for (unsigned int i = 0; i < cells_along_dim; ++i)
      for (unsigned int j = 0; j < cells_along_dim; ++j)
      {
        std::vector<unsigned int> cell_nodes;
        unsigned int start_index = i*(celltype + 1) + j*(grid_size)*(celltype + 1);
        for (unsigned int k = 0; k <= celltype +1; ++k)
          for (unsigned int l = 0; l <= celltype + 1; ++l)
          {
            cell_nodes.push_back(grid[start_index + k*grid_size + l]);
          }
        this->cells.push_back(Cell<dim, spacedim, celltype>(this, this->cells.size(), cell_nodes));
      }
    break;
  case 3:
    break;
  }
}

template <unsigned int dim, unsigned int spacedim, unsigned int celltype>
void dromon::Mesh<dim, spacedim, celltype>::insert_cell(
    const std::array<unsigned int,
                     GeometryInfo<dim, spacedim, celltype>::nodes_per_cell>
        &cell)
{
  this->cells.push_back(Cell<dim, spacedim, celltype>(this, this->cells.size(), cell));
}

template<unsigned int dim, unsigned int spacedim, unsigned int celltype>
unsigned int dromon::Mesh<dim, spacedim, celltype>::n_cells() const {
  return this->cells.size();
}

template<unsigned int dim, unsigned int spacedim, unsigned int celltype>
unsigned int dromon::Mesh<dim, spacedim, celltype>::n_faces() const {
  return this->half_faces.size();
}

template<unsigned int dim, unsigned int spacedim, unsigned int celltype>
void dromon::Mesh<dim, spacedim, celltype>::deactivate_cell(const unsigned int cell_index) {\
  assert(cell_index < this->cells.size() && "Cell_index out of range!");
  this->cells[cell_index].deactivate();
}

template<unsigned int dim, unsigned int spacedim, unsigned int celltype>
void dromon::Mesh<dim, spacedim, celltype>::activate_cell(const unsigned int cell_index) {
  assert(cell_index < this->cells.size() && "Cell_index out of range!");
  this->cells[cell_index].activate();
}

template<unsigned int dim, unsigned int spacedim, unsigned int celltype>
dromon::GeomBase<dim, dim, spacedim> *
dromon::Mesh<dim, spacedim, celltype>::get_pointer_to_cell(const unsigned int &index) const {
  //return (index < this->cells.size()) ? const_cast<Cell<dim, spacedim, celltype>*>(&this->cells[index]) : nullptr;
  return const_cast<Cell<dim, spacedim, celltype>*>(&this->cells[index]);
}

template<unsigned int dim, unsigned int spacedim, unsigned int celltype>
dromon::GeomBase<dim-1, dim, spacedim> *
dromon::Mesh<dim, spacedim, celltype>::get_pointer_to_face(const unsigned int &index) const {
  return const_cast<HalfFace<dim, spacedim>*>(&this->half_faces[index]);
}

template<unsigned int dim, unsigned int spacedim, unsigned int celltype>
void dromon::Mesh<dim, spacedim, celltype>::spawn_and_assign_faces() {
  assert(cells.size() != 0 && "The mesh contains no cells!");
  for (auto& cell : cells)
  {
    std::vector<int> temp_faces(2*dim);
    //for the 2*dim faces
    for (unsigned int i = 0; i < 2*dim; ++i)
    {
      //spawn new face
      //...
      std::vector<unsigned int> current_face;
      if (dim != 3)
      {

        int start = i%2 == 0 ? 0 : i == 2*dim -1 ? (celltype+2)*(celltype+2) - (celltype+2) : i*(celltype+1);
        int factor = (i % (2*dim - 1) == 0) ? celltype+1 : (celltype+1)*(celltype + 2);
        for (unsigned int j = 0; j < dim; ++j)
        {
          //add to the current face
          current_face.push_back(cell.get_node_index(start + j*factor));
        }
        HalfFace<dim, spacedim> temp_face(this, current_face);
        this->push_back_face(&temp_face, cell.index);
      }
      else{
        //dim == 3
        for (unsigned int j = 0; j <= dim; ++j)
        {
          assert(false && "Not implemented yet!");
        }
      }

    }
  }
}


//template<unsigned int dim, unsigned int spacedim, unsigned int celltype>
//typename dromon::IteratorRanger<dromon::Mesh<dim, spacedim, celltype>::ActiveCellIterator>
//dromon::Mesh<dim, spacedim, celltype>::active_cell_iterators() const {
//    return IteratorRanger<ActiveCellIterator>(this->begin_active_cell(), this->end_cell());
//}

// -------------------------------------------------------------------
// -- Explicit specializations

//dim-linear elements
template class dromon::Mesh<1,1,0>;
template class dromon::Mesh<1,2,0>;
template class dromon::Mesh<1,3,0>;

template class dromon::Mesh<2,2,0>;
template class dromon::Mesh<2,3,0>;

template class dromon::Mesh<3,3,0>;

//dim-quadratic elements
template class dromon::Mesh<1,1,1>;
template class dromon::Mesh<1,2,1>;
template class dromon::Mesh<1,3,1>;

template class dromon::Mesh<2,2,1>;
template class dromon::Mesh<2,3,1>;

template class dromon::Mesh<3,3,1>;

//dim-cubic elements
template class dromon::Mesh<1,1,2>;
template class dromon::Mesh<1,2,2>;
template class dromon::Mesh<1,3,2>;

template class dromon::Mesh<2,2,2>;
template class dromon::Mesh<2,3,2>;

template class dromon::Mesh<3,3,2>;

