//
// Created by Jake J. Harmon (jake.harmon@ieee.org) on 7/20/21.
//

#include "MeshBase.h"

//template<unsigned int dim, unsigned int spacedim>
//dromon::Node<dim, spacedim, double> *dromon::MeshBase<dim, spacedim>::get_node(unsigned int i) {
//    assert(i < this->nodes.size() && "Index i out of bounds for number of nodes!");
//    return &*this->nodes[i];
//}

//template<unsigned int dim, unsigned int spacedim>
//unsigned int dromon::MeshBase<dim, spacedim>::n_nodes() const {
//    return this->nodes.size();
//}

//template<unsigned int dim, unsigned int spacedim>
//void dromon::MeshBase<dim, spacedim>::push_back_node(dromon::Node<dim, spacedim, double> *node_pb) {
//    this->nodes.push_back(std::make_unique<Node<dim, spacedim, double>>(*node_pb));
//}



// -------------------------------------------------------------------
// -- Explicit specializations

template class dromon::MeshBase<1,1>;
template class dromon::MeshBase<1,2>;
template class dromon::MeshBase<1,3>;

template class dromon::MeshBase<2,2>;
template class dromon::MeshBase<2,3>;

template class dromon::MeshBase<3,3>;