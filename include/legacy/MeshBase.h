//
// Created by Jake J. Harmon (jake.harmon@ieee.org) on 7/20/21.
//

#ifndef DROMON_MESHBASE_H
#define DROMON_MESHBASE_H

#include "config.h"
#include "GeomBase.h"
#include <vector>
#include <memory>
#include <iostream>
#include <ostream>
#include "IteratorRanger.h"
#include "Materials.h"

DROMON_NAMESPACE_OPEN

namespace IteratorState
{
    enum IteratorStates
    {
        valid,
        past_the_end,
        invalid
    };
}
/**
 *
 * @tparam spacedim
 * The node is the simplest geometric unit, consisting of $spacedim$ scalars
 */
    template <unsigned int dim, unsigned int spacedim, class Real = double>
    class Node : public GeomBase<0, dim, spacedim>
    {
    public:
        Node(MeshBase<dim, spacedim>* mesh, Point<spacedim, Real>& point);
        Point<spacedim, Real>& get_point();
        Real operator()(unsigned int pos) const;
        Real& operator[](const unsigned int& pos);
        template <unsigned int dimS, unsigned int spacedimS, class RealS>
        friend std::ostream& operator<<(std::ostream& os, const Node<dimS, spacedimS, RealS>& nt);
        Node operator*(Real scalar) const;
        void operator*=(Real scalar);
        void operator+=(const Node<dim, spacedim, Real>& n2);

    private:
        Point<spacedim, Real> coord;
    };

    template<unsigned int dim, unsigned int spacedim, class Real>
    Node<dim, spacedim, Real>::Node(MeshBase<dim, spacedim> *mesh, Point<spacedim, Real> &point) : GeomBase<0, dim, spacedim>(mesh, 0) {
        coord = point;
    }

    template<unsigned int dim, unsigned int spacedim, class Real>
    Point<spacedim, Real> &Node<dim, spacedim, Real>::get_point() {
        return coord;
    }



    template<unsigned int dim, unsigned int spacedim, class Real>
    Real Node<dim, spacedim, Real>::operator()(unsigned int pos) const {
        return coord(pos);
    }

    template<unsigned int dim, unsigned int spacedim, class Real>
    Real &Node<dim, spacedim, Real>::operator[](const unsigned int & pos) {
        return coord[pos];
    }


    template<unsigned int dimS, unsigned int spacedimS, class RealS>
    std::ostream &operator<<(std::ostream &os, const Node<dimS, spacedimS, RealS> &nt) {

        for (unsigned int i = 0; i < spacedimS; ++i) {
            os << nt(i) << " ";
        }
        return os;
    }

    template<unsigned int dim, unsigned int spacedim, class Real>
    void Node<dim, spacedim, Real>::operator*=(Real scalar) {
        this->coord*=scalar;
    }

    template<unsigned int dim, unsigned int spacedim, class Real>
    Node<dim, spacedim, Real> Node<dim, spacedim, Real>::operator*(Real scalar) const {
            Node<dim, spacedim, Real> out = *this;
            out*=scalar;
            return out;
    }

    template<unsigned int dim, unsigned int spacedim, class Real>
    void Node<dim, spacedim, Real>::operator+=(const Node<dim, spacedim, Real>& n2) {
        this->coord += n2.coord;
    }







//    template <unsigned int dim, unsigned int spacedim, class Real>
//    std::ostream &operator<<(std::ostream &os, const Node<dim, spacedim, Real> &nt) {
//        for (unsigned int i = 0; i < spacedim; ++i)
//            os << nt.get_point()[i] << " ";
//        return os;
//    }

    template <unsigned int dim, unsigned int spacedim>
class MeshBase
{
public:
    virtual Node<dim, spacedim, double>* get_node(const unsigned int& i);
    virtual unsigned int n_nodes() const;
    virtual void push_back_node(Node<dim, spacedim, double>* node_pb);
    std::vector<std::unique_ptr<Node<dim, spacedim, double>>> nodes;
    unsigned int material_domain_id = 0;


    virtual unsigned int n_cells() const;
    virtual unsigned int n_faces() const;

    virtual unsigned int cell_type() const;

    virtual GeomBase<dim, dim, spacedim>* get_pointer_to_cell(const unsigned int& index) const;

    virtual GeomBase<dim-1, dim, spacedim>* get_pointer_to_face(const unsigned int& index) const;
    /*
    * Iterators
    * */
    struct ActiveCellIterator
    {
        using iterator_category = std::forward_iterator_tag;
        using difference_type = std::ptrdiff_t;
        using value_type = GeomBase<dim, dim, spacedim>;
        using pointer = GeomBase<dim, dim, spacedim>*;
        using reference = GeomBase<dim, dim, spacedim>&;
        ActiveCellIterator(MeshBase<dim, spacedim>* mesh, pointer ptr, const unsigned int& present_index = 0) : mesh(mesh), m_ptr(ptr), present_index(present_index) {
            iterator_state = IteratorState::valid;
        }
        ActiveCellIterator()
        {
            iterator_state = IteratorState::invalid;
        }

        const reference operator*() const {return *m_ptr;}
        const pointer operator->() const {return m_ptr;}

        reference operator*() {return *m_ptr;}
        pointer operator->() {return m_ptr;}
        //Prefix increment
        ActiveCellIterator& operator++()
        {
            do
            {
                ++present_index;
               // ++m_ptr;
                m_ptr = mesh->get_pointer_to_cell(present_index);
            }while(present_index < mesh->n_cells() && !m_ptr->is_active());
            return *this;
        }
        //Postfix increment
        ActiveCellIterator operator++(int){
            auto tmp = *this;
            ++(*this);
            return tmp;
        }

        friend bool operator==(const ActiveCellIterator& a, const ActiveCellIterator& b) {return a.m_ptr == b.m_ptr;}
        friend bool operator!=(const ActiveCellIterator& a, const ActiveCellIterator& b){ return a.m_ptr != b.m_ptr;}


    private:
        //The present index allows (since this iterator is exclusively a forward iterator), keeping track of when
        //we have exceeded the range of the iterator...
        unsigned int present_index;
        pointer m_ptr;
        const MeshBase<dim, spacedim>* mesh;
        IteratorState::IteratorStates iterator_state;


    };

    //virtual ActiveCellIterator end_cell() const {return ActiveCellIterator(const_cast<MeshBase<dim, spacedim>*>(this), const_cast<GeomBase<dim, dim, spacedim>*>(&cells[cells.size()]), n_cells());}
    virtual ActiveCellIterator end_cell() const {return ActiveCellIterator();}
    virtual ActiveCellIterator begin_active_cell() const {return ActiveCellIterator();}
//    virtual ActiveCellIterator begin_active_cell() const {
//        for (unsigned int i = 0; i < this->cells.size(); ++i)
//            if (cells[i].is_active())
//                return ActiveCellIterator(const_cast<Mesh<dim, spacedim, celltype>*>(this), const_cast<Cell<dim, spacedim, celltype>*>(&cells[i]), i);
//        return end_cell();
//    }


    virtual IteratorRanger<ActiveCellIterator> active_cell_iterators() const
    {
        return IteratorRanger<ActiveCellIterator>(begin_active_cell(), end_cell());
    }


private:
    bool active = true;
};

    template<unsigned int dim, unsigned int spacedim>
    Node<dim, spacedim, double> *MeshBase<dim, spacedim>::get_node(const unsigned int& i) {
#ifdef DEBUG
        assert(i < this->nodes.size() && "Index i out of bounds for number of nodes!");
#endif
        return &*this->nodes[i];
    }

    template<unsigned int dim, unsigned int spacedim>
    unsigned int MeshBase<dim, spacedim>::n_nodes() const {
        return this->nodes.size();
    }

    template<unsigned int dim, unsigned int spacedim>
    void MeshBase<dim, spacedim>::push_back_node(Node<dim, spacedim, double> *node_pb) {
        this->nodes.push_back(std::make_unique<Node<dim, spacedim, double>>(*node_pb));
    }

    template<unsigned int dim, unsigned int spacedim>
    unsigned int MeshBase<dim, spacedim>::n_cells() const {
        return 0;
    }

    template<unsigned int dim, unsigned int spacedim>
    unsigned int MeshBase<dim, spacedim>::n_faces() const {
        return 0;
    }

    template<unsigned int dim, unsigned int spacedim>
    unsigned int MeshBase<dim, spacedim>::cell_type() const {
        return 0;
    }

    template<unsigned int dim, unsigned int spacedim>
    GeomBase<dim, dim, spacedim> *MeshBase<dim, spacedim>::get_pointer_to_cell(const unsigned int& index) const {
        return nullptr;
    }

    template<unsigned int dim, unsigned int spacedim>
    GeomBase<dim - 1, dim, spacedim> *MeshBase<dim, spacedim>::get_pointer_to_face(const unsigned int &index) const {
        return nullptr;
    }


DROMON_NAMESPACE_CLOSE

#endif //DROMON_MESHBASE_H
