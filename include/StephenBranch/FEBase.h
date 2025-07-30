//
// Created by Jake J. Harmon (jake.harmon@ieee.org) on 8/6/21.
//

#ifndef DROMON_FEBASE_H
#define DROMON_FEBASE_H

#include "config.h"
#include "Point.h"
#include "MultiIndex.h"
#include "DoFBase.h"

DROMON_NAMESPACE_OPEN


    enum Conformity {
        Hdiv = 0,
        Hcurl = 1
    };

    template<unsigned int dim, class Real = double>
    class FiniteElement {
    public:


        static const unsigned int dimension = dim;
        // static const unsigned int space_dimension = spacedim;

        FiniteElement() = default;

        FiniteElement(const unsigned int deg, const unsigned int components, const Conformity conforming_type, const CurrentType current_type);

        virtual Real evaluate_shape_function(const DoFDirection& direction, const MultiIndex<dim> &indices,
                                             const Point<dim, Real> &uv){
            return Real(0);
        }
        virtual Real evaluate_shape_function_divergence(const DoFDirection& direction, const MultiIndex<dim> &indices,
                                                        const Point<dim, Real> &uv){
          return Real(0);
        }

        virtual Real evaluate_shape_function(const DoFBase<dim>& dof, const Point<dim, Real>& uv){
          return Real(0);
        }

        virtual Real evaluate_shape_function_divergence(const DoFBase<dim>& dof,
                                                        const Point<dim, Real> &uv){
          return Real(0);
        }

        virtual const CurrentType &get_type() const;

        virtual const unsigned int &get_degree() const;

        virtual const unsigned int &get_n_components() const;

        virtual const Conformity &get_conformity() const;

        virtual std::string get_FE_name() const;

        //In the event of anisotropic orders, max_expansion_order should be a MultiIndex as well...
        virtual bool is_subset_of_max(const MultiIndex<dim> &indices, DoFDirection dir, const unsigned int max_expasion_order) const;

        friend bool operator==(const FiniteElement<dim, Real>& m1, const FiniteElement<dim, Real>& m2)
        {
            if (m1.get_FE_name() == m2.get_FE_name()
            && m1.get_type() == m2.get_type()
            && m1.get_degree() == m2.get_degree()
            && m1.get_n_components() == m2.get_n_components()
            && m1.get_conformity() == m2.get_conformity())
                return true;
            else return false;
        }
        friend bool operator!=(const FiniteElement<dim, Real>& m1, const FiniteElement<dim, Real>& m2)
        {
            return !(m1 == m2);
        }

        virtual bool contains_this_dof(const DoFBase<dim>& dof) const;

        //virtual std::pair<MultiIndex<dim, double>, bool> determine_pairing(const unsigned int& local_face_index, const unsigned int & neighbor_face_index, const bool& matched_face_orientation) const;
        virtual void determine_pairing(const unsigned int& local_face_index, const unsigned int & neighbor_face_index, const bool& matched_face_orientation, MultiIndex<dim, double>* coordinate_orientation, bool* invert_value) const;

    private:
        const CurrentType type;
        const Conformity conforming_space;
        const unsigned n_components;
        const unsigned int degree;
    };

    template<unsigned int dim,class Real>
    FiniteElement<dim, Real>::FiniteElement(const unsigned int deg, const unsigned int components,
                                                      const Conformity conforming_type, const CurrentType current_type) : degree(deg),
                                                                                          n_components(components),
                                                                                          conforming_space(
                                                                                                  conforming_type),
                                                                                                  type(current_type){}

    template<unsigned int dim, class Real>
    const CurrentType &FiniteElement<dim, Real>::get_type() const {
        return type;
    }

    template<unsigned int dim, class Real>
    const unsigned int &FiniteElement<dim, Real>::get_degree() const {
        return degree;
    }

    template<unsigned int dim, class Real>
    const unsigned int &FiniteElement<dim, Real>::get_n_components() const {
        return n_components;
    }

    template<unsigned int dim, class Real>
    const Conformity &FiniteElement<dim, Real>::get_conformity() const {
        return conforming_space;
    }

    template<unsigned int dim,  class Real>
    std::string FiniteElement<dim, Real>::get_FE_name() const {
        return std::string();
    }

    template<unsigned int dim,  class Real>
    bool FiniteElement<dim, Real>::is_subset_of_max(const MultiIndex<dim> &indices, DoFDirection dir,
                                                              const unsigned int max_expasion_order) const {
        return true;
    }
    template <unsigned int dim,  class Real>
    bool FiniteElement<dim, Real>::contains_this_dof(
        const DoFBase<dim> &dof) const {
        std::cout << "Should not be getting here!" << std::endl;
        return false;
    }
//    template <unsigned int dim,  class Real>
//    std::pair<MultiIndex<dim, double>, bool>
//    FiniteElement<dim, spacedim, Real>::determine_pairing(
//        const unsigned int &local_face_index,
//        const unsigned int &neighbor_face_index,
//        const bool& matched_face_orientation) const {
//      return std::make_pair<MultiIndex<dim, double>, bool>(MultiIndex<dim, double>(0), false);
//    }

    template <unsigned int dim,  class Real>
    void
    FiniteElement<dim, Real>::determine_pairing(
        const unsigned int &local_face_index,
        const unsigned int &neighbor_face_index,
        const bool& matched_face_orientation, MultiIndex<dim, double>* coordinate_orientation, bool* invert_value) const {
      return;
    }

    DROMON_NAMESPACE_CLOSE
#endif //DROMON_FEBASE_H
