//
// Created by Jake J. Harmon (jake.harmon@ieee.org) on 10/5/21.
//

#ifndef DROMON_DOFBASE_H
#define DROMON_DOFBASE_H

#include "config.h"
#include "MultiIndex.h"
#include "DoFGeomBase.h"
#include <iostream>
#include "Point.h"


DROMON_NAMESPACE_OPEN

//Forward declarations
template<unsigned int dim, class Real>
class FiniteElement;

enum DoFType
{
  CellDoF, FaceDoF
};

inline std::ostream& operator<<(std::ostream& os, const DoFType& dof_type)
{
  if (dof_type == CellDoF)
    return os << "CellDoF";
  else if (dof_type == FaceDoF)
    return os << "FaceDoF";
  else {
    assert(false && "Operator not (currently) valid for choice of DoFType!");
    return os;
  }
}


template <unsigned int dim, class Real = double>
class DoFBase{
public:
  DoFBase() = default;
  DoFBase(DoFGeomBase<dim> *host, int index_on_host, int global_index, DoFDirection direction, DoFType dof_type, CurrentType cur_type,
          MultiIndex<dim> orders);
  //Pointer to DoFGeomBase that the DoFBase belongs to
  DoFGeomBase<dim>* host;
  int index_on_host;

  const DoFDirection direction;
  const DoFType dof_type;
  const CurrentType cur_type;

  const MultiIndex<dim> orders;

  MultiIndex<dim, double> coordinate_orientation;
  bool invert_value;

  // Note: The global index is SIGNED, so that in the event of mismatched directions
  // (i.e., in terms
  // of the unit vectors on two cells), the correct results can be achieved
  int global_index;

  // Note: The sibling global index is SIGNED, so that in the event of mismatched directions
  // (i.e., in terms
  // of the unit vectors on two cells), the correct results can be achieved
  int sibling_global_index; //If the DoF actually exists over two cells, include the index of its sibling

  // Note: The active index is SIGNED, so that in the event of mismatched directions
  // (i.e., in terms
  // of the unit vectors on two cells), the correct results can be achieved
  int active_index;
  bool is_active;
  bool is_new = true;

  template <unsigned int dimY, class RealY>
  friend std::ostream& operator<<(std::ostream& os, const DoFBase<dimY, RealY>& dof);

  const unsigned int get_order(const unsigned int& i) const;
  const double get_coordinate_orientation(const unsigned int& i) const;

  FiniteElement<dim, Real>*  get_active_fe_pointer();
  void set_active_fe_pointer(FiniteElement<dim, Real>* fe_ptr);

  Real evaluate_shape_function(const Point<dim, Real>& uv) const;
  Real evaluate_shape_function_divergence(const Point<dim, Real>& uv) const;

private:
  FiniteElement<dim, Real>* active_fe_pointer = nullptr;
};


template<unsigned int dim, class Real>
DoFBase<dim, Real>::DoFBase(DoFGeomBase<dim> *host, int index_on_host, int global_index, DoFDirection direction, DoFType dof_type, CurrentType cur_type,
                      MultiIndex<dim> orders) : direction(direction), dof_type(dof_type), cur_type(cur_type), orders(orders), invert_value(false), coordinate_orientation(1) {
  this->host = host;
  this->index_on_host = index_on_host;
  this->global_index = global_index;
}

template<unsigned int dimY>
std::ostream &operator<<(std::ostream &os, const DoFBase<dimY> &dof) {
  os << " Global index (active index): " << dof.global_index << " (" << dof.active_index << "), Activity: "
     << (dof.is_active == true ? "true" : "false") << ", DoF Type: " << dof.dof_type << ", DoFDirection: " << dof.direction << ", CurrentType: " << dof.cur_type
     << ", Orders: (";
  for (unsigned int i = 0; i < dimY-1; ++i)
    os << dof.get_order(i) << ", ";
  os << dof.get_order(dimY-1) << ")";

  os << ", Coordinate Orientation: (";
  for (unsigned int i = 0; i < dimY-1; ++i)
    os << dof.get_coordinate_orientation(i) << ", ";
  os << dof.get_coordinate_orientation(dimY-1) << ")";

  return os << ", Invert value: " << dof.invert_value << ", Sibling Global Index: " << dof.sibling_global_index;
}

template <unsigned int dim, class Real>
const unsigned int DoFBase<dim, Real>::get_order(const unsigned int &i) const {
  #ifdef DEBUG
  assert(i < dim && "Index i out of range!");
  #endif
  return orders(i);
}
template <unsigned int dim, class Real>
const double
DoFBase<dim, Real>::get_coordinate_orientation(const unsigned int &i) const {
  #ifdef DEBUG
  assert(i < dim && "Index i out of range!");
  #endif
  return coordinate_orientation(i);
}

template <unsigned int dim, class Real>
FiniteElement<dim, Real> *
DoFBase<dim, Real>::get_active_fe_pointer() {
  return active_fe_pointer;
}
template <unsigned int dim, class Real>
void DoFBase<dim, Real>::set_active_fe_pointer(
    FiniteElement<dim, Real> *fe_ptr) {
  this->active_fe_pointer = fe_ptr;
}
template <unsigned int dim, class Real>
Real DoFBase<dim, Real>::evaluate_shape_function(const Point<dim, Real> &uv) const {
  return (this->active_fe_pointer)->evaluate_shape_function(*this, uv);
}
template <unsigned int dim, class Real>
Real DoFBase<dim, Real>::evaluate_shape_function_divergence(const Point<dim, Real> &uv) const {
  return (this->active_fe_pointer)->evaluate_shape_function_divergence(*this, uv);
}

DROMON_NAMESPACE_CLOSE
#endif //DROMON_DOFBASE_H
