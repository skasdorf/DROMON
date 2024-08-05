//
// Created by Jake J. Harmon (jake.harmon@ieee.org) on 10/5/21.
//

#include "DoFBase.h"
#include "FEBase.h"


//template<unsigned int dim, class Real>
//dromon::DoFBase<dim, Real>::DoFBase(DoFGeomBase<dim> *host, int index_on_host, int global_index, DoFDirection direction, DoFType dof_type, CurrentType cur_type,
//                            MultiIndex<dim> orders) : direction(direction), dof_type(dof_type), cur_type(cur_type), orders(orders), invert_value(false), coordinate_orientation(1) {
//  this->host = host;
//  this->index_on_host = index_on_host;
//  this->global_index = global_index;
//}


//template <unsigned int dim, class Real>
//const unsigned int dromon::DoFBase<dim, Real>::get_order(const unsigned int &i) const {
//  assert(i < dim && "Index i out of range!");
//  return orders(i);
//}
//template <unsigned int dim, class Real>
//const double
//dromon::DoFBase<dim, Real>::get_coordinate_orientation(const unsigned int &i) const {
//  assert(i < dim && "Index i out of range!");
//  return coordinate_orientation(i);
//}
//template <unsigned int dim, class Real>
//dromon::FiniteElement<dim, Real> *
//dromon::DoFBase<dim, Real>::get_active_fe_pointer() {
//  return active_fe_pointer;
//}
//template <unsigned int dim, class Real>
//void dromon::DoFBase<dim, Real>::set_active_fe_pointer(
//    dromon::FiniteElement<dim, Real> *fe_ptr) {
//  this->active_fe_pointer = fe_ptr;
//}
//template <unsigned int dimY, class RealY>
//std::ostream &operator<<(std::ostream &os,
//                                   const dromon::DoFBase<dimY, RealY> &dof) {
//  os << " Global index (active index): " << dof.global_index << " (" << dof.active_index << "), Activity: "
//     << (dof.is_active == true ? "true" : "false") << ", DoF Type: " << dof.dof_type << ", DoFDirection: " << dof.direction << ", CurrentType: " << dof.cur_type
//     << ", Orders: (";
//  for (unsigned int i = 0; i < dimY-1; ++i)
//    os << dof.get_order(i) << ", ";
//  os << dof.get_order(dimY-1) << ")";
//
//  os << ", Coordinate Orientation: (";
//  for (unsigned int i = 0; i < dimY-1; ++i)
//    os << dof.get_coordinate_orientation(i) << ", ";
//  os << dof.get_coordinate_orientation(dimY-1) << ")";
//
//  return os << ", Invert value: " << dof.invert_value << ", Sibling Global Index: " << dof.sibling_global_index;
//}

//template<unsigned int dimY, class RealY>
//inline std::ostream &operator<<(std::ostream &os, const dromon::DoFBase<dimY, RealY> &dof) {
//  os << " Global index (active index): " << dof.global_index << " (" << dof.active_index << "), Activity: "
//     << (dof.is_active == true ? "true" : "false") << ", DoF Type: " << dof.dof_type << ", DoFDirection: " << dof.direction << ", CurrentType: " << dof.cur_type
//     << ", Orders: (";
//  for (unsigned int i = 0; i < dimY-1; ++i)
//    os << dof.get_order(i) << ", ";
//  os << dof.get_order(dimY-1) << ")";
//
//  os << ", Coordinate Orientation: (";
//  for (unsigned int i = 0; i < dimY-1; ++i)
//    os << dof.get_coordinate_orientation(i) << ", ";
//  os << dof.get_coordinate_orientation(dimY-1) << ")";
//
//  return os << ", Invert value: " << dof.invert_value << ", Sibling Global Index: " << dof.sibling_global_index;
//}
// -------------------------------------------------------------------
// -- Explicit specializations

//In the event of needing to use other floating-point representations, e.g., float or quad,
//add the specializations as below...
template class dromon::FiniteElement<1, double>;
template class dromon::FiniteElement<1, float>;

template class dromon::FiniteElement<2, double>;
template class dromon::FiniteElement<2, float>;

template class dromon::FiniteElement<3, double>;
template class dromon::FiniteElement<3, float>;


