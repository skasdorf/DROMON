//
// Created by Jake J. Harmon (jake.harmon@ieee.org) on 8/7/21.
//

#ifndef DROMON_FECOLLECTION_H
#define DROMON_FECOLLECTION_H

#include "FEBase.h"
#include "config.h"
#include <map>
#include <memory>
#include <vector>

DROMON_NAMESPACE_OPEN

// TODO: Support anisotropy in $p$
template <unsigned int dim, unsigned int spacedim, class Real = double>
class FECollection {
public:
  template <class FE> void push_back(FE *entry);
  // std::map<CurrentType, bool> types;
  CurrentType type = Nothing;
  unsigned int max(unsigned int dim_index) const;
  FiniteElement<dim, Real> *get_fe(unsigned int index);
  const unsigned int get_fe_degree(unsigned int index);
  unsigned int next_fe(const unsigned int& current_index) const;
  unsigned int next_fe(const unsigned int& current_index, const unsigned int& num_p_refinements) const;
  unsigned int prev_fe(const unsigned int& current_index) const;
  unsigned int prev_fe(const unsigned int& current_index, const unsigned int& num_p_coarsenings) const;
  friend bool operator==(const FECollection<dim, spacedim, Real> &m1,
                         const FECollection<dim, spacedim, Real> &m2) {
    if (m1.fe_collection.size() != m2.fe_collection.size())
      return false;
    for (unsigned int i = 0; i < m1.fe_collection.size(); ++i)
      if (m1.fe_collection[i] != m2.fe_collection[i])
        return false;
    return true;
  }
  friend bool operator!=(const FECollection<dim, spacedim, Real> &m1,
                         const FECollection<dim, spacedim, Real> &m2) {
    return !(m1 == m2);
  }

private:
  unsigned int max_degree = 0;
  std::vector<std::shared_ptr<FiniteElement<dim,  Real>>>
      fe_collection;
};
// TODO: Support anisotropy in $p$
template <unsigned int dim, unsigned int spacedim, class Real>
template <class FE>
void FECollection<dim, spacedim, Real>::push_back(FE *entry) {
  static_assert(std::is_base_of<FiniteElement<dim, Real>, FE>::value,
                "T1 must derive from Base");
  fe_collection.push_back(std::make_shared<FE>(*entry));
  // types[entry.get_type()] = true;
  if (type == Nothing)
    type = entry->get_type();
  else
    assert(type == entry->get_type() &&
           "Illegal assigment of a Finite Element to the collection! All "
           "entries must be of the same type!");

  assert(entry->get_degree() > this->max_degree && "Finite Elements MUST be inserted in ascending order!");
  if (entry->get_degree() > this->max_degree)
    this->max_degree = entry->get_degree();
}
// TODO: Support anisotropy in $p$
template <unsigned int dim, unsigned int spacedim, class Real>
unsigned int
FECollection<dim, spacedim, Real>::max(unsigned int dim_index) const {
  return dim_index < dim ? max_degree : 0;
}

template <unsigned int dim, unsigned int spacedim, class Real>
FiniteElement<dim, Real> *
FECollection<dim, spacedim, Real>::get_fe(unsigned int index) {
  #ifdef DEBUG
  assert(index < fe_collection.size() && "Index out of range!");
  #endif
  return fe_collection[index].get();
}

template <unsigned int dim, unsigned int spacedim, class Real>
const unsigned int
FECollection<dim, spacedim, Real>::get_fe_degree(unsigned int index) {
#ifdef DEBUG
  assert(index < fe_collection.size() && "Index out of range!");
#endif
  return fe_collection[index]->get_degree();
}
template <unsigned int dim, unsigned int spacedim, class Real>
unsigned int FECollection<dim, spacedim, Real>::next_fe(const unsigned int& current_index) const {
  if (current_index + 1 < fe_collection.size())
    return current_index+1;
  else
    return current_index;
}
template <unsigned int dim, unsigned int spacedim, class Real>
unsigned int FECollection<dim, spacedim, Real>::next_fe(
    const unsigned int &current_index,
    const unsigned int& num_p_refinements) const {
  unsigned int out = current_index;
  for (unsigned int i =0; i < num_p_refinements; ++i)
    out = next_fe(out);
  return out;
}
template <unsigned int dim, unsigned int spacedim, class Real>
unsigned int FECollection<dim, spacedim, Real>::prev_fe(
    const unsigned int &current_index) const {
  if (current_index - 1 >= 0)
    return current_index-1;
  else
    return current_index;
}
template <unsigned int dim, unsigned int spacedim, class Real>
unsigned int FECollection<dim, spacedim, Real>::prev_fe(
    const unsigned int &current_index,
    const unsigned int &num_p_coarsenings) const {
  unsigned int out = current_index;
  for (unsigned int i =0; i < num_p_coarsenings; ++i)
    out = prev_fe(out);
  return out;
}

// The FECollectionCollector is for storing FECollections for the different
// types of current expansions (e.g., electric or magnetic)
// TODO: Support anisotropy in $p$
template <unsigned int dim, unsigned int spacedim, class Real = double>
class FECollectionCollector {
public:
  void push_back(const FECollection<dim, spacedim, Real> &entry);
  unsigned int get_number_of_types();
  CurrentType get_type(unsigned int &index);
  std::vector<CurrentType> get_current_types();
  FECollection<dim, spacedim, Real> *
  get_fe_collection(const unsigned int &index);
  FECollection<dim, spacedim, Real> *
  get_fe_collection(const CurrentType &cur_type);

  friend bool operator==(const FECollectionCollector<dim, spacedim, Real> &m1,
                         const FECollectionCollector<dim, spacedim, Real> &m2) {
    if (m1.fe_collection_collector.size() != m2.fe_collection_collector.size())
      return false;
    for (unsigned int i = 0; i < m1.fe_collection_collector.size(); ++i)
      if (m1.fe_collection_collector[i] != m2.fe_collection_collector[i])
        return false;
    return true;
  }
  friend bool operator!=(const FECollectionCollector<dim, spacedim, Real> &m1,
                         const FECollectionCollector<dim, spacedim, Real> &m2) {
    return !(m1 == m2);
  }

private:
  std::vector<FECollection<dim, spacedim, Real>> fe_collection_collector;
  // std::vector<CurrentType> types;
  std::map<CurrentType, unsigned int> type_to_index;
};

template <unsigned int dim, unsigned int spacedim, class Real>
void FECollectionCollector<dim, spacedim, Real>::push_back(
    const FECollection<dim, spacedim, Real> &entry) {
  fe_collection_collector.push_back(entry);
  //types.push_back(entry.type);
  assert(type_to_index.find(entry.type) == type_to_index.end() && "Adding an FECollection of a current type that already belongs to the FECollectionCollector is illegal!");
  type_to_index[entry.type] = fe_collection_collector.size()-1;
}

template <unsigned int dim, unsigned int spacedim, class Real>
unsigned int FECollectionCollector<dim, spacedim, Real>::get_number_of_types() {
  return type_to_index.size();
}

template <unsigned int dim, unsigned int spacedim, class Real>
CurrentType
FECollectionCollector<dim, spacedim, Real>::get_type(unsigned int &index) {
  assert(index < type_to_index.size() && "Index out of range!");
  return fe_collection_collector[index].type;
}

template <unsigned int dim, unsigned int spacedim, class Real>
FECollection<dim, spacedim, Real> *
FECollectionCollector<dim, spacedim, Real>::get_fe_collection(
    const unsigned int &index) {
#ifdef DEBUG
  assert(index < fe_collection_collector.size() && "Index out of range!");
#endif
  return &(fe_collection_collector.at(index));
}

template <unsigned int dim, unsigned int spacedim, class Real>
FECollection<dim, spacedim, Real> *
FECollectionCollector<dim, spacedim, Real>::get_fe_collection(
    const CurrentType &cur_type)
{
  return (get_fe_collection(type_to_index.at(cur_type)));
}
template <unsigned int dim, unsigned int spacedim, class Real>
std::vector<CurrentType>
FECollectionCollector<dim, spacedim, Real>::get_current_types() {
  std::vector<CurrentType> out;
  for (auto it = type_to_index.begin(); it != type_to_index.end(); ++it)
  {
    out.push_back(it->first);
  }
  return out;
}

DROMON_NAMESPACE_CLOSE
#endif // DROMON_FECOLLECTION_H
