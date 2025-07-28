//
// Created by Jake J. Harmon (jake.harmon@ieee.org) on 8/6/21.
//

#ifndef DROMON_MULTIINDEX_H
#define DROMON_MULTIINDEX_H

#include "config.h"
#include <vector>

#include <cstdarg>

DROMON_NAMESPACE_OPEN

template <unsigned int dim, class T = unsigned int>
class MultiIndex
{
public:
    //explicit MultiIndex(const std::vector<unsigned int>& indices);
    MultiIndex() = default;


    MultiIndex(const T& i) /*: indices(dim)*/ {
        for (auto& index : indices)
          index = i;
    }
    MultiIndex(const T& i, const T& j) /*: indices(dim)*/ {
        #ifdef DEBUG
        assert(dim == 2 && "Incorrect constructor for the dim of the MultiIndex!");
        #endif
        indices[0] = i;
        indices[1] = j;
    }
    MultiIndex(const T& i, const T& j, const T& k) /*: indices(dim)*/ {
        #ifdef DEBUG
        assert(dim == 3 && "Incorrect constructor for the dim of the MultiIndex!");
        #endif
        indices[0] = i;
        indices[1] = j;
        indices[2]= k;
    }


    const T& operator()(const unsigned int& component) const
    {
#ifdef DEBUG
        assert(component < dim && "Component out of range!");
#endif
        return indices[component];
    }

    T& operator[](const unsigned int& component)
    {
#ifdef DEBUG
        assert(component < dim && "Component out of range!");
#endif
        return indices[component];
    }

    const T& at(const unsigned int& component) const
    {
#ifdef DEBUG
      assert(component < dim && "Component out of range!");
#endif
      return indices[component];
    }

    friend bool operator==(const MultiIndex<dim, T>& m1, const MultiIndex<dim, T>& m2)
    {
        for (unsigned int i = 0; i < dim; ++i)
            if (m1.indices[i] != m2.indices[i])
                return false;
        return true;
    }
    friend bool operator!=(const MultiIndex<dim, T>& m1, const MultiIndex<dim, T>& m2)
    {
        return !(m1 == m2);
    }

private:
 // std::vector<T> indices;
  T indices[dim];
};





DROMON_NAMESPACE_CLOSE
#endif //DROMON_MULTIINDEX_H
