//
// Created by Jake J. Harmon (jake.harmon@ieee.org) on 8/3/21.
//

/*
 * Based on the implementation in Deal.II of the IteratorRange and IteratorOverIterators.
 * Provides class templates for taking standard iterators (i.e., the active cell iterator)
 * and permitting range-based for loops -- for (auto& x : items) ...
 */
#ifndef DROMON_ITERATORRANGER_H
#define DROMON_ITERATORRANGER_H

#include "config.h"

DROMON_NAMESPACE_OPEN
//Forward declaration
template <class Iterator>
        class IteratorOverIterators;

template <class Iterator>
class IteratorRanger{
public:
    using IteratorOverIteratorsInteral = IteratorOverIterators<Iterator>;
    using iterator = Iterator;
    IteratorRanger();
    IteratorRanger(const iterator begin, const iterator end);

    IteratorOverIteratorsInteral begin();
    IteratorOverIteratorsInteral begin() const;

    IteratorOverIteratorsInteral end();
    IteratorOverIteratorsInteral end() const;

private:
    const IteratorOverIteratorsInteral it_begin;
    const IteratorOverIteratorsInteral it_end;
};

template <class Iterator>
class IteratorOverIterators
{
public:
    using BaseIterator = Iterator;

    explicit IteratorOverIterators(const BaseIterator& iterator);
    IteratorOverIterators()= default;
    const BaseIterator  &operator*() const;
    const BaseIterator *operator->() const;

    IteratorOverIterators& operator++();
    IteratorOverIterators operator++(int);

    bool operator!=(const IteratorOverIterators& i_o_i) const;

    operator const BaseIterator &() const;

    using iterator_category = std::forward_iterator_tag;
    using value_type = Iterator;
    using difference_type = typename Iterator::difference_type;
    using point = Iterator*;
    using reference = Iterator&;

private:
    BaseIterator element_of_iterator_collection;

};


    template <typename Iterator>
    inline IteratorOverIterators<Iterator>::IteratorOverIterators(
            const BaseIterator &iterator)
            : element_of_iterator_collection(iterator)
    {}



    template <typename Iterator>
    inline const typename IteratorOverIterators<Iterator>::BaseIterator &
    IteratorOverIterators<Iterator>::operator*() const
    {
        return element_of_iterator_collection;
    }



    template <typename Iterator>
    inline const typename IteratorOverIterators<Iterator>::BaseIterator *
    IteratorOverIterators<Iterator>::operator->() const
    {
        return &element_of_iterator_collection;
    }



    template <typename Iterator>
    inline IteratorOverIterators<Iterator> &
    IteratorOverIterators<Iterator>::operator++()
    {
        ++element_of_iterator_collection;
        return *this;
    }



    template <typename Iterator>
    inline IteratorOverIterators<Iterator>
    IteratorOverIterators<Iterator>::operator++(int)
    {
        const IteratorOverIterators old_value = *this;
        ++element_of_iterator_collection;
        return *old_value;
    }



    template <typename Iterator>
    inline bool
    IteratorOverIterators<Iterator>::
    operator!=(const IteratorOverIterators &i_o_i) const
    {
        return element_of_iterator_collection != i_o_i.element_of_iterator_collection;
    }



    template <typename Iterator>
    inline IteratorOverIterators<Iterator>::operator const BaseIterator &() const
    {
        return element_of_iterator_collection;
    }

    template <typename Iterator>
    inline IteratorRanger<Iterator>::IteratorRanger()
    {}



    template <typename Iterator>
    inline IteratorRanger<Iterator>::IteratorRanger(const iterator b,
                                                  const iterator e)
            : it_begin(b)
            , it_end(e)
    {}


    template <typename Iterator>
    inline typename IteratorRanger<Iterator>::IteratorOverIteratorsInteral
    IteratorRanger<Iterator>::begin()
    {
        return it_begin;
    }


    template <typename Iterator>
    inline typename IteratorRanger<Iterator>::IteratorOverIteratorsInteral
    IteratorRanger<Iterator>::begin() const
    {
        return it_begin;
    }


    template <typename Iterator>
    inline typename IteratorRanger<Iterator>::IteratorOverIteratorsInteral
    IteratorRanger<Iterator>::end()
    {
        return it_end;
    }


    template <typename Iterator>
    inline typename IteratorRanger<Iterator>::IteratorOverIteratorsInteral
    IteratorRanger<Iterator>::end() const
    {
        return it_end;
    }
DROMON_NAMESPACE_CLOSE
#endif //DROMON_ITERATORRANGER_H
