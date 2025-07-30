//
// Created by Jake J. Harmon (jake.harmon@ieee.org) on 12/13/21.
//

#ifndef DROMON_DOFGEOMBASE_H
#define DROMON_DOFGEOMBASE_H
#include "config.h"
#include <map>
#include <unordered_map>

DROMON_NAMESPACE_OPEN

namespace DoFGeom {
    enum DoFGeomType
    {
        Cell, Face
    };
}

template <unsigned int dim>
class DoFGeomBase {
public:
    DoFGeomBase(unsigned int index, bool active, unsigned int geom_level, DoFGeom::DoFGeomType type)
    {
        this->index = index;
        this->active = active;
        this->geom_level = geom_level;
        this->type = type;

    }
    const unsigned int& get_index() const {
        return index;
    }
    unsigned int& get_fe_index(unsigned int index){
        return fe_indices.at(index);
    }
    unsigned int& get_fe_index(const CurrentType& cur_type)
    {
      return fe_indices.at(get_relevant_fe_index(cur_type));
    }

//    void set_fe_index(unsigned int fe_index, const CurrentType& cur_type){
//        // assert(type < fe_indices.size() && "Type is out of range!");
//        auto index = get_relevant_fe_index(cur_type);
//        this->fe_indices[index] = fe_index;
//      //  this->active_fe_degrees[index] = active_degree;
//    }
    const  DoFGeom::DoFGeomType& get_type() const{
        return type;
    }
    virtual void initialize_fe_indices(unsigned int starting_index, const std::vector<CurrentType>& cur_types);
    virtual void set_fe_index(const unsigned int& new_index);
    virtual void set_fe_index(const unsigned int& new_index, const CurrentType& cur_type);
    virtual const unsigned int get_neighbor_index() const
    {
      return std::numeric_limits<unsigned int>::max();
    }

    virtual const unsigned int parent_cell_index() const{
      return std::numeric_limits<unsigned int>::max();
    }


    virtual void clear_dofs() = 0;
    virtual void deactivate_dofs() = 0;
    virtual void activate_dofs() = 0;

    //TODO: Generalize for arbitrary dim, i.e. not just for i, j
    virtual const bool is_compatible(const DoFDirection& dir, const unsigned int& i, const unsigned int& j) const = 0;
    virtual unsigned int n_dofs() const = 0;
    virtual unsigned int n_active_dofs() const = 0;


    virtual const unsigned int get_relevant_fe_index(const CurrentType& cur_type) const;
    virtual void mark_as_populated(unsigned int populated_degree);
    virtual const unsigned int active_degree(unsigned int type) const;

    //The list of the fe_indices associated with this DoFGeom
    //This is a vector to support multiple types of expansions on the same object
    //(e.g., magnetic and electric current density vectors)

    std::vector<unsigned int> fe_indices;
    std::vector<unsigned int> active_fe_degrees;



    unsigned int index;
    bool active;
    unsigned int geom_level;
    DoFGeom::DoFGeomType type;
    unsigned n_supported_expansions;
    unsigned int populated_expansion_order = 0;

    // Anytime the expansion order is changed, we must check
    // whether new DoFs actually have to be added (as opposed to just activated/deactivated)
    bool fully_populated = true;

    private:
    std::unordered_map<CurrentType, unsigned int> cur_type_to_fe_index;



};

    template<unsigned int dim>
    void DoFGeomBase<dim>::initialize_fe_indices(unsigned int starting_index, const std::vector<CurrentType>& cur_types) {
           // assert(n_types == cur_types.size() && "N_types and the size of cur_types do not match!");
            const unsigned int n_types = cur_types.size();
            fe_indices = std::vector<unsigned int>(n_types, starting_index);
            //active_fe_degrees =  std::vector<unsigned int>(n_types, active_degree);
            fully_populated = false;
            for (unsigned int i = 0; i < cur_types.size(); ++i)
              cur_type_to_fe_index[cur_types[i]] = i;
    }

    template<unsigned int dim>
    void DoFGeomBase<dim>::mark_as_populated(unsigned int populated_expansion_order) {
        this->fully_populated = true;
        this->populated_expansion_order = std::max(this->populated_expansion_order, populated_expansion_order);
    }
    template <unsigned int dim>
    const unsigned int
    DoFGeomBase<dim>::active_degree(unsigned int type) const {
      return active_fe_degrees[type];
    }
    template <unsigned int dim>
    const unsigned int DoFGeomBase<dim>::get_relevant_fe_index(
        const CurrentType &cur_type) const {
      auto temp = cur_type_to_fe_index.find(cur_type);
      #ifdef DEBUG
      if (cur_type_to_fe_index.find(cur_type) == cur_type_to_fe_index.end())
        assert(false && "The CurrentType chosen does not exist in this map!");
      #endif

      return temp->second;
    }
    template <unsigned int dim>
    void DoFGeomBase<dim>::set_fe_index(const unsigned int &new_index,
                                  const CurrentType &cur_type)
    {
      auto index = get_relevant_fe_index(cur_type);
      if (index < new_index)
        fully_populated = false;
      this->fe_indices[index] = new_index;
    }

    // This function simply refines all current types
    template <unsigned int dim>
    void DoFGeomBase<dim>::set_fe_index(const unsigned int &new_index) {
      for (unsigned int i = 0; i < this->fe_indices.size(); ++i)
      {
        if (this->fe_indices[i] < new_index)
          fully_populated = false;
        this->fe_indices[i] = new_index;
      }
    }

    DROMON_NAMESPACE_CLOSE
#endif //DROMON_DOFGEOMBASE_H
