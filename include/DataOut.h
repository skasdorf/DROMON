//
// Created by Jake J. Harmon (jake.harmon@ieee.org) on 7/19/21.
//

#ifndef DROMON_DATAOUT_H
#define DROMON_DATAOUT_H

#include "config.h"
#include "mesh.h"
#include "DoFHandler.h"
#include "DoFBase.h"

#include <fstream>
#include <chrono>
#include <ctime>
unsigned int vtk_point_index_from_ijk(unsigned int i, unsigned int j, unsigned int order);
DROMON_NAMESPACE_OPEN

template <unsigned int dim, unsigned int spacedim, unsigned int celltype, class CoefficientType = std::complex<double>, class Real = double>
class DataOut {
public:
  void attach_mesh(Mesh<dim, spacedim, celltype> *mesh);
  void attach_dof_handler(DoFHandler<dim, spacedim, CoefficientType, Real> *dof_handler);

  void vtk_out();
  DataOut();
  DataOut(Mesh<dim, spacedim, celltype> *mesh, std::string file_path, UpdateFlags update_flags = UpdateFlags::update_default);
  DataOut(Mesh<dim, spacedim, celltype> *mesh, DoFHandler<dim, spacedim, CoefficientType, Real> *dof_handler, std::string file_path, UpdateFlags update_flags = UpdateFlags::update_default);

  void add_cell_data(std::vector<double> &cell_data, std::string name = "default");
  void add_node_data(std::vector<double> &node_data, std::string name = "default");
  void set_flags(UpdateFlags update_flags);

private:
  Mesh<dim, spacedim, celltype> *mesh = nullptr;
  DoFHandler<dim, spacedim, CoefficientType, Real> *dof_handler = nullptr;

  std::string file_path;

  std::vector<std::vector<double>> cell_data_doubles;
  std::vector<std::vector<int>> cell_data_int;

  std::vector<std::string> cell_data_doubles_names;
  std::vector<std::string> cell_data_int_names;

  UpdateFlags update_flags;
};

    /**
     *
     * @tparam dim
     * @tparam spacedim
     * @tparam celltype
     *
     * It is important to note how VTK files treat higher order cells. For dim=2, for example, with a quadratic
     * cell, we have that
     * 3---------6----------2
     * |         |          |
     * |         |          |
     * 7---------8----------5
     * |         |          |
     * |         |          |
     * 0---------4----------1
     * This orientation holds for the biquadratic approximation, as well as the Lagrange cells.
     */
    template<unsigned int dim, unsigned int spacedim, unsigned int celltype, class CoefficientType, class Real>
    void DataOut<dim, spacedim, celltype, CoefficientType, Real>::vtk_out() {
        auto timenow = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
        std::ofstream vtk_out(file_path + ".vtk");
        vtk_out << "# vtk DataFile Version 3.0" << std::endl;
        //vtk_out << "# By Jake J. Harmon " << std::endl;
        vtk_out << "# File generated by the BrokkrIE Library (By Jake J. Harmon)";
        if (!((this->update_flags & UpdateFlags::suppress_comments) == UpdateFlags::suppress_comments))
          vtk_out << " on " << ctime(&timenow) << std::endl;
        else
          vtk_out << std::endl;

        vtk_out << "ASCII" << std::endl << "DATASET UNSTRUCTURED_GRID\n" << std::endl;

        vtk_out << "POINTS " << (std::round(std::pow(mesh->get_cell_type()+2, dim)))*this->mesh->n_cells() << " double\n";
        //Print out of the Points

        for (const auto& cell : mesh->cells)
        {
            const auto global_node_indices = cell.get_node_indices();
            for (const auto& node_index : global_node_indices)
            {
                vtk_out << *mesh->get_node(node_index) << std::endl;
            }
        }

        //Print out the number of cells
        vtk_out << "CELLS " << mesh->n_cells() << " " << (std::round(std::pow(mesh->get_cell_type()+2, dim))+1)*this->mesh->n_cells() << "\n" << std::endl;

        //Match the points with the cells
        unsigned int start = 0;
        for (const auto& cell : mesh->cells)
        {
            std::vector<unsigned int> connectivity(cell.n_nodes());
            unsigned int n = celltype + 1;
            switch (dim)
            {
                case 2:
                    for (unsigned int j= 0; j <= n; ++j)
                        for (unsigned int i = 0; i <= n; ++i)
                        {
                            unsigned int local_index = i + (n+1)*j;
                            unsigned int vtk_index = vtk_point_index_from_ijk(i,j,n);
                            connectivity[vtk_index] = local_index;
                        }
                    break;
            }
            vtk_out << cell.n_nodes();
            for (const auto index : connectivity)
                 vtk_out << "\t" << index + start;
            vtk_out << std::endl;
            start += cell.n_nodes();
        }

        //The cell types for the cells
        unsigned int vtk_cell_type;
        switch (mesh->get_cell_type())
        {
            case 0:
                switch (dim)
                {
                    case 1:
                        vtk_cell_type =3;
                        break;
                    case 2:
                        vtk_cell_type=9;
                        break;
                    case 3:
                        vtk_cell_type = 12;
                        break;
                }
                break;
            case 1:
                switch (dim)
                {
                    case 1:
                        vtk_cell_type = 21;
                        break;
                    case 2:
                        vtk_cell_type= 70; //28;
                        break;
                    case 3:
                        vtk_cell_type = 29;
                        break;
                }
                break;
            default:
                switch (dim)
                {
                    case 1:
                        vtk_cell_type = 68;
                        break;
                    case 2:
                        vtk_cell_type=70;
                        break;
                    case 3:
                        vtk_cell_type = 72;
                        break;
                }
                break;
        }
        vtk_out << "\nCELL_TYPES " << mesh->n_cells() << std::endl;
        for (unsigned int i = 0; i < mesh->n_cells(); ++i)
            vtk_out << " " << vtk_cell_type << std::endl;

        // Check if there is any scalar CELL_DATA and any POINTS_DATA to write to file...
        if (cell_data_doubles.size() > 0 || cell_data_int.size() > 0)
        {
          vtk_out << "CELL_DATA " << mesh->n_cells() << std::endl;

          //Get Scalar data first
          for (unsigned int i = 0; i < this->cell_data_doubles.size(); ++i)
          {
            vtk_out << "SCALARS " << this->cell_data_doubles_names[i] << " double 1" << std::endl;
            vtk_out << "LOOKUP_TABLE default" << std::endl;
            for (const auto& entry : this->cell_data_doubles[i])
              vtk_out << " " << entry;

            vtk_out << std::endl;
          }

        }

        vtk_out << std::endl;
        vtk_out.close();


        std::ofstream global_data_out(file_path + ".global.dat");
        global_data_out << "# Global Data File" << std::endl;
        global_data_out << "# File generated by the BrokkrIE Library (By Jake J. Harmon)";
        if (!((this->update_flags & UpdateFlags::suppress_comments) == UpdateFlags::suppress_comments))
          global_data_out << " on " << ctime(&timenow) << std::endl;
        else
          global_data_out << std::endl;
        // If we have any other update flags, then we have to print out to
        // the general mesh data

        if ((this->update_flags & UpdateFlags::verbose_output) == UpdateFlags::verbose_output)
        {
          // Output the General Statistics for the Mesh
          global_data_out << "# Number of cells (active)" << std::endl;
          global_data_out << mesh->n_cells() /* << " (" << mesh->n_active_cells() << ")" */ << std::endl;

          if (dof_handler != nullptr)
          {
            global_data_out << "# Number of DoFs (active)" << std::endl;
            global_data_out << dof_handler->get_n_dofs() << " (" << dof_handler->get_n_active_dofs() << ")" << std::endl;

            // Now print out the DoF info themselves
            for (auto& cell : dof_handler->get_dof_cells())
            {
              for (auto& dof : cell.get_dofs())
                global_data_out << dof << std::endl;
            }
            for (auto& face : dof_handler->get_dof_faces())
            {
              for (auto& dof : face.get_dofs())
                global_data_out << dof << std::endl;
            }
          }
        }

        if ((this->update_flags & UpdateFlags::output_active) == UpdateFlags::output_active)
        {
          // Output the General Statistics for the Mesh
          global_data_out << "# Number of cells (active)" << std::endl;
          global_data_out << mesh->n_cells() /* << " (" << mesh->n_active_cells() << ")" */ << std::endl;

          if (dof_handler != nullptr)
          {
            global_data_out << "# Number of DoFs (active)" << std::endl;
            global_data_out << dof_handler->get_n_dofs() << " (" << dof_handler->get_n_active_dofs() << ")" << std::endl;

            // Now print out the DoF info themselves
            for (auto& cell : dof_handler->get_dof_cells())
            {
              for (auto& dof : cell.get_dofs()) {
                if (dof.is_active)
                  global_data_out << dof << std::endl;
              }
            }
            for (auto& face : dof_handler->get_dof_faces())
            {
              for (auto& dof : face.get_dofs())
                if (dof.is_active)
                  global_data_out << dof << std::endl;
            }
          }
        }

        global_data_out.close();
    }

    template<unsigned int dim, unsigned int spacedim, unsigned int celltype, class CoefficientType, class Real>
    void DataOut<dim, spacedim, celltype, CoefficientType, Real>::attach_mesh(Mesh<dim, spacedim, celltype> *mesh) {
            this->mesh = mesh;
    }

    template<unsigned int dim, unsigned int spacedim, unsigned int celltype, class CoefficientType, class Real>
    void DataOut<dim, spacedim, celltype, CoefficientType, Real>::attach_dof_handler(DoFHandler<dim, spacedim, CoefficientType, Real> *dof_handler) {
      this->dof_handler = dof_handler;
    }

    template<unsigned int dim, unsigned int spacedim, unsigned int celltype, class CoefficientType, class Real>
    DataOut<dim, spacedim, celltype, CoefficientType, Real>::DataOut() {}

    template<unsigned int dim, unsigned int spacedim, unsigned int celltype, class CoefficientType, class Real>
    DataOut<dim, spacedim, celltype, CoefficientType, Real>::DataOut(Mesh<dim, spacedim, celltype> *mesh, std::string file_path, UpdateFlags update_flags) : update_flags(update_flags) {
            this->mesh = mesh;
            this->file_path = file_path;
    }

    template<unsigned int dim, unsigned int spacedim, unsigned int celltype, class CoefficientType, class Real>
    DataOut<dim, spacedim, celltype, CoefficientType, Real>::DataOut(Mesh<dim, spacedim, celltype> *mesh, DoFHandler<dim, spacedim, CoefficientType, Real> *dof_handler, std::string file_path, UpdateFlags update_flags) : update_flags(update_flags) {
      this->mesh = mesh;
      this->dof_handler = dof_handler;
      this->file_path = file_path;
    }

    template <unsigned int dim, unsigned int spacedim, unsigned int celltype, class CoefficientType, class Real>
    void
    DataOut<dim, spacedim, celltype, CoefficientType, Real>::add_cell_data(std::vector<double> &cell_data, std::string name) {
        assert(mesh != nullptr && "No mesh is attached to DataOut!");
        assert(cell_data.size() == mesh->n_cells() && "Size of cell_data does not match the number of cells!");

          cell_data_doubles.push_back(cell_data);
          if (name == "default")
          {
            name += cell_data_doubles.size()-1;
          }
          cell_data_doubles_names.push_back(name);

    }
    template <unsigned int dim, unsigned int spacedim, unsigned int celltype, class CoefficientType, class Real>
    void DataOut<dim, spacedim, celltype, CoefficientType, Real>::set_flags(UpdateFlags update_flags) {
      this->update_flags = update_flags;
    }

    DROMON_NAMESPACE_CLOSE
#endif //DROMON_DATAOUT_H