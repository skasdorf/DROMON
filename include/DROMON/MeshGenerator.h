//
// Created by Jake J. Harmon (jake.harmon@ieee.org) on 7/18/21.
//

#ifndef DROMON_MESHGENERATOR_H
#define DROMON_MESHGENERATOR_H

#include "config.h"
#include "mesh.h"
#include "Point.h"
template<unsigned int dim, class Real>
bool is_same_coord(Point<dim, Real> p1, Point<dim, Real> p2) {
        if ((p1-p2).norm() < 1e-12)
            return true;
    return false;
}
DROMON_NAMESPACE_OPEN

    namespace MeshGenerator {
        /**
         *
         * @tparam dim
         * @tparam spacedim
         * @tparam Real
         * @param mesh
         * @param left
         * @param right
         * Generates a hypercube, i.e., a line in dim=1, a square in dim=2, and a cube in dim = 3.
         * The size of the hypercube is described by the tensor product internal $[left, right]^{dim}$.
         */
        template<unsigned int dim, unsigned int spacedim, unsigned int surf_type>
        void hyper_cube(Mesh<dim, spacedim, surf_type> &mesh, Point<spacedim, double> center, double sidelength, unsigned int n_cells_per_dim);

        template<unsigned int dim, unsigned int spacedim, unsigned int surf_type>
        void
        dromon::MeshGenerator::hyper_cube(Mesh<dim, spacedim, surf_type> &mesh, Point<spacedim, double> center, double sidelength,
                                            unsigned int n_cells_per_dim) {
//            if (mesh_order > 0)
//                assert(mesh_order + 2 + ((n_subdiv + 1) / (mesh_order + 1) - 1) * (mesh_order + 1) == n_subdiv + 1
//                    && "Insufficient n_subdiv for the desired mesh order!");
            unsigned int n_subdiv = n_cells_per_dim*(surf_type+1);
            // We first build a reference cube centered at the origin of size [-1,1]^dim.
            // The hypercube is then scaled and translated to match the desired $center$ and $sidelength$.
            std::vector<double> tensor_points;
            for (unsigned int i = 0; i <= n_subdiv; ++i)
                tensor_points.push_back((double) -1 + (double) i * ((double) 2) / ((double) n_subdiv));

            switch (dim) {
                case 1:
                    //In this case, nothing needs to be done except divided into cells, so we first add all of the points
                    for (unsigned int i = 0; i < tensor_points.size(); ++i) {
                        Point<spacedim, double> temp_point;
                        temp_point[0] = tensor_points[i];
                        Node<dim, spacedim, double> temp_node(&mesh, temp_point);
                        mesh.push_back_node(&temp_node);
                    }
                    for (unsigned int cell_index = 0; cell_index < n_subdiv / (surf_type + 1); ++cell_index) {
                        std::vector<int> temp_vertex_indices;
                        for (unsigned int i = cell_index * (surf_type + 1);
                             i <= cell_index * (surf_type + 1) + surf_type + 1; ++i)
                            temp_vertex_indices.push_back(i);
                        //TODO: Fix cell
//                        Cell<dim, dim, spacedim> temp_cell(&mesh, temp_vertex_indices);
//                        mesh.push_back_cell(&temp_cell);
                    }
                    break;
                case 2:
                    //In this case, we build the shell of the hypercube
                    if (spacedim == 3) {
                        //fix x
                        std::vector<unsigned int> face1, face2, face3, face4, face5, face6;
                        // std::vector<Point_3> coords;
                        //fix x
                        int i = 0;
                        bool skip = false;
                        for (int k = 0; k < tensor_points.size(); ++k) {
                            for (int j = tensor_points.size() -1 ; j >= 0; --j) {
                                skip = false;
                                Point<spacedim, double> temp = {tensor_points[i], tensor_points[j], tensor_points[k]};
                                for (int l = 0; l < mesh.n_nodes(); ++l) {

                                    if (is_same_coord(temp, static_cast<Node<dim, spacedim, double>*>(mesh.get_node(l))->get_point())) {
                                        face1.push_back(l);
                                        skip = true;
                                        break;
                                    }
                                }
                                if (!skip) {
                                    Node<dim, spacedim, double> temp_node(&mesh, temp);
                                    mesh.push_back_node(&temp_node);
                                    face1.push_back(mesh.n_nodes() - 1);

                                }
                            }
                        }
                        i = tensor_points.size() - 1;
                        skip = false;
                        for (int k = 0; k < tensor_points.size(); ++k) {
                            for (int j = 0; j < tensor_points.size(); ++j) {
                                skip = false;
                                Point<spacedim, double> temp = {tensor_points[i], tensor_points[j], tensor_points[k]};
                                for (int l = 0; l < mesh.n_nodes(); ++l) {
                                    if (is_same_coord(temp, static_cast<Node<dim, spacedim, double>*>(mesh.get_node(l))->get_point())) {
                                        face2.push_back(l);
                                        skip = true;
                                        break;
                                    }

                                }
                                if (!skip) {
                                    Node<dim, spacedim, double> temp_node(&mesh, temp);
                                    mesh.push_back_node(&temp_node);
                                    face2.push_back(mesh.n_nodes() - 1);

                                }
                            }
                        }
                        //fix y
                        int j = 0;
                        for (int k = 0; k < tensor_points.size(); ++k){
                             for (i = 0; i < tensor_points.size(); ++i) {
                                skip = false;
                                Point<spacedim, double> temp = {tensor_points[i], tensor_points[j], tensor_points[k]};
                                for (int l = 0; l < mesh.n_nodes(); ++l) {
                                    if (is_same_coord(temp, static_cast<Node<dim, spacedim, double>*>(mesh.get_node(l))->get_point())) {
                                        face3.push_back(l);
                                        skip = true;
                                        break;
                                    }

                                }
                                if (!skip) {
                                    Node<dim, spacedim, double> temp_node(&mesh, temp);
                                    mesh.push_back_node(&temp_node);
                                    face3.push_back(mesh.n_nodes() - 1);

                                }
                            }
                        }

                        j = tensor_points.size() - 1;
                        for (int k = 0; k < tensor_points.size(); ++k) {
                            for (i = tensor_points.size()-1; i >= 0; --i) {
                                skip = false;
                                Point<spacedim, double> temp = {tensor_points[i], tensor_points[j], tensor_points[k]};
                                for (int l = 0; l < mesh.n_nodes(); ++l) {
                                    if (is_same_coord(temp, static_cast<Node<dim, spacedim, double>*>(mesh.get_node(l))->get_point())) {
                                        face4.push_back(l);
                                        skip = true;
                                        break;
                                    }

                                }
                                if (!skip) {
                                    Node<dim, spacedim, double> temp_node(&mesh, temp);
                                    mesh.push_back_node(&temp_node);
                                    face4.push_back(mesh.n_nodes() - 1);

                                }
                            }
                        }

                        int k = 0;
                        for (j = tensor_points.size()-1; j >= 0; --j) {
                             for (i = 0; i < tensor_points.size(); ++i) {
                                skip = false;
                                Point<spacedim, double> temp = {tensor_points[i], tensor_points[j], tensor_points[k]};
                                for (int l = 0; l < mesh.n_nodes(); ++l) {
                                    if (is_same_coord(temp, static_cast<Node<dim, spacedim, double>*>(mesh.get_node(l))->get_point())) {
                                        face5.push_back(l);
                                        skip = true;
                                        break;
                                    }

                                }
                                if (!skip) {
                                    Node<dim, spacedim, double> temp_node(&mesh, temp);
                                    mesh.push_back_node(&temp_node);
                                    face5.push_back(mesh.n_nodes() - 1);

                                }
                            }
                        }

                        k = tensor_points.size() - 1;
                        for (j = 0; j < tensor_points.size(); ++j){
                             for (i = 0; i < tensor_points.size(); ++i) {
                                skip = false;
                                Point<spacedim, double> temp = {tensor_points[i], tensor_points[j], tensor_points[k]};
                                for (int l = 0; l < mesh.n_nodes(); ++l) {
                                    if (is_same_coord(temp, static_cast<Node<dim, spacedim, double>*>(mesh.get_node(l))->get_point())) {
                                        face6.push_back(l);
                                        skip = true;
                                    }

                                }
                                if (!skip) {
                                    Node<dim, spacedim, double> temp_node(&mesh, temp);
                                    mesh.push_back_node(&temp_node);
                                    face6.push_back(mesh.n_nodes() - 1);

                                }
                            }
                        }
                        //Now insert the faces
                        mesh.cells_from_grid(face1,n_subdiv + 1);
                        mesh.cells_from_grid(face2,n_subdiv + 1);
                        mesh.cells_from_grid(face3,n_subdiv + 1);
                        mesh.cells_from_grid(face4,n_subdiv + 1);
                        mesh.cells_from_grid(face5,n_subdiv + 1);
                        mesh.cells_from_grid(face6,n_subdiv + 1);
                    }

                    break;
            }
            std::cout << "Done!" << std::endl;
            mesh.spawn_and_assign_faces();
            //mesh.finalize_geometry(); //This function should make the normals of each patch consistent, etc.

        }

        template<unsigned int dim, unsigned int spacedim, unsigned int surf_type>
        void hyper_sphere(Mesh<dim, spacedim, surf_type> &mesh, Point<spacedim, double> center, double sidelength, unsigned int n_subdiv);

        template<unsigned int dim, unsigned int spacedim, unsigned int surf_type>
        void
        dromon::MeshGenerator::hyper_sphere(Mesh<dim, spacedim, surf_type> &mesh, Point<spacedim, double> center, double sidelength, unsigned int n_subdiv)
        {
            //first, generate a hypercube
            hyper_cube(mesh, center, sidelength, n_subdiv);
            //now, loop through every vertex and project them such that the radius is sidelength
            for (auto& node : mesh.nodes)
            {
                auto& pt = node->get_point();
                const auto length = pt.norm()/sidelength;
                pt/=length;

            }
        }

        template<unsigned int dim, unsigned int spacedim, unsigned int surf_type>
        void square_plate(Mesh<dim, spacedim, surf_type> &mesh, Point<spacedim, double> center, double sidelength, unsigned int n_cells_per_dim);

        template<unsigned int dim, unsigned int spacedim, unsigned int surf_type>
        void
        dromon::MeshGenerator::square_plate(Mesh<dim, spacedim, surf_type> &mesh, Point<spacedim, double> center, double sidelength, unsigned int n_cells_per_dim)
        {

          // We first build a reference cube centered at the origin of size [-1,1]^dim.
          // The hypercube is then scaled and translated to match the desired $center$ and $sidelength$.
          const unsigned int n_subdiv = n_cells_per_dim*(surf_type+1);
          std::vector<double> tensor_points;
          for (unsigned int i = 0; i <= n_subdiv; ++i)
            tensor_points.push_back((double) -1 + (double) i * ((double) 2) / ((double) n_subdiv));

          for (auto& tensor_point : tensor_points)
            tensor_point *= sidelength/double(2);

          switch (dim) {
          case 1:
            //In this case, nothing needs to be done except divided into cells, so we first add all of the points
            for (unsigned int i = 0; i < tensor_points.size(); ++i) {
              Point<spacedim, double> temp_point;
              temp_point[0] = tensor_points[i];
              Node<dim, spacedim, double> temp_node(&mesh, temp_point);
              mesh.push_back_node(&temp_node);
            }
            for (unsigned int cell_index = 0; cell_index < n_subdiv / (surf_type + 1); ++cell_index) {
              std::vector<int> temp_vertex_indices;
              for (unsigned int i = cell_index * (surf_type + 1);
                   i <= cell_index * (surf_type + 1) + surf_type + 1; ++i)
                temp_vertex_indices.push_back(i);
              //TODO: Fix cell
//                        Cell<dim, dim, spacedim> temp_cell(&mesh, temp_vertex_indices);
//                        mesh.push_back_cell(&temp_cell);
            }
            break;
          case 2:
            //In this case, we build the shell of the hypercube
            if (spacedim == 3) {
              //fix x
              std::vector<unsigned int> face1;
              // std::vector<Point_3> coords;
              //fix x
              int i = 0;
              bool skip = false;
              for (int k = 0; k < tensor_points.size(); ++k) {
                for (int j = 0 ; j < tensor_points.size(); ++j) {
                  skip = false;
                  Point<spacedim, double> temp = {tensor_points[i], tensor_points[j], tensor_points[k]};
                  for (int l = 0; l < mesh.n_nodes(); ++l) {

                    if (is_same_coord(temp, static_cast<Node<dim, spacedim, double>*>(mesh.get_node(l))->get_point())) {
                      face1.push_back(l);
                      skip = true;
                      break;
                    }
                  }
                  if (!skip) {
                    Node<dim, spacedim, double> temp_node(&mesh, temp);
                    mesh.push_back_node(&temp_node);
                    face1.push_back(mesh.n_nodes() - 1);

                  }
                }
              }


              //Now insert the faces
              mesh.cells_from_grid(face1,n_subdiv + 1);

            }

            break;
          }
          std::cout << "Done!" << std::endl;
          mesh.spawn_and_assign_faces();
          //mesh.finalize_geometry(); //This function should make the normals of each patch consistent, etc.
        }

        template<unsigned int dim, unsigned int spacedim, unsigned int surf_type>
        void create_mesh(Mesh<dim, spacedim, surf_type> &mesh, std::vector<Point<spacedim, double>> points, std::vector<std::array<unsigned int, GeometryInfo<dim, spacedim, surf_type>::nodes_per_cell>> cells)
        {
            // Loop through all nodes and add them to the mesh...
            for (auto& temp : points)
            {
              Node<dim, spacedim, double> temp_node(&mesh, temp);
              mesh.push_back_node(&temp_node);
            }

            // Add all the cells
            for (auto& cell : cells)
              mesh.insert_cell(cell);

            std::cout << "Done!" << std::endl;
            mesh.spawn_and_assign_faces();
        }
    }

DROMON_NAMESPACE_CLOSE

#endif //DROMON_MESHGENERATOR_H
