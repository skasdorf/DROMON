//
// Created by Jake J. Harmon (jake.harmon@ieee.org) on 7/18/21.
// Edited by Christopher A. Erickson (Christopher.Erickson@ieee.org) on 7/3/25
//

#ifndef DROMON_MESHGENERATOR_H
#define DROMON_MESHGENERATOR_H

#include "config.h"
#include "mesh.h"
#include "Point.h"
/**
 * @brief Utility function to compare whether two points have the same coordinates within a tolerance.
 * 
 * @tparam dim Dimension of the point.
 * @tparam Real Numeric type.
 * @param p1 First point.
 * @param p2 Second point.
 * @return True if points are considered equal, false otherwise.
 */
template<unsigned int dim, class Real>
bool is_same_coord(Point<dim, Real> p1, Point<dim, Real> p2) {
        if ((p1-p2).norm() < 1e-12)
            return true;
    return false;
}
DROMON_NAMESPACE_OPEN

    namespace MeshGenerator {
        /**
         * @brief Generates a hypercube mesh (line, square, or cube) subdivided into cells.
         *
         * @tparam dim Topological dimension (1D, 2D, 3D).
         * @tparam spacedim Embedding space dimension.
         * @tparam surf_type Type of surface (e.g., number of nodes per edge).
         * @param mesh Mesh object to be filled.
         * @param center Center point of the hypercube.
         * @param sidelength Length of the hypercube side.
         * @param n_cells_per_dim Number of subdivisions per dimension.
         *
         * Generates a hypercube, i.e.,
         * - a line for dim=1,
         * - a square for dim=2,
         * - a cube for dim=3.
         * Generates a hypercube, i.e., a line in dim=1, a square in dim=2, and a cube in dim = 3.
         * The size of the hypercube is described by the tensor product internal $$[left, right]^{dim}$$.
         */
        template<unsigned int dim, unsigned int spacedim, unsigned int surf_type>
        void hyper_cube(Mesh<dim, spacedim, surf_type> &mesh, 
                        Point<spacedim, double> center,
                        double sidelength, 
                        unsigned int n_cells_per_dim);
        /**
         * @brief Generates a NASA Almond geometry mesh using parametric equations.
         *
         * This function constructs a mesh representing the NASA Almond test shape,
         * which is widely used in computational electromagnetics and aerodynamics benchmarks.
         * The shape is parameterized by two coordinates:
         * - Longitudinal coordinate t ∈ [-0.4167, 0.5833]
         * - Angular coordinate ψ ∈ [-π, π]
         *
         * The mesh is generated as a structured grid of points computed from the NASA Almond
         * parametric equations:
         * 
         * For t < 0:
         * \f[
         * x = d \cdot t,\quad
         * y = 0.13333 \cdot d \cdot \sqrt{1 - \left( \frac{t}{0.416667} \right)^2} \cos\psi,\quad
         * z = 0.064444 \cdot d \cdot \sqrt{1 - \left( \frac{t}{0.416667} \right)^2} \sin\psi
         * \f]
         *
         * For t ≥ 0:
         * \f[
         * x = d \cdot t,\quad
         * y = 4.83345 \cdot d \cdot \left[ \sqrt{1 - \left( \frac{t}{2.08335} \right)^2} - 0.96 \right] \cos\psi,\quad
         * z = 1.9115 \cdot d \cdot \left[ \sqrt{1 - \left( \frac{t}{2.08335} \right)^2} - 0.96 \right] \sin\psi
         * \f]
         *
         * @tparam dim       Topological dimension of the mesh (e.g., 2 for a surface mesh).
         * @tparam spacedim  Embedding space dimension (e.g., 3 for 3D space).
         * @tparam surf_type Type of surface representation (related to polynomial order).
         *
         * @param mesh    Mesh object to populate with nodes and cells.
         * @param center  Center point of the almond geometry.
         * @param d       Scaling parameter, default is 9.936 inches (NASA standard).
         * @param n_t     Number of subdivisions along the longitudinal axis (default: 50).
         * @param n_psi   Number of subdivisions along the angular direction (default: 100).
         *
         * The generated mesh consists of nodes laid out on a structured parametric grid
         * and quadrilateral cells formed by connecting neighboring points.
         *
         * @note This function is primarily intended for generating the NASA Almond benchmark geometry
         * for validation of computational electromagnetics or CFD solvers.
         *
         * @see create_mesh(), hyper_sphere(), hyper_cube()
         *
         * @warning Ensure dim=2 and spacedim=3 for correct geometry. Other configurations are not supported.
         *
         * **Example usage:**
         * @code
         * Mesh<2, 3, CUBICP> mesh;
         * Point<3, double> center = {0.0, 0.0, 0.0};
         * MeshGenerator::nasa_almond(mesh, center);
         * @endcode
         */
        template<unsigned int dim, unsigned int spacedim, unsigned int surf_type>
        void nasa_almond(Mesh<dim, spacedim, surf_type> &mesh,
                        Point<spacedim, double> center,
                        double d = 9.936,
                        unsigned int N_samples = 10);
        /**
         * @brief Implementation of nasa_almond.
         * @see nasa_almond()
         */
        template<unsigned int dim, unsigned int spacedim, unsigned int surf_type>
        void nasa_almond(Mesh<dim, spacedim, surf_type>& mesh,
                        Point<spacedim, double> center,
                        double d,
                        unsigned int N_samples)
        {
            // ---------- Helper Lambdas ----------
            // logspace: Generate a vector of logarithmically spaced values between a and b
            auto logspace = [](double a, double b, unsigned int n) {
                std::vector<double> v(n);
                double log_a = std::log10(a), log_b = std::log10(b);
                double step = (log_b - log_a) / (n - 1);
                for (unsigned int i = 0; i < n; ++i)
                    v[i] = std::pow(10, log_a + i * step);
                return v;
            };

            // rescale: Rescale an array to a new range [new_min, new_max]
            auto rescale = [](std::vector<double>& arr, double new_min, double new_max) {
                double old_min = *std::min_element(arr.begin(), arr.end());
                double old_max = *std::max_element(arr.begin(), arr.end());
                for (auto& x : arr)
                    x = new_min + (x - old_min) * (new_max - new_min) / (old_max - old_min);
            };

            // ---------- Generate Back Section ----------
            // The back (negative x-direction) is approximated from parameter a to tip_back
            double a = -0.41667; // Starting parametric coordinate for the back
            double tip_back = a - a / (5 * N_samples + 1); // Slightly beyond a for tip refinement
            std::vector<double> t_back = {a, a + (tip_back - a) / 3, a + 2 * (tip_back - a) / 3, tip_back};

            // Generate logarithmic spacing for refinement
            auto temp = logspace(1, 1.9, N_samples);
            rescale(temp, tip_back, 0); // Map values into [tip_back, 0]
            std::vector<double> t_new = {tip_back};
            // Insert intermediate values using sub-interval divisions
            for (unsigned int i = 1; i < temp.size(); ++i) {
                double hs = (temp[i] - temp[i - 1]) / 3;
                for (int j = 0; j < 3; ++j)
                    t_new.push_back(t_new.back() + hs);
            }
            t_back.insert(t_back.end(), t_new.begin(), t_new.end());
            std::sort(t_back.begin(), t_back.end());

            // ---------- Generate Tip Section ----------
            // The front (positive x-direction) is parameterized similarly
            double b = 0.58333; // Tip end coordinate
            double tip_end = b - b / (3 * N_samples + 1);
            auto temp_tip = logspace(1, 1.47, N_samples);
            rescale(temp_tip, 0, tip_end);
            std::vector<double> t_tip = {0};
            for (unsigned int i = 1; i < temp_tip.size(); ++i) {
                double hs = (temp_tip[i] - temp_tip[i - 1]) / 3;
                for (int j = 0; j < 3; ++j)
                    t_tip.push_back(t_tip.back() + hs);
            }
            t_tip.push_back(tip_end);
            t_tip.push_back(b);

            // Combine both sections (back + tip)
            std::vector<double> t_values = t_back;
            t_values.insert(t_values.end(), t_tip.begin(), t_tip.end());

            // Sort and remove duplicates within a tolerance (1e-12)
            const double tol = 1e-12;
            std::sort(t_values.begin(), t_values.end());
            auto new_end = std::unique(t_values.begin(), t_values.end(),
                [tol](double x, double y){ return std::abs(x-y) < tol; });
            t_values.erase(new_end, t_values.end());

            // ---------- Generate Points ----------
            std::vector<std::vector<unsigned int>> rings;  // Stores index rings
            std::vector<Point<spacedim, double>> points;   // Stores generated points

            // Loop through each longitudinal section
            for (unsigned int i = 0; i < t_values.size(); ++i) {
                std::vector<double> psi; // Azimuthal angles for each ring
                // Determine number of points in ring based on position
                if (i == 0 || i == t_values.size() - 1)
                    psi = {M_PI / 2}; // Single point at poles
                else if (i == 1 || i == t_values.size() - 2)
                    for (int k = 0; k < 9; ++k) psi.push_back(2 * M_PI * k / 8);
                else if (i == 2 || i == t_values.size() - 3)
                    for (int k = 0; k < 17; ++k) psi.push_back(2 * M_PI * k / 16);
                else
                    for (int k = 0; k < 25; ++k) psi.push_back(2 * M_PI * k / 24);

                std::vector<unsigned int> ring_indices;
                double t = t_values[i];
                std::cout << t << std::endl;

                // Generate coordinates for each angle in the ring
                for (unsigned int j = 0; j < psi.size(); ++j) {
                    double x = d * t;
                    double y, z;
                    if (t < 0) {
                        // Back section profile equation
                        double ratio = t / 0.416667;
                        double factor = std::sqrt(std::max(0.0, 1.0 - ratio * ratio));
                        y = 0.193333 * d * factor * std::cos(psi[j]);
                        z = 0.064444 * d * factor * std::sin(psi[j]);
                    } else {
                        // Front section profile equation
                        double ratio = t / 2.08335;
                        double factor = std::sqrt(std::max(0.0, 1.0 - ratio * ratio));
                        y = 4.83345 * d * (factor - 0.96) * std::cos(psi[j]);
                        z = 1.61115 * d * (factor - 0.96) * std::sin(psi[j]);
                    }

                    // Apply center offset and store the point
                    Point<spacedim, double> pt = {x + center[0], y + center[1], z + center[2]};
                    points.push_back(pt);
                    ring_indices.push_back(points.size() - 1);
                }
                rings.push_back(ring_indices);
            }

            // ---------- Add Nodes ----------
            for (auto& pt : points) {
                Node<dim, spacedim, double> node(&mesh, pt);
                mesh.push_back_node(&node);
            }

            // ---------- Generate Connectivity (quads only) ----------
            // Create quadrilateral cells connecting successive rings
            for (std::size_t i = 0; i + 1 < rings.size(); ++i) {
                const auto &R1 = rings[i];
                const auto &R2 = rings[i+1];
                std::size_t m1 = R1.size();
                std::size_t m2 = R2.size();
                std::size_t max_m = std::max(m1, m2);

                for (std::size_t k = 0; k < max_m; ++k) {
                    unsigned i0 = R1[k % m1];
                    unsigned i1 = R1[(k + 1) % m1];
                    unsigned i2 = R2[(k + 1) % m2];
                    unsigned i3 = R2[k % m2];

                    // Create cell connectivity array
                    std::array<unsigned,
                        GeometryInfo<dim, spacedim, surf_type>::nodes_per_cell> cell;

                    // Assign corner nodes
                    cell[0] = i0;
                    cell[1] = i1;
                    cell[2] = i2;
                    cell[3] = i3;

                    // Fill remaining nodes (for high-order) with duplicates
                    for (unsigned idx = 4; idx < cell.size(); ++idx)
                        cell[idx] = i0;

                    mesh.insert_cell(cell);
                }
            }

            // Generate faces and finalize the mesh
            mesh.spawn_and_assign_faces();
            std::cout << "NASA Almond mesh generated with " << mesh.n_nodes()
                    << " nodes and " << mesh.n_cells() << " cells." << std::endl;
        }


        /**
         * @brief Implementation of hyper_cube.
         * @see hyper_cube()
         */
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
            std::cout << "Hyper Cube mesh generated with " << mesh.n_nodes() << " nodes and " << mesh.n_cells() << " cells." << std::endl;
            mesh.spawn_and_assign_faces();
            //mesh.finalize_geometry(); //This function should make the normals of each patch consistent, etc.

        }
        /**
         * @brief Generates a hypersphere mesh by projecting a hypercube onto the sphere surface.
         *
         * @tparam dim Topological dimension.
         * @tparam spacedim Embedding space dimension.
         * @tparam surf_type Type of surface.
         * @param mesh Mesh object to be filled.
         * @param center Center point of the sphere.
         * @param sidelength Radius of the sphere.
         * @param n_subdiv Number of subdivisions per dimension.
         *
         * The hypersphere is approximated by projecting the hypercube grid onto a sphere.
         */
        template<unsigned int dim, unsigned int spacedim, unsigned int surf_type>
        void hyper_sphere(Mesh<dim, spacedim, surf_type> &mesh, Point<spacedim, double> center, double sidelength, unsigned int n_subdiv);
        /**
         * @brief Implementation of hyper_sphere.
         * @see hyper_sphere()
         */
        template<unsigned int dim, unsigned int spacedim, unsigned int surf_type>
        void
        dromon::MeshGenerator::hyper_sphere(Mesh<dim, spacedim, surf_type> &mesh, Point<spacedim, double> center, double sidelength, unsigned int n_subdiv)
        {   //CAE 7/3/25 I dont know if I like this kind of implmentation, why not an isosphere?
            //first, generate a hypercube
            hyper_cube(mesh, center, sidelength, n_subdiv);
            //now, loop through every vertex and project them such that the radius is sidelength
            for (auto& node : mesh.nodes)
            {
                auto& pt = node->get_point();
                const auto length = pt.norm()/sidelength;
                pt/=length;

            }
            std::cout << "Hyper Sphere mesh generated with " << mesh.n_nodes() << " nodes and " << mesh.n_cells() << " cells." << std::endl;
        }
        /**
         * @brief Generates a flat square plate mesh.
         *
         * @tparam dim Topological dimension.
         * @tparam spacedim Embedding space dimension.
         * @tparam surf_type Type of surface.
         * @param mesh Mesh object to be filled.
         * @param center Center point of the plate.
         * @param sidelength Length of the plate side.
         * @param n_cells_per_dim Number of subdivisions per dimension.
         */
        template<unsigned int dim, unsigned int spacedim, unsigned int surf_type>
        void square_plate(Mesh<dim, spacedim, surf_type> &mesh, Point<spacedim, double> center, double sidelength, unsigned int n_cells_per_dim);
        /**
         * @brief Implementation of square_plate.
         * @see square_plate()
         */
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
          std::cout << "Square Plate mesh generated with " << mesh.n_nodes() << " nodes and " << mesh.n_cells() << " cells." << std::endl;
          mesh.spawn_and_assign_faces();
          //mesh.finalize_geometry(); //This function should make the normals of each patch consistent, etc.
        }
        /**
         * @brief Creates a mesh by explicitly providing node coordinates and cell connectivity.
         *
         * @tparam dim Topological dimension.
         * @tparam spacedim Embedding space dimension.
         * @tparam surf_type Type of surface.
         * @param mesh Mesh object to be filled.
         * @param points List of node coordinates.
         * @param cells List of cells, each defined by a set of node indices.
         *
         * This function is useful for custom meshes built from arbitrary point and cell definitions.
         */
        template<unsigned int dim, unsigned int spacedim, unsigned int surf_type>
        void create_mesh(
            Mesh<dim, spacedim, surf_type>& mesh,
            std::vector<Point<spacedim, double>> points,
            std::vector<std::array<unsigned int, GeometryInfo<dim, spacedim, surf_type>::nodes_per_cell>> cells);

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

            std::cout << "Mesh generated with " << mesh.n_nodes() << " nodes and " << mesh.n_cells() << " cells." << std::endl;
            mesh.spawn_and_assign_faces();
        }
    }

DROMON_NAMESPACE_CLOSE

#endif //DROMON_MESHGENERATOR_H
