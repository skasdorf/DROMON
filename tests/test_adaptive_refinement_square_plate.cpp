#include "DROMON/FE_HdivMaxOrtho.h"
#include "DROMON/Materials.h"
#include "DROMON/MeshGenerator.h"
#include "DROMON/Point.h"
#include "DROMON/config.h"
#include "DROMON/mesh.h"


#include <iostream>
#include "programs.h"



int main() {
  using namespace dromon;

  // Make the MaterialData
  // In this case, a PEC object embedded in air
  MaterialData<double> mat_dom;
  mat_dom.push_back(dromon::MaterialDomain<double>(
      dromon::Material<double>(1.0, 1.0, true, false),
      dromon::Material<double>(1.0, 1.0, false, false)));

  // Now let's define the mesh

  //  dromon::Mesh<2, 3, CUBICP> mesh;
  //  Point<3, double> center = {0.0, 0.0, 0.0};
  //  double sidelength = 1.0;
  //  unsigned int n_cells_per_dim = 1;
  //  dromon::MeshGenerator::hyper_sphere(mesh, center, sidelength,
  //                                        n_cells_per_dim);
  Mesh<2, 3, LINEARP> mesh;
  Point<3, double> center = {0.0, 0.0, 0.0};
    double sidelength = 2.0;
    unsigned int n_cells_per_dim = 5;
    MeshGenerator::square_plate(mesh, center, sidelength,
                                          n_cells_per_dim);
  // CAE: Not sure when this was written or by who, but the adaptive solver needs two meshes, so I just decided to use the same mesh twice. Probably should come back to this later and aactually look at how adaptive solver works and adjust the meshes as needed. Perhpas it is expecting the PEC and Air to be seperate meshes?
  Problems::AdaptiveSolver solver(&mesh, &mesh, &mat_dom, 1, 10 ,15,15,15,15);
  solver.set_plane_wave_excitation(300e6, {1.0,0.0});
  std::cout << std::setprecision(16);
  double theta_sc = constants<double>::PI / 2.0, phi_sc = 0.0;
  double R_dist_scalar = 200000.0;
  solver.set_scattering_parameters(theta_sc,phi_sc,R_dist_scalar,{0.0,0.0,1.0});
  // CAE: This does not compile, as it is asking for 11 inputs. specifically all the various solves of the matrices. I think this is all for the adjoint solving. Still new to the code so I gotta figure out how to get all those.
  // CAE: Testing out placeholders?
  std::vector<std::complex<double>> forward_matrix;
  std::vector<std::complex<double>> forward_excitation;
  std::vector<std::complex<double>> forward_solution;
  std::vector<std::complex<double>> adjoint_solution;
  std::complex<double> gradient_solution;
  // CAE: I also re wrote the refinement call below, to show the 11 inputs, and what they mean as it was very vague before.
  solver.execute_refinement(
    1,                  // adjoint_starting_index
    0.01,               // reltol
    forward_matrix,
    forward_excitation,
    forward_solution,
    adjoint_solution,
    gradient_solution,
    10.0,               // sidelength
    0,                  // HOPSflag
    0,                  // perturbVar
    0.01                // perturbSize
  );
  std::cout << "hello";
  return 0;
}

