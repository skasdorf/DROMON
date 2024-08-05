#include "DROMON/DataOut.h"
#include "DROMON/DoFHandler.h"
#include "DROMON/EFIEIntegrator.h"
#include "DROMON/Excitations.h"
#include "DROMON/FECollection.h"
#include "DROMON/FE_HdivMaxOrtho.h"
#include "DROMON/GalerkinSystem.h"
#include "DROMON/IntegratorBase.h"
#include "DROMON/IteratorRanger.h"
#include "DROMON/Materials.h"
#include "DROMON/MatrixSolving.h"
#include "DROMON/MatrixSolvingPolicies.h"
#include "DROMON/MeshGenerator.h"
#include "DROMON/Point.h"
#include "DROMON/PostProcessing.h"
#include "DROMON/UniformIntegrator.h"
#include "DROMON/config.h"
#include "DROMON/mesh.h"
#include <iostream>
#include <map>
#include "DROMON/Refinement.h"
#include "DROMON/ErrorEstimation.h"
#include "DROMON/AdjointExcitations.h"

int main()
{  
  const unsigned int dim = 2;
  const unsigned int spacedim = 3;

  // dromon::Mesh<dim, spacedim, CUBICP> mesh;
  dromon::Mesh<dim, spacedim, LINEARP> mesh;

  std::vector<Point<spacedim, double>> points = {
      {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0., 1., 0.},
      {1., 1., 0.},    {2, 0., 0.},     {2., 1., 0.}};

    std::vector<std::array<
        unsigned int,
        dromon::GeometryInfo<dim, spacedim,
        mesh.cell_type>::nodes_per_cell>> cells = {{0, 1, 2, 3},
                 {1, 4, 3, 5}};


  dromon::MeshGenerator::create_mesh(mesh, points, cells);

  Point<dim, double> uv_1 = {-1, -1};
  auto jacobian_1 = mesh.cells[0].jacobian(uv_1);
  Point<dim, double> uv_2 = {1, -1};
  auto jacobian_2 = mesh.cells[0].jacobian(uv_1);
  Point<dim, double> uv_3 = {-1, 1};
  auto jacobian_3 = mesh.cells[0].jacobian(uv_1);
  Point<dim, double> uv_4 = {1, 1};
  auto jacobian_4 = mesh.cells[0].jacobian(uv_1);

  auto r1 = mesh.cells[0].r(uv_1);
  auto r2 = mesh.cells[0].r(uv_2);

  auto au = mesh.cells[0].drdu(uv_1);
  auto av = mesh.cells[0].drdv(uv_1);

  // Check the jacobian values
  const double reference_jacobian = 1./4.;
  if (jacobian_1 != reference_jacobian || jacobian_2 != reference_jacobian || jacobian_3 != reference_jacobian || jacobian_4 != reference_jacobian)
  {
  	std::cout << "Error: Jacobian values do not match the reference!" << std::endl;
  	return 1;
  }
  // Check the r values
  const Point<spacedim, double> r1_reference = points[0];
  const Point<spacedim, double> r2_reference = points[1];
  if (r1 != r1_reference || r2 != r2_reference)
  {
  	std::cout << "Error: r1 or r2 values do not match the reference!" << std::endl;
  	return 1;
  }
  const Point<spacedim, double> au_reference = {0.5,0.0,0.0};
  const Point<spacedim, double> av_reference = {0.0,0.5,0.0};
  if (au != au_reference || av != av_reference)
  {
  	std::cout << "Error: au or av values do not match the reference!" << std::endl;
  	return 1;
  }
  return 0;
 }
