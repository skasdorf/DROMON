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

  dromon::FECollectionCollector<dim, spacedim, double>
      fe_collection_collector;
  // Create a signle FECollection for the Electric Currents...
  dromon::FECollection<dim, spacedim, double> fe_collection_EFIE;

  const unsigned int starting_expansion_order = 1;
  const unsigned int ending_expansion_order = 3;

  for (unsigned exps = starting_expansion_order; exps <= ending_expansion_order;
       ++exps) {

    auto current_fe = dromon::FE_HdivMaxOrtho<dim, spacedim, double>(
        exps, dromon::CurrentType::Electric);
    fe_collection_EFIE.push_back(&current_fe);
  }

  fe_collection_collector.push_back(fe_collection_EFIE);

  dromon::DoFHandler<dim, spacedim, double> dof_handler(mesh);
  dof_handler.distribute_dofs(&fe_collection_collector);

  dromon::DataOut<dim, spacedim, LINEARP, double, double> data_out(
      &mesh, "test_mesh_generation", dromon::verbose_output | dromon::suppress_comments);
  data_out.attach_dof_handler(&dof_handler);

  std::vector<double> fe_degree_vector, n_dofs_per_cell;

  for (const auto &dof_cell : dof_handler.get_dof_cells()) {
    fe_degree_vector.push_back(dof_cell.active_degree(0));
    n_dofs_per_cell.push_back(
        dof_handler.get_n_active_dofs_on_cell(dof_cell.index));
  }
  data_out.add_cell_data(fe_degree_vector, "fe_degrees");
  data_out.add_cell_data(n_dofs_per_cell, "dofs_per_cell");

  data_out.vtk_out();
  
  return 0;
 }
