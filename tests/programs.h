#include "DROMON/AdjointExcitations.h"
#include "DROMON/DataOut.h"
#include "DROMON/DoFHandler.h"
#include "DROMON/EFIEIntegrator.h"
#include "DROMON/ErrorEstimation.h"
#include "DROMON/Excitations.h"
#include "DROMON/FECollection.h"
#include "DROMON/FE_HdivMaxOrtho.h"
#include "DROMON/GalerkinSystem.h"
#include "DROMON/Materials.h"
#include "DROMON/MatrixSolving.h"
#include "DROMON/MatrixSolvingPolicies.h"
#include "DROMON/MeshGenerator.h"
#include "DROMON/Point.h"
#include "DROMON/PostProcessing.h"
#include "DROMON/Refinement.h"
#include "DROMON/UniformIntegrator.h"
#include "DROMON/config.h"
#include "DROMON/mesh.h"

#include <complex>
#include "mkl.h"
#include "mkl_lapack.h"
#include <chrono>

#include <iostream>
#include <map>

namespace Problems {
using namespace dromon;
template <unsigned int dim, unsigned int spacedim, unsigned int patch_order>
class RegularSolver {
public:
  RegularSolver(Mesh<dim, spacedim, patch_order> *mesh,
                MaterialData<double> *mat_dom,
                const unsigned int &min_expansion_order = 1,
                const unsigned int &max_expansion_order = 1,
                const unsigned int &ngl_regular = 5,
                const unsigned int &ngl_vertex = 5,
                const unsigned int &ngl_edge = 5,
                const unsigned int &ngl_self = 5)
      : min_expansion_order(min_expansion_order),
        max_expansion_order(max_expansion_order), ngl_regular(ngl_regular),
        ngl_vertex(ngl_vertex), ngl_edge(ngl_edge), ngl_self(ngl_self) {
    this->mesh = mesh;
    this->mat_dom = mat_dom;
    this->initialize_fe_collections();
    this->initialize_dof_handler();
  }
  void conduct_solution_step();

  Point<spacedim, std::complex<double>> post_process_scattered_field();
  double post_process_RCS();
  void initialize_DoF_systems(const unsigned int &initial_fe_index);
  // void initialize_DoF_systems(const std::vector<unsigned int>&
  // initial_fe_indices);

  void set_plane_wave_excitation(const double &frequency,
                                 const Point<dim, std::complex<double>> &E_mag);

private:
  void fill_system();
  void assemble_system_matrix();
  void solve_system();
  void initialize_fe_collections() {
    for (unsigned int exps = min_expansion_order; exps <= max_expansion_order;
         ++exps) {
      auto current_fe =
          FE_HdivMaxOrtho<dim, spacedim, double>(exps, CurrentType::Electric);
      fe_collection_EFIE.push_back(&current_fe);
    }
    fe_collection_collector.push_back(fe_collection_EFIE);
  }

  void initialize_dof_handler() {
    this->dof_handler = std::make_unique<DoFHandler<dim, spacedim>>(
        DoFHandler<dim, spacedim>(*mesh));
  }
  double freq;
  const unsigned int min_expansion_order;
  const unsigned int max_expansion_order;

  const unsigned int ngl_regular;
  const unsigned int ngl_self;
  const unsigned int ngl_vertex;
  const unsigned int ngl_edge;

  Mesh<dim, spacedim, patch_order> *mesh;
  MaterialData<double> *mat_dom;

  FECollectionCollector<dim, spacedim, double> fe_collection_collector;
  FECollection<dim, spacedim, double> fe_collection_EFIE;

  std::unique_ptr<DoFHandler<dim, spacedim>> dof_handler;
  std::unique_ptr<GalerkinSystem<dim, spacedim>> cg_sys;
  std::unique_ptr<GalerkinSystem<dim, spacedim>> cg_sys_HOPS;


  const double theta_inc = constants<double>::PI / 2.0;
  const double phi_inc = 0.0;

  std::unique_ptr<Excitations::PlaneWave<spacedim, double>> excitation;

  EFIEIntegrator<DoFParent<dim, spacedim>, std::complex<double>> integrator;
  std::unique_ptr<MatrixSolving::Matrix<
      std::complex<double>, dromon::MatrixSolvingPolicies::MKLPolicy>>
      matrix;
  std::vector<std::complex<double>> current_forward_solution;

  //  std::unique_ptr<MatrixSolving::Matrix<
  //      std::complex<double>, dromon::MatrixSolvingPolicies::MKLPolicy>>
  //      matrix_hi;
  //  std::unique_ptr<MatrixSolving::Matrix<
  //      std::complex<double>, dromon::MatrixSolvingPolicies::MKLPolicy>>
  //      matrix_lo;
};
template <unsigned int dim, unsigned int spacedim, unsigned int patch_order>
void RegularSolver<dim, spacedim, patch_order>::initialize_DoF_systems(
    const unsigned int &initial_fe_index) {
  this->dof_handler->distribute_dofs(&fe_collection_collector,
                                     initial_fe_index);
  this->cg_sys = std::make_unique<GalerkinSystem<dim, spacedim>>(
      GalerkinSystem<dim, spacedim>(this->dof_handler.get()));
  this->cg_sys->update_system_size();
}
template <unsigned int dim, unsigned int spacedim, unsigned int patch_order>
void RegularSolver<dim, spacedim, patch_order>::fill_system() {
  // First, we must make sure that the DoFHandler is up-to-date
  dof_handler->distribute_dofs(&fe_collection_collector);
  std::cout << "N active DoFs: " << dof_handler->get_n_active_dofs()
            << std::endl;
  cg_sys->update_system_size();

  // Now, since we do not want to perform unnecessary integrals, we set the DoF
  // Mask
  cg_sys->set_mask_from_dof_activity();

  // Now we need to fill the system
  cg_sys->fill_system(integrator, *mat_dom, *excitation, ngl_regular,
                      ngl_vertex, ngl_edge, ngl_self);
  cg_sys->fill_forward_excitation(integrator, *excitation, *mat_dom, ngl_regular);
}
template <unsigned int dim, unsigned int spacedim, unsigned int patch_order>
void RegularSolver<dim, spacedim, patch_order>::assemble_system_matrix() {
  this->matrix.reset();
  this->matrix = std::make_unique<MatrixSolving::Matrix<
      std::complex<double>, dromon::MatrixSolvingPolicies::MKLPolicy>>(
      MatrixSolving::Matrix<std::complex<double>,
                            dromon::MatrixSolvingPolicies::MKLPolicy>(
          this->dof_handler->get_n_active_dofs(),
          this->dof_handler->get_n_active_dofs()));

  matrix->matrix_from_galerkin_system(this->cg_sys.get(),
                                      this->dof_handler.get());
  matrix->RHS_from_galerkin_system(this->cg_sys.get(), this->dof_handler.get());
}
template <unsigned int dim, unsigned int spacedim, unsigned int patch_order>
void RegularSolver<dim, spacedim, patch_order>::solve_system() {
  matrix->solve_forward(false);
  this->current_forward_solution = matrix->get_solution_copy();
}
template <unsigned int dim, unsigned int spacedim, unsigned int patch_order>
void RegularSolver<dim, spacedim, patch_order>::conduct_solution_step() {
  // Runs through all the other options...
  this->fill_system();
  this->assemble_system_matrix();
  this->solve_system();
}

template <unsigned int dim, unsigned int spacedim, unsigned int patch_order>
Point<spacedim, std::complex<double>>
RegularSolver<dim, spacedim, patch_order>::post_process_scattered_field() {
  // For monostatic scattering
  return PostProcessing::compute_E_scattered(
      this->dof_handler.get(), current_forward_solution, *mat_dom, *excitation,
      theta_inc, phi_inc, ngl_regular);
}

template <unsigned int dim, unsigned int spacedim, unsigned int patch_order>
void RegularSolver<dim, spacedim, patch_order>::set_plane_wave_excitation(
    const double &frequency, const Point<dim, std::complex<double>> &E_mag) {
  this->excitation = std::make_unique<Excitations::PlaneWave<spacedim, double>>(
      Excitations::PlaneWave<spacedim, double>(frequency, E_mag, theta_inc,
                                               phi_inc));
  this->freq = frequency;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Adaptive Refinement Program
template <unsigned int dim, unsigned int spacedim, unsigned int patch_order>
class AdaptiveSolver {
public:
  AdaptiveSolver(Mesh<dim, spacedim, patch_order> *mesh, Mesh<dim, spacedim, patch_order> *mesh1,
                 MaterialData<double> *mat_dom, MaterialData<double> *mat_dom1,
                 const unsigned int &min_expansion_order = 1,
                 const unsigned int &max_expansion_order = 1,
                 const unsigned int &ngl_regular = 5,
                 const unsigned int &ngl_vertex = 5,
                 const unsigned int &ngl_edge = 5,
                 const unsigned int &ngl_self = 5)
      : min_expansion_order(min_expansion_order),
        max_expansion_order(max_expansion_order), ngl_regular(ngl_regular),
        ngl_vertex(ngl_vertex), ngl_edge(ngl_edge), ngl_self(ngl_self) {
    this->mesh = mesh;
    this->mat_dom = mat_dom;
    this->mesh1 = mesh1;
    this->mat_dom1 = mat_dom1;

    this->initialize_fe_collections();
    this->initialize_dof_handler();

    this->initialize_fe_collections_HOPS();
    this->initialize_dof_handler_HOPS();

  }
  auto execute_refinement(const unsigned int& adjoint_starting_index, const double& reltol, std::vector<std::complex<double>> &forward_matrix, std::vector<std::complex<double>> &forward_excitation,
   std::vector<std::complex<double>> &forward_solution, std::vector<std::complex<double>> &adjoint_solution,
   std::complex<double> &gradient_solution, std::complex<double> &gradient2_solution, double sidelength, int HOPSflag, int perturbVar, double perturbSize);
  void conduct_solution_step();
  double estimate_error();
  bool refine_mesh();
  void solve(const unsigned int& adjoint_starting_index);

  Point<spacedim, std::complex<double>> post_process_scattered_field();
  Point<spacedim, std::complex<double>> post_process_scattered_field_HO();
  double post_process_RCS();


  void set_plane_wave_excitation(const double &frequency,
                                 const Point<dim, std::complex<double>> &E_mag);
  void
  set_scattering_parameters(const double &theta_sc, const double &phi_sc,
                            const double &R_dist_scalar,
                            const Point<spacedim, double> &isolation_direction =
                                Point<spacedim, double>());

private:
  std::vector<std::complex<double>> error_contributions;
  std::complex<double> current_total_error;
  double abs_min_error;
  double abs_max_error;

  void output_data(const std::string file_name);
  void compute_QoI_using_adjoint_solution();
  void HOPS();
  void get_nHat();
  Point<3, double> get_faceCenter(int idx);
  void initialize_DoF_systems(const unsigned int &initial_fe_index);
  void fill_forward_system();
  void fill_HOPS_system();
  void fill_adjoint_excitation();
  void assemble_and_solve_system_matrix();
  void hops_system();

  void initialize_fe_collections() {
    for (unsigned int exps = min_expansion_order; exps <= max_expansion_order;
         ++exps) {
      auto current_fe =
          FE_HdivMaxOrtho<dim, spacedim, double>(exps, CurrentType::Electric);
      fe_collection_EFIE.push_back(&current_fe);
    }
    fe_collection_collector.push_back(fe_collection_EFIE);
  }

  void initialize_fe_collections_HOPS() {
    std::cout << "max expansion order: " << max_expansion_order << std::endl;
    for (unsigned int exps = min_expansion_order; exps <= max_expansion_order;
         ++exps) {
      auto current_fe =
          FE_HdivMaxOrtho<dim, spacedim, double>(exps, CurrentType::Electric);
      fe_collection_EFIE_HOPS.push_back(&current_fe);
    }
    fe_collection_collector_HOPS.push_back(fe_collection_EFIE_HOPS);

    for (unsigned int exps = min_expansion_order; exps <= max_expansion_order;
         ++exps) {
      auto current_fe =
          FE_HdivMaxOrtho<dim, spacedim, double>(exps, CurrentType::Electric);
      fe_collection_EFIE_HOPS2.push_back(&current_fe);
    }
    fe_collection_collector_HOPS2.push_back(fe_collection_EFIE_HOPS2);
  }

  void initialize_dof_handler() {
    this->dof_handler = std::make_unique<DoFHandler<dim, spacedim>>(
        DoFHandler<dim, spacedim>(*mesh));
    this->dof_handler_projection = std::make_unique<DoFHandler<dim, spacedim>>(
        DoFHandler<dim, spacedim>(*mesh));
  }

  void initialize_dof_handler_HOPS() {
    this->dof_handler_HOPS = std::make_unique<DoFHandler<dim, spacedim>>(
        DoFHandler<dim, spacedim>(*this->mesh1));
    this->dof_handler_projection_HOPS = std::make_unique<DoFHandler<dim, spacedim>>(
        DoFHandler<dim, spacedim>(*this->mesh1));

    this->dof_handler_HOPS2 = std::make_unique<DoFHandler<dim, spacedim>>(
        DoFHandler<dim, spacedim>(*this->mesh1));
    this->dof_handler_projection_HOPS2 = std::make_unique<DoFHandler<dim, spacedim>>(
        DoFHandler<dim, spacedim>(*this->mesh1));
  }
  bool check_error_estimation_correctness = true;

  const unsigned int min_expansion_order;
  const unsigned int max_expansion_order;

  const unsigned int ngl_regular;
  const unsigned int ngl_self;
  const unsigned int ngl_vertex;
  const unsigned int ngl_edge;

  Mesh<dim, spacedim, patch_order> *mesh;
  Mesh<dim, spacedim, patch_order> *mesh1;
  MaterialData<double> *mat_dom;
  MaterialData<double> *mat_dom1;

  FECollectionCollector<dim, spacedim, double> fe_collection_collector;
  FECollection<dim, spacedim, double> fe_collection_EFIE;

  FECollectionCollector<dim, spacedim, double> fe_collection_collector_HOPS;
  FECollection<dim, spacedim, double> fe_collection_EFIE_HOPS;

  FECollectionCollector<dim, spacedim, double> fe_collection_collector_HOPS2;
  FECollection<dim, spacedim, double> fe_collection_EFIE_HOPS2;

  std::unique_ptr<DoFHandler<dim, spacedim>> dof_handler;
  std::unique_ptr<DoFHandler<dim, spacedim>> dof_handler_projection;

  std::unique_ptr<DoFHandler<dim, spacedim>> dof_handler_HOPS;
  std::unique_ptr<DoFHandler<dim, spacedim>> dof_handler_projection_HOPS;

    std::unique_ptr<DoFHandler<dim, spacedim>> dof_handler_HOPS2;
  std::unique_ptr<DoFHandler<dim, spacedim>> dof_handler_projection_HOPS2;

  std::unique_ptr<GalerkinSystem<dim, spacedim>> cg_sys;
  std::unique_ptr<GalerkinSystem<dim, spacedim>> cg_sys2;

  std::unique_ptr<GalerkinSystem<dim, spacedim>> cg_sys_HOPS;
  std::unique_ptr<GalerkinSystem<dim, spacedim>> cg_sys_HOPS2;


  std::complex<double> out2;
  double sidelength;


  double freq;
  const double theta_inc = constants<double>::PI / 2.0;
  const double phi_inc = 0.0;

  std::unique_ptr<Excitations::PlaneWave<spacedim, double>> excitation;
  std::unique_ptr<Excitations::PlaneWave<spacedim, double>> excitation_HOPS;
  std::unique_ptr<Excitations::PlaneWave<spacedim, double>> excitation_HOPS2;

  // Here we define pointers for the implemented adjoint excitations,
  // though only one may be active at a time currently
  // In the event of multiple QoIs, these could be vectors of unique pointers
  // for different QoIs
  std::unique_ptr<Excitations::AdjointScatteredFieldExcitation<
      dromon::DoFParent<dim, spacedim, double>, std::complex<double>, double,
      double>>
      adjoint_excitation_Esc;

  std::unique_ptr<Excitations::AdjointScatteredFieldExcitation<
      dromon::DoFParent<dim, spacedim, double>, std::complex<double>, double,
      double>>
      adjoint_excitation_Esc2;

  std::unique_ptr<Excitations::AdjointRCSExcitation<
      dromon::DoFParent<dim, spacedim, double>, std::complex<double>, double,
      double>>
      adjoint_excitation_RCS;

  EFIEIntegrator<DoFParent<dim, spacedim>, std::complex<double>> integrator;
  EFIEIntegrator<DoFParent<dim, spacedim>, std::complex<double>> integrator_HOPS;
  EFIEIntegrator<DoFParent<dim, spacedim>, std::complex<double>> integrator_HOPS2;


  std::vector<std::complex<double>> current_forward_solution;
  std::vector<std::complex<double>> current_adjoint_solution;
  std::vector<std::complex<double>> current_adjoint_rhs;
  std::vector<std::complex<double>> current_adjoint_rhs2;
  std::vector<std::complex<double>> forward_excitation_HO;
  std::vector<std::complex<double>> forward_solution_ho;

  std::vector<std::complex<double>> current_forward_solution_HOPS;
  std::vector<std::complex<double>> current_forward_solution_HOPS2;

  std::vector<std::complex<double>> forward_excitation_HO_HOPS;
  std::vector<std::complex<double>> forward_excitation_HO_HOPS2;

  std::vector<std::complex<double>> forward_solution_ho_HOPS;
  std::vector<std::complex<double>> forward_solution_ho_HOPS2;

  std::vector<std::complex<double>> forward_system_HO_HOPS;
  std::vector<std::complex<double>> forward_system_HO_HOPS2;

  std::vector<std::complex<double>> solution_gradient;
  std::vector<std::complex<double>> matrix_gradient;

  std::vector<std::complex<double>> solution_gradient2;
  std::vector<std::complex<double>> matrix_gradient2;
  std::vector<std::complex<double>> excitation_gradient2;

  std::vector<std::complex<double>> forward_matrix;
  std::vector<std::complex<double>> HOPS_matrix;
  std::vector<std::complex<double>> inverse_matrix;
  int HOPSflag;
  int perturbVar;
  double perturbSize;


  /*
   * Parameters for the scattering QoI choices...
   */
  double reltol;
  double abstol;
  unsigned int max_refinement_iters = 12;

  bool use_Esc_QoI = false;
  bool use_RCS_QoI = false;
  double theta_sc;
  double phi_sc;
  double R_dist_scalar;
  Point<spacedim, double> isolation_direction;
  Point<spacedim, std::complex<double>> Esc_value;
  std::complex<double> gradient;
  std::complex<double> gradient2;
  std::vector<double> nHat;
  double xDir;
  double yDir;
  double zDir;
  std::complex<double> qoi_out;
  std::complex<double> qoi_out_ho;
  std::complex<double> rcs_val;
};

//maybe need a dof_handler for HOPS
template <unsigned int dim, unsigned int spacedim, unsigned int patch_order>
void AdaptiveSolver<dim, spacedim, patch_order>::initialize_DoF_systems(
    const unsigned int &initial_fe_index) {
  this->dof_handler->distribute_dofs(&fe_collection_collector,
                                     initial_fe_index);
                                     
  this->dof_handler_HOPS->distribute_dofs(&fe_collection_collector_HOPS,
                                     initial_fe_index);

  this->dof_handler_HOPS2->distribute_dofs(&fe_collection_collector_HOPS2,
                                     initial_fe_index);                                     

  this->cg_sys = std::make_unique<GalerkinSystem<dim, spacedim>>(
      GalerkinSystem<dim, spacedim>(this->dof_handler.get()));
  this->cg_sys->update_system_size();

  this->cg_sys2 = std::make_unique<GalerkinSystem<dim, spacedim>>(
      GalerkinSystem<dim, spacedim>(this->dof_handler.get()));
  this->cg_sys2->update_system_size();

  this->cg_sys_HOPS = std::make_unique<GalerkinSystem<dim, spacedim>>(
      GalerkinSystem<dim, spacedim>(this->dof_handler_HOPS.get()));
  this->cg_sys_HOPS->update_system_size();

  this->cg_sys_HOPS2 = std::make_unique<GalerkinSystem<dim, spacedim>>(
      GalerkinSystem<dim, spacedim>(this->dof_handler_HOPS2.get()));
  this->cg_sys_HOPS2->update_system_size();  
}

//needs to be seperated out into fill_forward_system and fill_HOPS_system
//since excitations can only be done one at a time
template <unsigned int dim, unsigned int spacedim, unsigned int patch_order>
void AdaptiveSolver<dim, spacedim, patch_order>::fill_forward_system() {

  // First, we must make sure that the DoFHandler is up-to-date
  dof_handler->distribute_dofs(&fe_collection_collector);
  std::cout << "N active DoFs: " << dof_handler->get_n_active_dofs()
            << std::endl;
  cg_sys->update_system_size();
  cg_sys2->update_system_size();


  // Now, since we do not want to perform unnecessary integrals, we set the DoF
  // Mask
  cg_sys->set_mask_from_dof_activity();
  cg_sys2->set_mask_from_dof_activity();



  // Now we need to fill the system
  cg_sys->fill_system(integrator, *mat_dom, *excitation, ngl_regular,
                      ngl_vertex, ngl_edge, ngl_self);
  cg_sys->fill_forward_excitation(integrator, *excitation, *mat_dom, ngl_regular);


 

  // cg_sys->template fill_adjoint_excitation();
}

template <unsigned int dim, unsigned int spacedim, unsigned int patch_order>
void AdaptiveSolver<dim, spacedim, patch_order>::fill_HOPS_system() {


  dof_handler_HOPS->distribute_dofs(&fe_collection_collector_HOPS);
  cg_sys_HOPS->update_system_size();
  cg_sys_HOPS->set_mask_from_dof_activity();

  dof_handler_HOPS2->distribute_dofs(&fe_collection_collector_HOPS2);
  cg_sys_HOPS2->update_system_size();
  cg_sys_HOPS2->set_mask_from_dof_activity();

  // 1 for perturbation method, 2 for analytical
  int HOPSmode = 1;
  switch (HOPSmode)
  {
    case 1:
    {

        // // // This version is for mesh perturbations // // //
        // // // system is filled as normal             // // //
        // cg_sys_HOPS->fill_system(integrator_HOPS, *mat_dom, *excitation_HOPS, ngl_regular,
        //                   ngl_vertex, ngl_edge, ngl_self);
        // cg_sys_HOPS->fill_forward_excitation(integrator_HOPS, *excitation_HOPS, *mat_dom, ngl_regular);

        // // //  This version is for frequency/eps perturbations (theoretically r as well) // // // 
        // // //  perturbations here are multiplied directly within the matrix filling      // // // 
        cg_sys_HOPS->fill_system_perturb(integrator_HOPS, *mat_dom, *excitation_HOPS , this->perturbVar, this->perturbSize, ngl_regular, ngl_vertex, ngl_edge, ngl_self);
        cg_sys_HOPS->fill_forward_excitation_perturb(integrator_HOPS, *excitation_HOPS, *mat_dom, this->perturbVar, this->perturbSize, ngl_regular);

        // // // This version is for calculating second order terms                        // // //
        // // // HOPS would be for x(i+h), HOPS2 is x(i-h), and the forward system is x(i) // // //
        // // // then a second order central FD can be performed to calculated second order terms // //
        cg_sys_HOPS2->fill_system_perturb(integrator_HOPS2, *mat_dom, *excitation_HOPS2 , this->perturbVar, (2.0-this->perturbSize), ngl_regular, ngl_vertex, ngl_edge, ngl_self);
        cg_sys_HOPS2->fill_forward_excitation_perturb(integrator_HOPS2, *excitation_HOPS2, *mat_dom, this->perturbVar, (2.0-this->perturbSize), ngl_regular);        
      
      // cg_sys_HOPS->fill_forward_excitation(integrator_HOPS, *excitation_HOPS, ngl_regular);
      // cg_sys_HOPS->fill_system(integrator_HOPS, *mat_dom, *excitation_HOPS, ngl_regular, ngl_vertex, ngl_edge, ngl_self);

      break;
    }
    case 2:
    {

      cg_sys_HOPS->fill_system_HOPS(integrator_HOPS, *mat_dom, *excitation_HOPS, ngl_regular,
                          ngl_vertex, ngl_edge, ngl_self);
      cg_sys_HOPS->fill_forward_excitation_HOPS(integrator_HOPS, *excitation_HOPS, *mat_dom, ngl_regular);
      break;

    }

  }
}

template <unsigned int dim, unsigned int spacedim, unsigned int patch_order>
void AdaptiveSolver<dim, spacedim,
                    patch_order>::hops_system() {

      //build the hops system matrix
        dromon::MatrixSolving::Matrix<std::complex<double>,
                                    dromon::MatrixSolvingPolicies::MKLPolicy>
        matrix_hi_HOPS(dof_handler->get_n_active_dofs(),
                  dof_handler->get_n_active_dofs());

        dromon::MatrixSolving::Matrix<std::complex<double>,
                                    dromon::MatrixSolvingPolicies::MKLPolicy>
        matrix_hi_HOPS2(dof_handler->get_n_active_dofs(),
                  dof_handler->get_n_active_dofs());

    // this->excitation_HOPS.reset();
    this->excitation_HOPS = std::make_unique<Excitations::PlaneWave<spacedim, double>>(
      Excitations::PlaneWave<spacedim, double>(freq, {1.0, 0.0}, theta_inc,
                                               phi_inc, *mat_dom));

    this->excitation_HOPS2 = std::make_unique<Excitations::PlaneWave<spacedim, double>>(
      Excitations::PlaneWave<spacedim, double>(freq, {1.0, 0.0}, theta_inc,
                                               phi_inc, *mat_dom));                                               

    this->fill_HOPS_system();

    matrix_hi_HOPS.matrix_from_galerkin_system(cg_sys_HOPS.get(), dof_handler_HOPS.get());
    matrix_hi_HOPS.RHS_from_galerkin_system(cg_sys_HOPS.get(), dof_handler_HOPS.get());
    this->forward_excitation_HO_HOPS  = matrix_hi_HOPS.get_RHS_copy();
    std::vector<std::complex<double>> perturb_mat = matrix_hi_HOPS.get_matrix_copy();
    matrix_hi_HOPS.solve_forward(false);
    this->forward_solution_ho_HOPS  = matrix_hi_HOPS.get_solution_copy();

    matrix_hi_HOPS2.matrix_from_galerkin_system(cg_sys_HOPS2.get(), dof_handler_HOPS2.get());
    matrix_hi_HOPS2.RHS_from_galerkin_system(cg_sys_HOPS2.get(), dof_handler_HOPS2.get());
    this->forward_excitation_HO_HOPS2  = matrix_hi_HOPS2.get_RHS_copy();
    std::vector<std::complex<double>> perturb_mat2 = matrix_hi_HOPS2.get_matrix_copy();
    matrix_hi_HOPS.solve_forward(false);
    this->forward_solution_ho_HOPS2  = matrix_hi_HOPS2.get_solution_copy();    


  // 1 for perturbation method, 2 for analytical
  int HOPSmode = 1;
  switch (HOPSmode)
  {
    case 1:
    {
      int m = this->forward_excitation_HO_HOPS.size();
      std::vector<std::complex<double>> delta_mat(m*m);
      std::vector<std::complex<double>> delta_vec(m);
      std::vector<std::complex<double>> delta_solution(m);

      std::vector<std::complex<double>> delta_mat2(m*m);
      std::vector<std::complex<double>> delta_vec2(m);
      std::vector<std::complex<double>> delta_solution2(m);
      std::vector<std::complex<double>> delta_excitation2(m);

      double division1, division2;
      std::complex<double> normL;
      //take the difference between perturbed and forward matrix
      for (int i = 0; i < m*m; ++i){

        // //perturb material (this isn't quite right, should be eps*0.01 and not eps_0 *0.01)
        if (this->perturbVar == 1){
          double eps = mat_dom->get_material_domain(0).get_exterior().epsr;
          // delta_mat[i] = (perturb_mat[i] - this->forward_matrix[i]) / (8.8541878128e-12 * 0.01);
          delta_mat[i] = (perturb_mat[i] - this->forward_matrix[i]) / ((this->perturbSize - 1.0) * eps);
          delta_mat2[i] = (perturb_mat[i] - 2.0*this->forward_matrix[i] + perturb_mat2[i]) / (eps * (this->perturbSize - 1.0)) / (eps * (this->perturbSize - 1.0));

          // std::cout << "division value: " << ((this->perturbSize - 1.0) * eps) << std::endl;
        }
        //perturb frequency
        else if (this->perturbVar == 2){
          delta_mat[i] = (perturb_mat[i] - this->forward_matrix[i]) / (this->freq * 2 * 3.1415926535 * (this->perturbSize - 1.0));
          // std::cout << "division value f: " << (this->freq * 2 * 3.1415926535 * (this->perturbSize - 1.0)) << std::endl;

        }
        //perturb radius
        else if (this->perturbVar == 3){
          delta_mat[i] = (perturb_mat[i] - this->forward_matrix[i]) / ((this->perturbSize - 1.0)*this->sidelength);
          // std::cout << "division value r: " << ((this->perturbSize - 1.0) * this->sidelength) << std::endl;
          division1 = ((this->perturbSize - 1.0)*this->sidelength);

        }
      }

      std::complex<double> normLu;
      std::complex<double> normU;
      std::vector<std::complex<double>> product(m);
      //multiple by forward solution
      for (int i = 0; i < m; ++i){
        for (int j = 0; j < m; ++j){
          product[i] += delta_mat[i*m+j]*this->forward_solution_ho[j];
        }
      }

      for (int i = 0; i < m; ++i){

        // perturb material
        if (this->perturbVar == 1){
          double eps = mat_dom->get_material_domain(0).get_exterior().epsr;
          // delta_vec[i] = (this->forward_excitation_HO_HOPS[i] - this->forward_excitation_HO[i]) / (8.8541878128e-12 * (this->perturbSize - 1.0));
          delta_vec[i] = (this->forward_excitation_HO_HOPS[i] - this->forward_excitation_HO[i]) / (eps * (this->perturbSize - 1.0));

          delta_solution[i] = (this->forward_solution_ho_HOPS[i] - this->forward_solution_ho[i]) / (eps * (this->perturbSize - 1.0));
          delta_solution2[i] = (this->forward_solution_ho_HOPS[i] - 2.0 * this->forward_solution_ho[i] + this->forward_solution_ho_HOPS2[i]) / (eps * (this->perturbSize - 1.0)) / (eps * (this->perturbSize - 1.0));
          delta_excitation2[i] = (this->forward_excitation_HO_HOPS[i] - 2.0 * this->forward_excitation_HO[i] + this->forward_excitation_HO_HOPS2[i]) / (eps * (this->perturbSize - 1.0)) / (eps * (this->perturbSize - 1.0));

          // delta_solution[i] = 0.0;
        }
        //pertrub freq
        else if (this->perturbVar == 2){
          delta_vec[i] = (this->forward_excitation_HO_HOPS[i] - this->forward_excitation_HO[i]) / (this->freq * 2 * 3.1415926535 * (this->perturbSize - 1.0));

          delta_solution[i] = (this->forward_solution_ho_HOPS[i] - this->forward_solution_ho[i]) / (this->freq * 2 * 3.1415926535 * (this->perturbSize - 1.0));
        }
        // perturb radius
        else if (this->perturbVar == 3){
          delta_vec[i] = (this->forward_excitation_HO_HOPS[i] - this->forward_excitation_HO[i]) / ((this->perturbSize - 1.0)*this->sidelength);
          division2 = ((this->perturbSize - 1.0)*this->sidelength);

          delta_solution[i] = (this->forward_solution_ho_HOPS[i] - this->forward_solution_ho[i]) / ((this->perturbSize - 1.0)*this->sidelength);
        }
      }
      //store
      // auto forward_system_HO_HOPS = product;
      this->forward_system_HO_HOPS = product;
      this->forward_excitation_HO_HOPS = delta_vec;
      this->solution_gradient = delta_solution;
      this->solution_gradient2 = delta_solution2;
      this->matrix_gradient = delta_mat;
      this->matrix_gradient2 = delta_mat2;
      this->excitation_gradient2 = delta_excitation2;

      std::cout << "division1: " << division1 << "   dvision2: " << division2 << std::endl;

      break;
    }

    case 2:
    {
      auto forward_system_HO_HOPS = matrix_hi_HOPS.multiply(this->forward_solution_ho);
      this->forward_system_HO_HOPS = forward_system_HO_HOPS;
      this->HOPS_matrix = matrix_hi_HOPS.get_matrix_copy();
      this->matrix_gradient = this->HOPS_matrix;

      break;
    }
    
  }

  // for (int i = 0; i < this->forward_system_HO_HOPS.size(); ++i){
  // std::cout << "deltaL*u: " << forward_system_HO_HOPS[i] << "    deltaG[i]: " << forward_excitation_HO_HOPS[i] << std::endl;
  // }

}

template <unsigned int dim, unsigned int spacedim, unsigned int patch_order>
void AdaptiveSolver<dim, spacedim, patch_order>::fill_adjoint_excitation() {
  if (use_Esc_QoI){
    cg_sys->fill_adjoint_excitation(*adjoint_excitation_Esc, ngl_regular);
    cg_sys->fill_adjoint_excitation_perturb(*adjoint_excitation_Esc2, ngl_regular, this->perturbSize);
  }
  else if (use_RCS_QoI)
    cg_sys->fill_adjoint_excitation(*adjoint_excitation_RCS, ngl_regular);
  else {
    #ifdef DEBUG
      assert(false && "Case not implemented!");
    #endif
    return;
  }
}

template <unsigned int dim, unsigned int spacedim, unsigned int patch_order>
void AdaptiveSolver<dim, spacedim,
                    patch_order>::assemble_and_solve_system_matrix() {

  // Solve the lower order (forward) problem
  {
    this->dof_handler->spawn_p_reduced_dof_handler(
        dof_handler_projection.get());
    this->dof_handler_projection->distribute_dofs();
    dromon::MatrixSolving::Matrix<std::complex<double>,
        dromon::MatrixSolvingPolicies::MKLPolicy>
        matrix_lo(dof_handler_projection->get_n_active_dofs(),
        dof_handler_projection->get_n_active_dofs());

    matrix_lo.matrix_from_galerkin_system(cg_sys.get(),
                                          dof_handler_projection.get());
    matrix_lo.RHS_from_galerkin_system(cg_sys.get(),
                                       dof_handler_projection.get());

    matrix_lo.solve_forward(false);
    current_forward_solution = matrix_lo.get_solution_copy();

  }

  this->Esc_value = this->post_process_scattered_field();

  //jake:
  // Now that we have the lower order forward solution, we must compute the
  // adjoint solution First, reset the unique_ptr for the adjoint excitation

  //Me (maybe):
  //order here doesnt matter, the cg_sys has a different adjoint excitation than forward excitation.  
  //So the excitation pointer must be reset before filling, but this could be done directly after filling the forward excitation
  if (use_Esc_QoI) {
    this->qoi_out = isolation_direction.dot(Esc_value);
    this->adjoint_excitation_Esc.reset();
    this->adjoint_excitation_Esc =
          std::make_unique<Excitations::AdjointScatteredFieldExcitation<
          dromon::DoFParent<dim, spacedim, double>,
          std::complex<double>, double, double>>(Excitations::AdjointScatteredFieldExcitation<
          dromon::DoFParent<dim, spacedim, double>,
          std::complex<double>, double, double>(
          freq, theta_sc, phi_sc, isolation_direction, R_dist_scalar));
    // Now fill the adjoint excitation
    cg_sys->fill_adjoint_excitation(*adjoint_excitation_Esc,
                                    this->ngl_regular);

    this->adjoint_excitation_Esc2.reset();
    this->adjoint_excitation_Esc2 =
          std::make_unique<Excitations::AdjointScatteredFieldExcitation<
          dromon::DoFParent<dim, spacedim, double>,
          std::complex<double>, double, double>>(Excitations::AdjointScatteredFieldExcitation<
          dromon::DoFParent<dim, spacedim, double>,
          std::complex<double>, double, double>(
          freq, theta_sc, phi_sc, isolation_direction, R_dist_scalar));
    // Now fill the adjoint excitation
    cg_sys2->fill_adjoint_excitation_perturb(*adjoint_excitation_Esc2,
                                    this->ngl_regular, this->perturbSize);                                    

  } else if (use_RCS_QoI) {
    this->qoi_out = this->post_process_RCS();
    // We first have to use the forward solution to compute the scattered
    // field
    this->adjoint_excitation_RCS.reset();
    this->adjoint_excitation_RCS = std::make_unique<Excitations::AdjointRCSExcitation<
        dromon::DoFParent<dim, spacedim, double>, std::complex<double>,
        double, double>>(
        Excitations::AdjointRCSExcitation<
        dromon::DoFParent<dim, spacedim, double>, std::complex<double>,
        double, double>(freq, theta_sc, phi_sc, Esc_value, R_dist_scalar));
    // Now fill the adjoint excitation
    cg_sys->fill_adjoint_excitation(*adjoint_excitation_RCS,
                                    this->ngl_regular);
  }

  // Solve the higher order (adjoint) problem
  {
    dromon::MatrixSolving::Matrix<std::complex<double>,
        dromon::MatrixSolvingPolicies::MKLPolicy>
        matrix_hi(dof_handler->get_n_active_dofs(),
        dof_handler->get_n_active_dofs());

    dromon::MatrixSolving::Matrix<std::complex<double>,
        dromon::MatrixSolvingPolicies::MKLPolicy>
        matrix_hi2(dof_handler->get_n_active_dofs(),
        dof_handler->get_n_active_dofs());

    matrix_hi.matrix_from_galerkin_system(cg_sys.get(), dof_handler.get());
    matrix_hi.adjoint_RHS_from_galerkin_system(cg_sys.get(), dof_handler.get());
    this->current_adjoint_rhs  = matrix_hi.get_RHS_copy();

    // matrix_hi2.matrix_from_galerkin_system(cg_sys2.get(), dof_handler.get());
    matrix_hi2.adjoint_RHS_from_galerkin_system(cg_sys2.get(), dof_handler.get());
    this->current_adjoint_rhs2 = matrix_hi2.get_RHS_copy();

    auto t1 = std::chrono::high_resolution_clock::now();
    matrix_hi.solve_adjoint(false);
    auto t2 = std::chrono::high_resolution_clock::now();
    auto invert_time = std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1);
    // auto invert_time = t2-t1;

    std::cout << "time for matrix inverse:  " << invert_time.count() << std::endl;
    this->current_adjoint_solution  = matrix_hi.get_solution_copy();
  }

  if (this->check_error_estimation_correctness)
  {
    //this is the higher order forward problem
    dromon::MatrixSolving::Matrix<std::complex<double>,
        dromon::MatrixSolvingPolicies::MKLPolicy>
        matrix_hi(dof_handler->get_n_active_dofs(),
        dof_handler->get_n_active_dofs());

    //this entry is temporary.  Once the forward problem is solved, the LHS and RHS are changed
    //this is a copy of the forward system that will not be solve so that the forward matrix is unaltered
    //this will be used for checking the L J_S gradient
    dromon::MatrixSolving::Matrix<std::complex<double>,
        dromon::MatrixSolvingPolicies::MKLPolicy>
        matrix_hi_unsolved(dof_handler->get_n_active_dofs(),
        dof_handler->get_n_active_dofs());

    matrix_hi.matrix_from_galerkin_system(cg_sys.get(), dof_handler.get());
    matrix_hi.RHS_from_galerkin_system(cg_sys.get(), dof_handler.get());
    this->forward_excitation_HO  = matrix_hi.get_RHS_copy();
    this->forward_matrix = matrix_hi.get_matrix_copy();
    matrix_hi.solve_forward(false);
    this->forward_solution_ho  = matrix_hi.get_solution_copy();

    matrix_hi_unsolved.matrix_from_galerkin_system(cg_sys.get(), dof_handler.get());
    matrix_hi_unsolved.RHS_from_galerkin_system(cg_sys.get(), dof_handler.get());

    int m = forward_excitation_HO.size();
    std::vector<std::complex<double>> test_vec1 = matrix_hi_unsolved.multiply(this->forward_solution_ho);
    std::vector<std::complex<double>> test_vec(m);

    std::complex<double> matVal;
    std::complex<double> transpVal;
    std::vector<std::complex<double>> transpMat(m*m);
    std::vector<std::complex<double>> multMat(m*m);
    std::complex<double> matSum;


    // // // LAPACKE_zgetrs(is_row_major ? LAPACK_ROW_MAJOR : LAPACK_COL_MAJOR, solve_adjoint ? 'C' : 'N',
    // // //                                       int(m_rows), 1, matrix_data, lda, &ipiv[0], rhs_data, ldb)

    // std::cout << "m size: " << m << std::endl;
    int* ipiv = new int[m];
    lapack_int lda = m;
    std::vector<std::complex<double>> inverse_mat = forward_matrix;
    // std::vector<std::complex<double>> inverse_mat_conjTransp = forward_matrix;
    // // std::vector<std::complex<double>> test_vec1(m);
    
    // // std::complex<double>* mat = &forward_matrix[0];
    auto error = LAPACKE_zgetrf(LAPACK_ROW_MAJOR, int(m), int(m), &inverse_mat[0], lda, &ipiv[0]);
    auto ierror = LAPACKE_zgetri(LAPACK_ROW_MAJOR, m, &inverse_mat[0], m, &ipiv[0]);

    this->inverse_matrix = inverse_mat;


    // std::vector<std::complex<double>> inverse2_mat = inverse_mat;
    // auto error2 = LAPACKE_zgetrf(LAPACK_ROW_MAJOR, int(m), int(m), &inverse2_mat[0], lda, &ipiv[0]);
    // auto ierror2 = LAPACKE_zgetri(LAPACK_ROW_MAJOR, m, &inverse2_mat[0], m, &ipiv[0]);

  }

  this->Esc_value = this->post_process_scattered_field_HO();
  qoi_out_ho = isolation_direction.dot(Esc_value);


}

template <unsigned int dim, unsigned int spacedim, unsigned int patch_order>
void AdaptiveSolver<dim, spacedim, patch_order>::conduct_solution_step() {
  // Runs through all the other options...
  this->fill_forward_system();
  //this->fill_HOPS_system();
  this->assemble_and_solve_system_matrix();
  if (HOPSflag)
  this->hops_system();
}

template <unsigned int dim, unsigned int spacedim, unsigned int patch_order>
Point<spacedim, std::complex<double>>
AdaptiveSolver<dim, spacedim, patch_order>::post_process_scattered_field() {
  // For monostatic scattering
  return PostProcessing::compute_E_scattered(
      this->dof_handler_projection.get(), current_forward_solution, *mat_dom, *excitation,
      theta_sc, phi_sc, ngl_regular, this->R_dist_scalar);
}

template <unsigned int dim, unsigned int spacedim, unsigned int patch_order>
Point<spacedim, std::complex<double>>
AdaptiveSolver<dim, spacedim, patch_order>::post_process_scattered_field_HO() {
  // For monostatic scattering
  return PostProcessing::compute_E_scattered(
      this->dof_handler.get(), forward_solution_ho, *mat_dom, *excitation,
      theta_sc, phi_sc, ngl_regular, this->R_dist_scalar);
}

template <unsigned int dim, unsigned int spacedim, unsigned int patch_order>
void AdaptiveSolver<dim, spacedim, patch_order>::set_plane_wave_excitation(
    const double &frequency, const Point<dim, std::complex<double>> &E_mag) {
  this->excitation = std::make_unique<Excitations::PlaneWave<spacedim, double>>(
      Excitations::PlaneWave<spacedim, double>(frequency, E_mag, theta_inc,
                                               phi_inc, *mat_dom));
  this->freq = frequency;
}



template <unsigned int dim, unsigned int spacedim, unsigned int patch_order>
void AdaptiveSolver<dim, spacedim, patch_order>::set_scattering_parameters(
    const double &theta_sc, const double &phi_sc, const double &R_dist_scalar,
    const Point<spacedim, double> &isolation_direction) {

  this->theta_sc = theta_sc;
  this->phi_sc = phi_sc;
  this->R_dist_scalar = R_dist_scalar;
  this->isolation_direction = isolation_direction;
  if (isolation_direction == Point<spacedim, double>())
    use_RCS_QoI = true;
  else 
    use_Esc_QoI = true;

}

template <unsigned int dim, unsigned int spacedim, unsigned int patch_order>
double AdaptiveSolver<dim, spacedim, patch_order>::estimate_error()
{
  dromon::KahanVector<std::complex<double>> error_contributions_temp(mesh->n_cells());
  dromon::ErrorEstimation::EstimateDWR(*dof_handler_projection, *dof_handler,
      current_forward_solution, current_adjoint_solution, *cg_sys,
      error_contributions_temp);
  this->current_total_error = error_contributions_temp.sum();
  this->error_contributions = error_contributions_temp.data();
  std::cout << "Total error: " << current_total_error << std::endl;
  if (use_RCS_QoI) {
    utility::make_real(this->error_contributions);
    this->current_total_error.imag(0);
  }
  // Set max and min errors
  this->abs_max_error = std::abs(error_contributions[0]);
  this->abs_min_error = std::abs(error_contributions[1]);
  for (unsigned int i = 1; i < this->error_contributions.size(); ++i)
  {
    const auto& current_val = std::abs(error_contributions[i]);
    if (abs_max_error < current_val)
      abs_max_error = current_val;
    if (abs_min_error > current_val)
      abs_min_error = current_val;
  }
  // Return the estimated relative error
  return std::abs(this->current_total_error)/(std::abs(this->current_total_error + this->qoi_out));
}

template <unsigned int dim, unsigned int spacedim, unsigned int patch_order>
bool AdaptiveSolver<dim, spacedim, patch_order>::refine_mesh()
{
  this->abstol =
    this->reltol * std::abs(this->current_total_error + qoi_out) / (double(mesh->n_cells())) * 20.0;
  std::cout << "Abstol: " << abstol << std::endl;
  if (this->abs_max_error < abstol)
    return true;
  // Refinement::uniform_p_refinement(dof_handler.get());
  dromon::Refinement::prediction_p_refinement(dof_handler.get(), error_contributions, abstol, -1);
  //dromon::Refinement::prediction_p_refinement(dof_handler_HOPS.get(), error_contributions, abstol, -1);
  return false;
}

template <unsigned int dim, unsigned int spacedim, unsigned int patch_order>
void AdaptiveSolver<dim, spacedim, patch_order>::solve(const unsigned int& adjoint_starting_index)
{
  this->initialize_DoF_systems(adjoint_starting_index);
  this->conduct_solution_step();
  this->compute_QoI_using_adjoint_solution();
}

template <unsigned int dim, unsigned int spacedim, unsigned int patch_order>
auto AdaptiveSolver<dim, spacedim, patch_order>::execute_refinement(const unsigned int& adjoint_starting_index, const double &reltol, std::vector<std::complex<double>> &forward_matrix, 
std::vector<std::complex<double>> &forward_excitation, std::vector<std::complex<double>> &forward_solution, 
std::vector<std::complex<double>> &adjoint_solution, std::complex<double> &gradient_solution, std::complex<double> &gradient2_solution, double sidelength, int HOPSflag, int perturbVar, double perturbSize)
{

  this->perturbVar = perturbVar;
  this->HOPSflag = HOPSflag;
  this->perturbSize = perturbSize;
  // this->get_nHat();
  this->reltol = reltol;
  double estimated_rel_error = 1;
  this->initialize_DoF_systems(adjoint_starting_index);
  bool is_refinement_complete = false;
  this->sidelength = sidelength;
  this->mesh1 = mesh1;

  // for (unsigned int refinement_iter = 0; refinement_iter < max_refinement_iters; ++refinement_iter)
  for (unsigned int refinement_iter = 0; refinement_iter < 1; ++refinement_iter)
  {
    this->conduct_solution_step();
    this->compute_QoI_using_adjoint_solution();
    auto current_relative_error = this->estimate_error();
    std::cout << "Relative error: " << current_relative_error << std::endl;
    this->output_data("AdaptiveSolve-" + std::to_string(refinement_iter));
    is_refinement_complete = this->refine_mesh();
    if (is_refinement_complete)
      break;
  }

  forward_matrix = this->forward_matrix;
  forward_excitation = this->forward_excitation_HO;
  forward_solution = this->forward_solution_ho;
  adjoint_solution = this->current_adjoint_solution;
  if (HOPSflag){
    gradient_solution = this->gradient;
    gradient2_solution = this->gradient2;
  }

  return this->out2;
}

template <unsigned int dim, unsigned int spacedim, unsigned int patch_order>
double AdaptiveSolver<dim, spacedim, patch_order>::post_process_RCS() {
  return PostProcessing::compute_RCS(
      this->dof_handler_projection.get(), current_forward_solution, *mat_dom, *excitation,
      theta_sc, phi_sc, ngl_regular,this->R_dist_scalar);
}

template <unsigned int dim, unsigned int spacedim, unsigned int patch_order>
void AdaptiveSolver<dim, spacedim,
                    patch_order>::compute_QoI_using_adjoint_solution()
{

  std::complex<double> complexj(0, 1.0);
  std::complex<double> out1(0.0);
  double wave = 299792458/this->freq;
  double omega = 2*3.14159*this->freq;
  double m = this->current_adjoint_solution.size();
  
  std::vector<std::complex<double>> test(m);
  for (int i = 0; i < m; ++i){
    test[i] = std::complex<double> {-1.0*current_adjoint_rhs[i].imag(), current_adjoint_rhs[i].real()};
  }

  for (unsigned int i = 0; i < this->current_adjoint_solution.size(); ++i)
  {
    out1 += this->forward_excitation_HO[i]*std::conj(this->current_adjoint_solution[i]);// * wave;

  }

  std::complex<double> out2(0.0);
  for (unsigned int i =0; i < this->current_adjoint_rhs.size(); ++i)
  {
    out2 += this->forward_solution_ho[i]*std::conj(this->current_adjoint_rhs[i]);// * wave;
    // out2 += this->forward_solution_ho[i]*std::conj(test[i]);
  }

  if (HOPSflag){
    std::complex<double> gradient(0.0);
    for (unsigned int i =0; i < this->current_adjoint_solution.size(); ++i)
    {
      //// forward_excitation_HO_HOPS is \delta G, forward_system_HO_HOPS is \deltaL * JS
      gradient += (this->forward_excitation_HO_HOPS[i] - this->forward_system_HO_HOPS[i])*std::conj(this->current_adjoint_solution[i]);// * wave;
      // gradient += -this->forward_system_HO_HOPS[i]*(std::conj(this->current_adjoint_solution[i]));

      // std::cout << "deltaG: " << this->forward_excitation_HO_HOPS[i] << "   deltaL*u: " << this->forward_system_HO_HOPS[i] << std::endl;
    }
    this->gradient = gradient;
    std::cout << "gradient: " << gradient << std::endl;


    //checking gradient alternate formulation
    std::vector<std::complex<double>> vecProd(m);
    for (int i = 0; i < m; ++i){
      for (int j = 0; j < m; ++j){     
        vecProd[i] += this->inverse_matrix[i*m+j] * (this->forward_excitation_HO_HOPS[j] - this->forward_system_HO_HOPS[j]);
      }
    }
    std::complex<double> check;
    for (int i = 0; i < m; ++i){
      check += vecProd[i] * std::conj(this->current_adjoint_rhs[i]);
    }

    std::cout << "gradient check: ----- " << check << std::endl;

    //second order // ------------------------------------------------
    std::vector<std::complex<double>> product(m);
    std::vector<std::complex<double>> gradSol(m);
    std::vector<std::complex<double>> adjoint_excitation_gradient(m);
    std::vector<std::complex<double>> vecProd2(m);
    std::vector<std::complex<double>> vec2(m);
    double eps = mat_dom->get_material_domain(0).get_exterior().epsr;


    for (int i = 0; i < m; ++i){
      adjoint_excitation_gradient[i] = (this->current_adjoint_rhs2[i] - this->current_adjoint_rhs[i]) / (eps * (this->perturbSize - 1.0));
    }

    std::vector<std::complex<double>> mat(m);
    std::vector<std::complex<double>> mat2(m);
    for (int i = 0; i < m; ++i){
      for (int j = 0; j < m; ++j){
        
        mat[i] += matrix_gradient[i*m+j] * vecProd[j];
        // mat2[i] += matrix_gradient2[i*m+j] * this->current_forward_solution[j];
        mat2[i] = 0.0;
      }
    }

    for (int i = 0; i < m; ++i){
      vec2[i] = this->excitation_gradient2[i] -2.0*mat[i] - mat2[i];
    }

    for (int i = 0; i < m; ++i){
      for (int j = 0; j < m; ++j){

        vecProd2[i] += this->inverse_matrix[i*m+j] * vec2[j];

      }
    }


    for (int i = 0; i < m; ++i){
      gradient2 += /*std::conj(adjoint_excitation_gradient[i])*vecProd[i]*vecProd[i]*/ + std::conj(this->current_adjoint_rhs[i]) * vecProd2[i];

    }


    // for (int i = 0; i < m; ++i){
    //   for (int j = 0; j < m; ++j){

    //     gradSol[i] += this->inverse_matrix[i*m+j] * (this->forward_excitation_HO_HOPS[j] - this->forward_system_HO_HOPS[j]);

    //   }
    // }

    // for (unsigned int i =0; i < m; ++i)
    // {
    //   for (unsigned int j =0; j < m; ++j)
    //   {
    //     // product[i] += this->matrix_gradient[i*m+j]*this->solution_gradient[j];
    //     product[i] += this->matrix_gradient[i*m+j]*gradSol[j];
    //   }
    // }

    // for (int i = 0; i < m; ++i){
    //   gradient2 += (-product[i]) * std::conj(this->current_adjoint_solution[i]) * 2.0;
    // }
    this->gradient2 = gradient2;
    std::cout << "second order gradient: " << gradient2 << std::endl << std::endl;
  }


  // for (int i =0; i<m; ++i){
  //   std::cout << this->forward_excitation_HO_HOPS[i] << "   " << this->forward_system_HO_HOPS[i] << std::endl;
  // }


  this->out2 = out2;

  std::complex<double> RCS = 4*3.14159*(R_dist_scalar)*(R_dist_scalar)*out1*std::conj(out1);
  std::complex<double> RCSlam = 4*3.14159*(R_dist_scalar)*(R_dist_scalar)*out1*std::conj(out1) / wave / wave;

  this->rcs_val = RCS;

  std::cout << "Esc: " << out1 << std::endl;
  std::cout << "Esc: " << out2 << std::endl;
  std::cout << "Esc (non-adjoint version): " << qoi_out << std::endl;
  std::cout << "Esc HO (non-adjoint version): " << qoi_out_ho << std::endl;
  std::cout << "Esc/lambda: " << out1/wave << std::endl;
  std::cout << "RCS value: " << RCS << std::endl;
  std::cout << "RCS/(lambda^2) value: " << RCSlam << std::endl;
  
  

  if (this->check_error_estimation_correctness)
  {
    std::cout << "Post-Process Computed Error: " << out2 - this->qoi_out << std::endl;
  }

}

template <unsigned int dim, unsigned int spacedim, unsigned int patch_order>
void AdaptiveSolver<dim, spacedim, patch_order>::output_data(const std::string file_name)
{

  dromon::DataOut<dim, spacedim, patch_order> data_out(
      mesh, file_name, dromon::verbose_output);
  data_out.attach_dof_handler(dof_handler.get());
  std::vector<double> fe_degree_vector, n_dofs_per_cell, abs_error_cell;
  for (const auto &dof_cell : dof_handler->get_dof_cells()) {
    fe_degree_vector.push_back(dof_cell.active_degree(0)-1);
    n_dofs_per_cell.push_back(
        dof_handler->get_n_active_dofs_on_cell(dof_cell.index));
    abs_error_cell.push_back(std::abs(this->error_contributions[dof_cell.index]));

  }
  data_out.add_cell_data(fe_degree_vector, "fe_degrees");
  data_out.add_cell_data(n_dofs_per_cell, "dofs_per_cell");
  data_out.add_cell_data(abs_error_cell, "abs_qoi_error");

  data_out.vtk_out();
}

} // namespace Problems
