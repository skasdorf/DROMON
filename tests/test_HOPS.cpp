#include "DROMON/FE_HdivMaxOrtho.h"
#include "DROMON/Materials.h"
#include "DROMON/MeshGenerator.h"
#include "DROMON/Point.h"
#include "DROMON/config.h"
#include "DROMON/mesh.h"


#include <iostream>
#include "programs.h"
#include <chrono>
#include <fstream>



int main() {
  using namespace dromon;

  //PED
  // mat_dom.push_back(dromon::MaterialDomain<double>(
  //   dromon::Material<double>(5.0, 1.0, false, false),
  //   dromon::Material<double>(1.0, 1.0, false, false)));

  //for gradient finite differences
  std::complex<double> out2;
  std::complex<double> out1;

  std::vector<std::complex<double>> forward_matrix0;
  std::vector<std::complex<double>> forward_matrix1;
  std::vector<std::complex<double>> forward_excitation0;
  std::vector<std::complex<double>> forward_excitation1;
  std::vector<std::complex<double>> forward_solution0;
  std::vector<std::complex<double>> forward_solution1;
  std::vector<std::complex<double>> adjoint_solution0;
  std::vector<std::complex<double>> adjoint_solution1;
  std::complex<double> gradient0;
  std::complex<double> gradient1;
  std::complex<double> gradient2_0;
  std::complex<double> gradient2_1;
  double wave;
  double freq;



  //1 for running a MC, 2 for testing/FD comparison
  int runMode = 1;

  //1 for material params, 2 for frequency, 3 for radius
  int perturbVar = 1;

  //perturbation size
  // double perturbSize = 1.0000001;
  double perturbSize = 1.0001;


  unsigned int n_cells_per_dim = 7;
  freq = 500.0e6;



  switch (runMode){
  case 1:
  {
    //Read in the distribution for MC
    std::vector<double> normalDist;
    switch (perturbVar)
    {
      case 1:
      {
        std::ifstream inFile("./inputDist/normalDist_eps.txt");
        double a;
        while(inFile >> a){
          normalDist.push_back(a);
        }
        break;
      }
      case 2:
      {
        std::ifstream inFile("./inputDist/normalDist_f.txt");
        double a;
        while(inFile >> a){
          normalDist.push_back(a);
        }
        break;
      }
      case 3:
      {
        std::ifstream inFile("./inputDist/normalDist_R.txt");
        double a;
        while(inFile >> a){
          normalDist.push_back(a);
        }
        break;
      }
    }


    for (int i = 0; i < normalDist.size(); ++i){


      std::cout << "New Run _________________________________________________________\n";

      MaterialData<double> mat_dom;
      MaterialData<double> mat_dom1;
      double sidelength;
      std::string saveName;

      switch (perturbVar)
      {
        case 1:
        {
          // Make the MaterialData
          // In this case, a PEC object embedded in air
          //PEC

          mat_dom.push_back(dromon::MaterialDomain<double>(
          dromon::Material<double>(1.0, 1.0, true, false),
          dromon::Material<double>(normalDist[i], 1.0, false, false)));

          mat_dom1.push_back(dromon::MaterialDomain<double>(
          dromon::Material<double>(1.0, 1.0, true, false),
          dromon::Material<double>(normalDist[i]*perturbSize, 1.0, false, false)));
          sidelength = 1.0;
          saveName = "./output/MC_mu80_sig3_plate_noSecondTerms.txt";
          
          break;
        }

        case 2:
        {

          mat_dom.push_back(dromon::MaterialDomain<double>(
          dromon::Material<double>(1.0, 1.0, true, false),
          dromon::Material<double>(1.0, 1.0, false, false)));

          mat_dom1.push_back(dromon::MaterialDomain<double>(
          dromon::Material<double>(1.0, 1.0, true, false),
          dromon::Material<double>(1.0, 1.0, false, false)));

          sidelength = 1.0;
          freq = normalDist[i];
          saveName = "./output/MC_mu400_sig25_Gcheck.txt";

          break;
        }

        case 3:
        {
          // Make the MaterialData
          // In this case, a PEC object embedded in air
          //PEC
          mat_dom.push_back(dromon::MaterialDomain<double>(
          dromon::Material<double>(1.0, 1.0, true, false),
          dromon::Material<double>(1.0, 1.0, false, false)));

          mat_dom1.push_back(dromon::MaterialDomain<double>(
          dromon::Material<double>(1.0, 1.0, true, false),
          dromon::Material<double>(1.0, 1.0, false, false)));

          sidelength = normalDist[i];
          saveName = "./output/MC_mu1_sig0p1_plate.txt";

          break;
        }


      }
      wave = 2.99792458e8 / freq;
      Mesh<2, 3, CUBICP> mesh;
      Point<3, double> center = {0.0, 0.0, 0.0};

      MeshGenerator::square_plate(mesh, center, sidelength,
                                      n_cells_per_dim);

      // MeshGenerator::hyper_sphere(mesh, center, sidelength,
      //                                     n_cells_per_dim);

      Mesh<2,3, CUBICP> mesh1;
      
      if (perturbVar == 3)
        sidelength *= perturbSize;

      // MeshGenerator::hyper_sphere(mesh1, center, sidelength,
      //                                     n_cells_per_dim);

      MeshGenerator::square_plate(mesh1, center, sidelength,
                                n_cells_per_dim);

      Problems::AdaptiveSolver solver(&mesh, &mesh1, &mat_dom, &mat_dom1, 1, 10 ,5,5,5,5);
      solver.set_plane_wave_excitation(freq, {1.0,0.0});



      if (perturbVar == 3)
        sidelength /= perturbSize;


      std::cout << std::setprecision(12);

      double theta_sc = constants<double>::PI / 2.0;
      double phi_sc = 0.0;
      double R_dist_scalar = 100.0;
      // R_dist_scalar *= wave;
      solver.set_scattering_parameters(theta_sc,phi_sc,R_dist_scalar,{0.0,0.0,1.0});

      auto t1 = std::chrono::high_resolution_clock::now();

      out1 = solver.execute_refinement(1, 0.01, forward_matrix0, forward_excitation0, forward_solution0, adjoint_solution0, gradient0, gradient2_0, sidelength, 1, perturbVar, perturbSize);

      auto t2 = std::chrono::high_resolution_clock::now();

      auto exec_time = std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1);

      std::cout << "total execution time: " << exec_time.count() << std::endl;
      std::cout << "variable-----------------: " << normalDist[i] << std::endl;

      //  out2.push_back(output);
      std::complex<double> RCS = 4.0*3.14159*R_dist_scalar*R_dist_scalar*out1*std::conj(out1);
      std::ofstream outFile(saveName, std::ios_base::app);
      outFile << normalDist[i] << "\t" << wave << "\t" << RCS << "\t" << out1 << "\t" << gradient0 << "\t" << gradient2_0 << "\n";


    }
    break;
  }
  case 2:
  {
    
    double sidelength;
    double mat;

    for (int i = 0; i < 2; ++i) 
    {
      MaterialData<double> mat_dom;
      MaterialData<double> mat_dom1;
      std::cout << "New Run _________________________________________________________\n";

      switch (perturbVar)
      {
        case 1:
        {

          // Make the MaterialData
          // In this case, a PEC object embedded in air
          //PEC
          // mat_dom.clear();
          mat = 1.0;
          mat += i*mat*(perturbSize - 1.0);
          mat_dom.push_back(dromon::MaterialDomain<double>(
          dromon::Material<double>(1.0, 1.0, true, false),
          dromon::Material<double>(mat, 1.0, false, false)));

          
          mat_dom1.push_back(dromon::MaterialDomain<double>(
          dromon::Material<double>(1.0, 1.0, true, false),
          dromon::Material<double>(mat, 1.0, false, false)));

          std::cout << "perturbing Material, eps = " << mat << std::endl;
          sidelength = 1.0;// / wave;
          break;
        }

        case 2:
        {
          
          // Make the MaterialData
          // In this case, a PEC object embedded in air
          //PEC
          mat_dom.push_back(dromon::MaterialDomain<double>(
          dromon::Material<double>(1.0, 1.0, true, false),
          dromon::Material<double>(1.0, 1.0, false, false)));

          mat_dom1.push_back(dromon::MaterialDomain<double>(
          dromon::Material<double>(1.0, 1.0, true, false),
          dromon::Material<double>(1.0, 1.0, false, false)));

          sidelength = 1.0;// / wave;
          freq = freq + i*freq*(perturbSize-1.0);
          std::cout << "perturbing freq, f = " << freq << std::endl;

          break;

        }

        case 3:
        {

          // Make the MaterialData
          // In this case, a PEC object embedded in air
          //PEC
          sidelength = 4.0;

          mat_dom.push_back(dromon::MaterialDomain<double>(
          dromon::Material<double>(1.0, 1.0, true, false),
          dromon::Material<double>(1.0, 1.0, false, false)));

          mat_dom1.push_back(dromon::MaterialDomain<double>(
          dromon::Material<double>(1.0, 1.0, true, false),
          dromon::Material<double>(1.0, 1.0, false, false)));

          sidelength += i*(perturbSize-1.0)*sidelength;
          std::cout << "perturbing radius, r = " << sidelength << std::endl;
          break;

        }


      }

      wave = 2.99792458e8 / freq;

      Mesh<2, 3, CUBICP> mesh;
      Point<3, double> center = {0.0, 0.0, 0.0};

      // MeshGenerator::hyper_sphere(mesh, center, sidelength,
      //                                     n_cells_per_dim);

      MeshGenerator::square_plate(mesh, center, sidelength,
                                n_cells_per_dim);

      Mesh<2,3, CUBICP> mesh1;

      if (perturbVar == 3)
        sidelength *= perturbSize;

      // MeshGenerator::hyper_sphere(mesh1, center, sidelength,
      //                                     n_cells_per_dim);

      MeshGenerator::square_plate(mesh1, center, sidelength,
                                n_cells_per_dim);    

      // for (int i = 0; i < mesh.n_nodes(); ++i)
      // std::cout << "mesh nodes: " << *mesh.get_node(i) << std::endl;

      if (perturbVar == 3)
        sidelength /= perturbSize;

      // double freq = 285e6;
      Problems::AdaptiveSolver solver(&mesh, &mesh1, &mat_dom, &mat_dom1, 1, 10 ,5,5,5,5);
      solver.set_plane_wave_excitation(freq, {1.0,0.0});
      std::cout << std::setprecision(12);
      double theta_sc = constants<double>::PI / 2.0;
      double phi_sc = 0.0;
      double R_dist_scalar = 1000.0;
      // R_dist_scalar *= wave;
      solver.set_scattering_parameters(theta_sc,phi_sc,R_dist_scalar,{0.0,0.0,-1.0});

      auto t1 = std::chrono::high_resolution_clock::now();

      if (i == 0)
        out1 = solver.execute_refinement(1, 0.01, forward_matrix0, forward_excitation0, forward_solution0, adjoint_solution0, gradient0, gradient2_0, sidelength, 1, perturbVar, perturbSize);
      if (i == 1)
        out2 = solver.execute_refinement(1, 0.01, forward_matrix1, forward_excitation1, forward_solution1, adjoint_solution1, gradient1, gradient2_1, sidelength, 1, perturbVar, perturbSize);

      auto t2 = std::chrono::high_resolution_clock::now();
      auto exec_time = std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1);
      std::cout << "total execution time: " << exec_time.count() << std::endl;

      std::cout << "wavelength: " << wave << std::endl;
      std::complex<double> RCS = 4.0*3.14159*R_dist_scalar*R_dist_scalar*out1*std::conj(out1);
    }

    int m = forward_excitation0.size();
    std::vector<std::complex<double>> forward_excitation_grad (m);
    std::vector<std::complex<double>> forward_matrix_grad (m*m);

    if (perturbVar == 1){

      for (int i = 0; i < m*m; ++i){
        forward_matrix_grad[i] = (forward_matrix1[i] - forward_matrix0[i])/((perturbSize - 1.0)*mat*8.8541878128e-12);//*8.8541878128e-12);
      }
      for (int i = 0; i < m; ++i){
        forward_excitation_grad[i] = (forward_excitation1[i] - forward_excitation0[i])/((perturbSize - 1.0)*mat*8.8541878128e-12);///(0.01*2.0);//*8.8541878128e-12);
      }
      std::cout << "FD division factor: " << ((perturbSize - 1.0)*mat*8.8541878128e-12) << std::endl;
      
    }
    else if (perturbVar == 2){

      for (int i = 0; i < m*m; ++i){
        forward_matrix_grad[i] = (forward_matrix1[i] - forward_matrix0[i])/((perturbSize - 1.0)*freq);///(0.01*2.0);//*8.8541878128e-12);
      }
      for (int i = 0; i < m; ++i){
        forward_excitation_grad[i] = (forward_excitation1[i] - forward_excitation0[i])/((perturbSize - 1.0)*freq);///(0.01*2.0);//*8.8541878128e-12);
      }
      std::cout << "FD division factor: " << ((perturbSize - 1.0)*freq) << std::endl;


    }
    else if (perturbVar == 3){

      for (int i = 0; i < m*m; ++i){
        forward_matrix_grad[i] = (forward_matrix1[i] - forward_matrix0[i])/((perturbSize - 1.0)*2.);///(0.01*2.0);//*8.8541878128e-12);
      }
      for (int i = 0; i < m; ++i){
        forward_excitation_grad[i] = (forward_excitation1[i] - forward_excitation0[i])/((perturbSize - 1.0)*2.);///(0.01*2.0);//*8.8541878128e-12);
      }
      std::cout << "FD division factor: " << ((perturbSize - 1.0)*sidelength) << std::endl;

    }

    std::vector<std::complex<double>> forward_system_grad (m);
    for (int i = 0; i < m; ++i){
      for (int j = 0; j < m; ++j){
        forward_system_grad[i] += forward_matrix_grad[i*m+j]*forward_solution0[j];
      }
    }

    std::complex<double> finiteGrad(0.0);

    for (unsigned int i =0; i < m; ++i)
      {
      //std::cout << "excitation: " << this->forward_excitation_HO_HOPS[i] << "   system: " << (this->forward_system_HO_HOPS[i]) << std::endl;
      finiteGrad += (forward_excitation_grad[i] - forward_system_grad[i])*std::conj(adjoint_solution0[i]);
      // std::cout << "deltaG: " << forward_excitation_grad[i] << "   deltaL*u: " << forward_system_grad[i] << std::endl;

    }

    

    if (perturbVar == 1){
      std::cout << "division: " << ((perturbSize - 1.0)*mat) << "  out2: " << out2 << "  out1: " << out1 << std::endl;
      std::cout << "Gradient from FD of G and L!!!:  " << finiteGrad << std::endl;
      std::cout << "Grad from FD of output: " << (out2 - out1)/((perturbSize - 1.0)*mat*8.8541878128e-12)<< std::endl;
      std::cout << "abs compare: " << abs((out2 - out1)/((perturbSize - 1.0)*mat*8.8541878128e-12)) << ", " << abs(finiteGrad) << std::endl;

    }
    else if (perturbVar == 2){
      std::cout << "Gradient from FD of G and L!!!:  " << finiteGrad << std::endl;
      std::cout << "Grad from FD of output: " << (out2 - out1)/((perturbSize - 1.0)*freq)<< std::endl;
      std::cout << "abs compare: " << abs((out2 - out1)/((perturbSize - 1.0)*freq)) << ", " << abs(finiteGrad) << std::endl;

    }
    else if (perturbVar == 3){
      std::cout << "Gradient from FD of G and L!!!:  " << finiteGrad << std::endl;
      std::cout << "Grad from FD of output: " << (out2 - out1)/((perturbSize - 1.0)*sidelength)<< std::endl;
      std::cout << "abs compare: " << abs((out2 - out1)/((perturbSize - 1.0)*sidelength)) << ", " << abs(finiteGrad) << std::endl;

    }


    break;
  }
}

  return 0;
}


  





    //     for (int i = 0; i < m; ++i){
    //   std::cout << forward_excitation_grad[i] << std::endl;
    // }

    // for (int i = 0; i < m; ++i){
    //   std::cout << "excitation: " << adjoint_solution0[i] << ", " << adjoint_solution1[i] << std::endl;
    // }

    // ///// //this version takes FD of LJ1 - LJ0
    // std::vector<std::complex<double>> forward_system_grad (m);
    // std::vector<std::complex<double>> forward_system0 (m);
    // std::vector<std::complex<double>> forward_system1 (m);
    // for (int i = 0; i < m; ++i){
    //   for (int j = 0; j < m; ++j){
    //     forward_system0[i] += forward_matrix0[i*m+j]*forward_solution0[j];
    //     forward_system1[i] += forward_matrix1[i*m+j]*forward_solution1[j];

    //   }
    // }

    // for (int i = 0; i < m; ++i){
    //   forward_system_grad[i] = (forward_system1[i] - forward_system0[i]) / (2*3.14159*10);
    // }
    



    // for (int i = 0; i < m; ++i){
    //   finiteGrad += forward_excitation_grad[i] * std::conj(adjoint_solution0[i]);
    //   for (int j = 0; j < m; ++j){
    //       finiteGrad += forward_matrix_grad[i*m+j]*forward_solution0[i] * std::conj(adjoint_solution0[j]);
    //   }
    // }




    // std::cout << "updated QoI: " << out1 + finiteGrad*(1000.0*2*3.14159) << std::endl;

  // for (int i = 0; i < out1.size(); ++i){
  //     std::cout << (out2[i] - out1[i]) / (2*3.14159 * 10) << std::endl;
  // }

  // std::cout << "FD gradient: " << std::endl;
  // for (int i = 0; i < out1.size(); ++i){
  //   std::cout << (out2[i] - out1[i])/(2*3.14159*10) << std::endl;
  // }

  // ///// this is finite differences calculation (grad G or grad L_EE J_S depending on the execute_refinement return)
  // std::ofstream outfileGrad;
  // outfileGrad.open("finiteGradients.txt", std::ios_base::app);
  // outfileGrad << "Grad from FD of output: \n";
  // for (int i = 0; i < 200; ++i){
  //   outfileGrad << (out2[i] - out1[i])/(2*3.14159*10) << std::endl;
  // }
  
  // std::ofstream outfile;
  // outfile.open("forwardL.txt", std::ios_base::app);
  // outfile << " grad from finite differences: \n";
  // for (int i = 0; i < out1.size(); ++i){
  //   outfile << out2[i] << std::endl;

  // }




      //  if (i == 0)
      //  out1 = solver.execute_refinement(1, 0.01, forward_matrix0, forward_excitation0, forward_solution0, adjoint_solution0);
      //  if (i == 1)
      //  out2 = solver.execute_refinement(1, 0.01, forward_matrix1, forward_excitation1, forward_solution1, adjoint_solution1);