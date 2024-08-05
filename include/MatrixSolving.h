//
// Created by Jake J. Harmon (jake.harmon@ieee.org) on 4/27/22.
//

#ifndef DROMON_MATRIXSOLVING_H
#define DROMON_MATRIXSOLVING_H
#include "config.h"

DROMON_NAMESPACE_OPEN
namespace MatrixSolving {
template <class CoefficientType, class MatrixSolvingPolicy>
class Matrix : private MatrixSolvingPolicy {
public:
  Matrix(const unsigned int &m_rows, const unsigned int& n_cols);

  using MatrixSolvingPolicy::solve_system;
  // At some point, when C++ concepts are supported, make this function based
  // on that
  template <class GalerkinSystem>
  void matrix_from_galerkin_system(GalerkinSystem *galerkinSystem, typename GalerkinSystem::DoFHandlerType* dof_handler);

  template <class GalerkinSystem>
  void RHS_from_galerkin_system(GalerkinSystem* galerkinSystem, typename GalerkinSystem::DoFHandlerType* dof_handler);

  template <class GalerkinSystem>
  void adjoint_RHS_from_galerkin_system(GalerkinSystem* galerkinSystem, typename GalerkinSystem::DoFHandlerType* dof_handler);

  void solve_forward(const bool& estimate_condition_number = false);
  void solve_adjoint(const bool& estimate_condition_number = false);
  // void solve_adjoint();
  typename norm_type<CoefficientType>::type norm_1() const;
  void purge();
  void zero_RHS();
  std::vector<CoefficientType> get_solution_copy();
  std::vector<CoefficientType> get_RHS_copy();
  unsigned int get_n_rows() const;
  unsigned int get_n_cols() const;
  CoefficientType &operator()(const unsigned int &i, const unsigned int &j);
  CoefficientType &operator()(const unsigned int &i);

  CoefficientType at(const unsigned int &i, const unsigned int &j) const;
  CoefficientType at(const unsigned int &i) const;

private:
  // Is_safe indicates whether the entries are safe to modify/access
  // After a solve, is_safe is false, and the RHS will be the solution in the
  // case of using LAPACK for matrix solving
  bool is_safe = true;
  const bool is_row_major = true;
  const unsigned int m;
  const unsigned int n;
  const unsigned int size;

  // A vector containing the actual entries of the left hand side matrix
  // which contains m rows and n columns ...
  // As part of the matrix solving procedure, the entries may be modified
  // std::vector<CoefficientType> LHS_data;
  std::unique_ptr<CoefficientType[]> LHS_data;

  // A vector containing the actual entries of the right hand side vector
  // which contains m rows
  // As part of the matrix solving procedure, the RHS_data may be replaced
  // by the solution data
  // std::vector<CoefficientType> RHS_data;
  std::unique_ptr<CoefficientType[]> RHS_data;




};
template <class CoefficientType, class MatrixSolvingPolicy>
void Matrix<CoefficientType, MatrixSolvingPolicy>::purge() {
  //  LHS_data.clear();
  //  RHS_data.clear();
  LHS_data.reset();
  RHS_data.reset();
}
template <class CoefficientType, class MatrixSolvingPolicy>
template <class GalerkinSystem>
void Matrix<CoefficientType, MatrixSolvingPolicy>::matrix_from_galerkin_system(
    GalerkinSystem *galerkinSystem, typename GalerkinSystem::DoFHandlerType* dof_handler) {
  // Loop through the subcells and subexcitations in the GalerkinSystem and
  // assemble the full global system
  const auto n_cells = galerkinSystem->size();

  // Get the DoFHandler attached to the GalerkinSystem
  //auto dof_handler = galerkinSystem->get_dof_handler();

  const bool is_symmetric = galerkinSystem->is_symmetric_system();
  //#ifndef DEBUG
  //#pragma omp parallel for
  //#endif
  for (unsigned int i = 0; i < n_cells; ++i) {
    // Get the subsystem and fill the global
    const auto& subsystem_test = galerkinSystem->subsystem_at(i);
    // const auto& subexcitation_test = galerkinSystem->subexcitation_at(i);

    // Loop through all the DoFs on this given cell
    const auto& dof_test_parent = dof_handler->get_dof_parent(i);
    for (unsigned int dof_test_index = 0;
         dof_test_index < dof_test_parent.n_dofs(); ++dof_test_index) {
      const auto &dof_test = dof_test_parent.get_dof(dof_test_index);
      const unsigned int &global_test_index = dof_test.global_index;
      const int &active_test_index = dof_test.active_index;
      if (!dof_test.is_active)
        continue;
      //      this->operator()(active_test_index) +=
      //          subexcitation_test.const_at(global_test_index);

      for (unsigned int j = 0; j < n_cells; ++j) {
        const auto& subsystem_trial = galerkinSystem->subsystem_at(j);
        //        if (is_symmetric && j > i)
        //          break;
        const auto& dof_trial_parent = dof_handler->get_dof_parent(j);

        for (unsigned int dof_trial_index = 0;
             dof_trial_index < dof_trial_parent.n_dofs(); ++dof_trial_index) {
          const auto &dof_trial = dof_trial_parent.get_dof(dof_trial_index);
          const unsigned int &global_trial_index = dof_trial.global_index;
          const int &active_trial_index = dof_trial.active_index;
          if (!dof_trial.is_active)
            continue;
          //          if (is_symmetric && active_trial_index > active_test_index && j < i)
          //            std::cout << "This situation does happen!" << std::endl;
          if (is_symmetric && active_trial_index > active_test_index)
            continue;
          //          if((active_trial_index == 0 && active_test_index == 7) || (active_trial_index == 7 && active_test_index == 0))
          //          {
          //            std::cout << "Made it here!" << std::endl;
          //          }
          if (active_trial_index == 4 && active_test_index == 6)
          {
            int pause_here =0;
          }
          if (is_symmetric && j < i)
          {
            //            #ifndef DEBUG
            //            #pragma omp critical
            //            #endif
            //  #pragma omp critical
            {
              auto temp1 =subsystem_test.const_at(global_test_index,
                                                   global_trial_index);
              auto temp2 = subsystem_trial.const_at(global_trial_index,
                                                    global_test_index);
              this->operator()(active_test_index, active_trial_index) +=
                  subsystem_trial.const_at(global_trial_index,
                                           global_test_index);
            }
          }
          else {
            //            #ifndef DEBUG
            //            #pragma omp critical
            //            #endif
            //   #pragma omp critical
            {
              auto temp1 =subsystem_test.const_at(global_test_index,
                                                   global_trial_index);
              auto temp2 = subsystem_trial.const_at(global_trial_index,
                                                    global_test_index);
              this->operator()(active_test_index, active_trial_index) +=
                  subsystem_test.const_at(global_test_index,
                                          global_trial_index);
            }
          }

        }
      }
    }
  }

  // Now that we filled the system and excitation matrices, we need to, if
  // the system is symmetric, fill the upper diagonal
  if (is_symmetric)
    for (unsigned int i = 0; i < m; ++i)
      for (unsigned int j = i + 1; j < n; ++j)
        this->operator()(i, j) = this->operator()(j, i);
}
template <class CoefficientType, class MatrixSolvingPolicy>
CoefficientType &Matrix<CoefficientType, MatrixSolvingPolicy>::operator()(
    const unsigned int &i, const unsigned int &j) {

  return LHS_data[this->n * i + j];
}
template <class CoefficientType, class MatrixSolvingPolicy>
CoefficientType &Matrix<CoefficientType, MatrixSolvingPolicy>::operator()(
    const unsigned int &i) {
  return RHS_data[i];
}
template <class CoefficientType, class MatrixSolvingPolicy>
CoefficientType Matrix<CoefficientType, MatrixSolvingPolicy>::at(
    const unsigned int &i, const unsigned int &j) const {

  return LHS_data[this->n * i + j];
}
template <class CoefficientType, class MatrixSolvingPolicy>
std::vector<CoefficientType>
Matrix<CoefficientType, MatrixSolvingPolicy>::get_solution_copy() {
  assert(!is_safe && "The system has not been solved, therefore the solution cannot be returned!");
  //return *RHS_data;
  std::vector<CoefficientType> out(this->m);
  for (unsigned int i = 0; i < this->m; ++i)
    out[i] = RHS_data[i];

  return out;
}
template <class CoefficientType, class MatrixSolvingPolicy>
std::vector<CoefficientType>
Matrix<CoefficientType, MatrixSolvingPolicy>::get_RHS_copy() {
  assert(is_safe && "The system has been solved, therefore the RHS cannot be returned!");
  //return *RHS_data;
  std::vector<CoefficientType> out(this->m);
  for (unsigned int i = 0; i < this->m; ++i)
    out[i] = RHS_data[i];

  return out;
}
template <class CoefficientType, class MatrixSolvingPolicy>
void Matrix<CoefficientType, MatrixSolvingPolicy>::solve_forward(const bool& estimate_condition_number)
{
  if (!estimate_condition_number)
    solve_system(is_row_major, m, n, LHS_data.get(), RHS_data.get(), false);
  else
    solve_system(is_row_major, m, n, LHS_data.get(), RHS_data.get(), false, this->norm_1());
  is_safe = false;
}

template <class CoefficientType, class MatrixSolvingPolicy>
void Matrix<CoefficientType, MatrixSolvingPolicy>::solve_adjoint(
    const bool &estimate_condition_number)
{
  if (!estimate_condition_number)
    solve_system(is_row_major, m, n, LHS_data.get(), RHS_data.get(), true);
  else
    solve_system(is_row_major, m, n, LHS_data.get(), RHS_data.get(), true, this->norm_1());
  is_safe = false;
}

template <class CoefficientType, class MatrixSolvingPolicy>
Matrix<CoefficientType, MatrixSolvingPolicy>::Matrix(const unsigned int& m_rows, const unsigned int& n_cols) : m(m_rows), n(n_cols), size(m_rows*n_cols)
{
  //  LHS_data.resize(m*n, CoefficientType(0.0));
  //  RHS_data.resize(m, CoefficientType(0.0));
  LHS_data.reset(new CoefficientType[m*n]);
  RHS_data.reset(new CoefficientType[m]);

  std::fill_n((LHS_data.get()), size, CoefficientType(0));
  std::fill_n((RHS_data.get()), m, CoefficientType(0));
}
template <class CoefficientType, class MatrixSolvingPolicy>
unsigned int Matrix<CoefficientType, MatrixSolvingPolicy>::get_n_rows() const {
  return m;
}
template <class CoefficientType, class MatrixSolvingPolicy>
unsigned int Matrix<CoefficientType, MatrixSolvingPolicy>::get_n_cols() const {
  return n;
}
template <class CoefficientType, class MatrixSolvingPolicy>
CoefficientType
Matrix<CoefficientType, MatrixSolvingPolicy>::at(const unsigned int &i) const {
  return this->RHS_data[i];
}
template <class CoefficientType, class MatrixSolvingPolicy>
typename norm_type<CoefficientType>::type Matrix<CoefficientType, MatrixSolvingPolicy>::norm_1() const {
  typename norm_type<CoefficientType>::type out(0);

  for (unsigned int j = 0; j < this->get_n_cols(); ++j)
  {
    typename norm_type<CoefficientType>::type col_sum(0);
    for (unsigned int i = 0; i < this->get_n_rows(); ++i)
    {
      col_sum += std::abs(this->at(i,j));
    }
    out = std::max(out, col_sum);
  }
  return out;
}
template <class CoefficientType, class MatrixSolvingPolicy>
void Matrix<CoefficientType, MatrixSolvingPolicy>::zero_RHS() {
  std::fill_n((RHS_data.get()), m, CoefficientType(0));
}

template <class CoefficientType, class MatrixSolvingPolicy>
template <class GalerkinSystem>
void Matrix<CoefficientType, MatrixSolvingPolicy>::
    RHS_from_galerkin_system(GalerkinSystem *galerkinSystem, typename GalerkinSystem::DoFHandlerType* dof_handler)
{
  // Loop through the subcells and subexcitations in the GalerkinSystem and
  // assemble the full global system
  const auto n_cells = galerkinSystem->size();

  // Get the DoFHandler attached to the GalerkinSystem
  // auto dof_handler = galerkinSystem->get_dof_handler();

  //#ifndef DEBUG
  //#pragma omp parallel for
  //#endif
  for (unsigned int i = 0; i < n_cells; ++i) {
    // Get the subsystem and fill the global

    const auto& subexcitation_test = galerkinSystem->subexcitation_at(i);

    // Loop through all the DoFs on this given cell
    const auto& dof_test_parent = dof_handler->get_dof_parent(i);
    const auto n_dofs_temp = dof_test_parent.n_dofs();
    for (unsigned int dof_test_index = 0;
         dof_test_index < dof_test_parent.n_dofs(); ++dof_test_index) {
      const auto &dof_test = dof_test_parent.get_dof(dof_test_index);
      const unsigned int &global_test_index = dof_test.global_index;
      const int &active_test_index = dof_test.active_index;
      if (!dof_test.is_active)
        continue;

      //      #ifndef DEBUG
      //      #pragma omp critical
      //      #endif
      this->operator()(active_test_index) +=
          subexcitation_test.const_at(global_test_index);
    }
  }
}

template <class CoefficientType, class MatrixSolvingPolicy>
template <class GalerkinSystem>
void Matrix<CoefficientType, MatrixSolvingPolicy>::
    adjoint_RHS_from_galerkin_system(GalerkinSystem *galerkinSystem, typename GalerkinSystem::DoFHandlerType* dof_handler)
{
  // Loop through the subcells and subexcitations in the GalerkinSystem and
  // assemble the full global system
  const auto n_cells = galerkinSystem->size();

  // Get the DoFHandler attached to the GalerkinSystem
  // auto dof_handler = galerkinSystem->get_dof_handler();

  //#ifndef DEBUG
  //#pragma omp parallel for
  //#endif
  for (unsigned int i = 0; i < n_cells; ++i) {
    // Get the subsystem and fill the global

    const auto& subexcitation_test = galerkinSystem->adjoint_subexcitation_at(i);

    // Loop through all the DoFs on this given cell
    const auto& dof_test_parent = dof_handler->get_dof_parent(i);
    for (unsigned int dof_test_index = 0;
         dof_test_index < dof_test_parent.n_dofs(); ++dof_test_index) {
      const auto &dof_test = dof_test_parent.get_dof(dof_test_index);
      const unsigned int &global_test_index = dof_test.global_index;
      const int &active_test_index = dof_test.active_index;
      if (!dof_test.is_active)
        continue;

      //      #ifndef DEBUG
      //      #pragma omp critical
      //      #endif
      this->operator()(active_test_index) +=
          subexcitation_test.const_at(global_test_index);
    }
  }
}


template <class CoefficientType, class MatrixSolvingPolicy>
inline std::ostream& operator<<(std::ostream& os, const Matrix<CoefficientType, MatrixSolvingPolicy>& mat)
{
  for (unsigned int i = 0; i < mat.get_n_rows(); ++i)
  {
    for (unsigned int j = 0; j < mat.get_n_cols(); ++j)
    {
      const auto& val = mat.at(i,j);
      if (std::is_same<CoefficientType, std::complex<double>>::value || std::is_same<CoefficientType, std::complex<float>>::value)
        os << val.real() << " " << val.imag()  << std::endl;
    }
  }
  return os;

}


} // namespace MatrixSolving

DROMON_NAMESPACE_CLOSE

#endif // DROMON_MATRIXSOLVING_H
