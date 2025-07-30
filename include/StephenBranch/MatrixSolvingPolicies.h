//
// Created by Jake J. Harmon (jake.harmon@ieee.org) on 4/27/22.
//

#ifndef DROMON_MATRIXSOLVINGPOLICIES_H
#define DROMON_MATRIXSOLVINGPOLICIES_H
#include <complex>
#include "config.h"
#include "mkl.h"
#include "mkl_lapack.h"

DROMON_NAMESPACE_OPEN

namespace MatrixSolvingPolicies
{
    class MKLPolicy
    {
    public:
        template <class CoefficientType>
        void solve_system(const bool& is_row_major, const unsigned int& m_rows, const unsigned int& n_cols, std::vector<CoefficientType>& matrix_data, std::vector<CoefficientType>& rhs_data,const bool& solve_adjoint, const double& anorm = 0.0);

        template <class CoefficientType>
        void solve_system(const bool& is_row_major, const unsigned int& m_rows, const unsigned int& n_cols, CoefficientType* matrix_data, CoefficientType* rhs_data,const bool& solve_adjoint, const double& anorm = 0.0);
      private:
        std::vector<lapack_int> ipiv;
        double a_norm;
    };


    template <class CoefficientType>
    void MKLPolicy::solve_system(const bool &is_row_major,
                                 const unsigned int &m_rows,
                                 const unsigned int &n_cols,
                                 std::vector<CoefficientType> &matrix_data,
                                 std::vector<CoefficientType> &rhs_data,
                                 const bool& solve_adjoint,
                                 const double& anorm)
    {
     // lapack_int ipiv [std::min(m_rows, n_cols)];
     if (ipiv.size() == 0) {
       ipiv.resize(std::min(m_rows, n_cols));
     }

      lapack_int lda = n_cols;
      lapack_int ldb = 1;
      lapack_complex_double test;

      std::cout << "is row major?: " << is_row_major << std::endl;

      // Now we have to decide which LAPACK call to use, based on the CoefficientType
      if (std::is_same<CoefficientType, std::complex<double>>::value)
      {
        auto error = LAPACKE_zgetrf(is_row_major ? LAPACK_ROW_MAJOR : LAPACK_COL_MAJOR, int(m_rows), int(n_cols),
                                    &matrix_data[0], lda, &ipiv[0]);
       assert(error == 0 && "Failure during matrix factorization!");

       auto error_solve = LAPACKE_zgetrs(is_row_major ? LAPACK_ROW_MAJOR : LAPACK_COL_MAJOR, solve_adjoint ? 'C' : 'N',
                                         int(m_rows), 1, &matrix_data[0], lda, &ipiv[0], &rhs_data[0], ldb);
       assert(error_solve == 0 && "Failure during matrix solve!");
       if (anorm > 0)
       {
         // Get an estimate of the condition number
         double rcond;
         auto error_cond = LAPACKE_zgecon(is_row_major ? LAPACK_ROW_MAJOR : LAPACK_COL_MAJOR, '1', int(m_rows), &matrix_data[0], lda, anorm, &rcond);
         assert(error_cond == 0 && "Failure during estimating condition number!");
         std::cout << "Estimated condition number: " << rcond << std::endl;
       }
      }
    }
    template <class CoefficientType>
    void MKLPolicy::solve_system(const bool &is_row_major,
                                 const unsigned int &m_rows,
                                 const unsigned int &n_cols,
                                 CoefficientType* matrix_data,
                                 CoefficientType* rhs_data,
                                 const bool& solve_adjoint,
                                 const double& anorm)
    {
      // lapack_int ipiv [std::min(m_rows, n_cols)];
      if (ipiv.size() == 0) {
        ipiv.resize(std::min(m_rows, n_cols));
      }

      lapack_int lda = n_cols;
      lapack_int ldb = 1;
      lapack_complex_double test;



      // Now we have to decide which LAPACK call to use, based on the CoefficientType
      if (std::is_same<CoefficientType, std::complex<double>>::value)
      {
        auto error = LAPACKE_zgetrf(is_row_major ? LAPACK_ROW_MAJOR : LAPACK_COL_MAJOR, int(m_rows), int(n_cols),
                                    &matrix_data[0], lda, &ipiv[0]);
        assert(error == 0 && "Failure during matrix factorization!");

        auto error_solve = LAPACKE_zgetrs(is_row_major ? LAPACK_ROW_MAJOR : LAPACK_COL_MAJOR, solve_adjoint ? 'C' : 'N',
                                          int(m_rows), 1, matrix_data, lda, &ipiv[0], rhs_data, ldb);
        assert(error_solve == 0 && "Failure during matrix solve!");

        if (anorm > 0)
        {
          // Get an estimate of the condition number
          double rcond;
          auto error_cond = LAPACKE_zgecon(is_row_major ? LAPACK_ROW_MAJOR : LAPACK_COL_MAJOR, '1', int(m_rows), matrix_data, lda, anorm, &rcond);
          assert(error_cond == 0 && "Failure during estimating condition number!");
          std::cout << "Estimated condition number: " << 1./rcond << std::endl;
        }

      }
    }

    }

DROMON_NAMESPACE_CLOSE
#endif // DROMON_MATRIXSOLVINGPOLICIES_H
