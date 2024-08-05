//
// Created by Jake J. Harmon (jake.harmon@ieee.org) on 8/9/21.
//

#ifndef DROMON_UTILITY_H
#define DROMON_UTILITY_H
#include "config.h"

DROMON_NAMESPACE_OPEN

// TODO: Memoization here...
namespace utility
{
    template <class Real>
    void make_real(std::vector<std::complex<Real>>& vec)
    {
      for (auto& entry : vec)
        entry.imag(0);
    }
    template <class Real>
    Point<3, Real> sph2cart(const Real& theta, const Real& phi, const Real& R)
    {
      Point<3, Real> out = {R*std::cos(phi)*std::sin(theta),
                            R*std::sin(phi)*std::sin(theta),
                            R*std::cos(theta)};
      return out;
    }

    template <class Real>
    Real lagrange_p(const unsigned int& K, const unsigned int &m, const Real& u);

    template <class Real>
    Real lagrange_p_prime(const unsigned int& K, const unsigned int &m, const Real& u);

    template<class Real>
    Real lagrange_p(const unsigned int &K, const unsigned int &m, const Real &u) {
        Real out=(Real)(1);
        const Real diff = (Real)(2)/((Real)(K));
        const Real um = (Real)(-1) + diff*m;
        for (unsigned int i = 0; i <= K; ++i)
        {
            if (i==m)
                continue;
            const Real ui = (Real)(-1) + diff*i;
            out *= (u-ui)/(um-ui);
        }
        return out;
    }

    template<class Real>
    Real lagrange_p_prime(const unsigned int &K, const unsigned int &m, const Real &u) {
      Real out = (Real)(0);
      const Real diff = (Real)(2) / ((Real)(K));
      const Real um = (Real)(-1) + diff * m;

      for (unsigned int j = 0; j <= K; ++j) {
        if (j == m)
          continue;
        Real uj = (Real)(-1) + diff * j;
        Real prod = (Real)(1);
        for (unsigned int i = 0; i <= K; ++i) {
          if (i == m || i == j)
            continue;
          const Real ui = (Real)(-1) + diff * i;
          prod *= (u - ui) / (um - ui);
        }
        out += prod/(um - uj);
      }
      return out;
    }

    template<class Fval, class Fderiv, class Real>
    Real newton_raphson(Fval fval, Fderiv fderiv, Real guess, Real abs_tolerance)
    {
      Real out = guess;
      Real prev = guess;
      Real current_error = std::numeric_limits<Real>::max();
      while (current_error > abs_tolerance)
      {
        out -= fval(out)/fderiv(out);
        current_error = std::abs(out-prev);
        prev = out;
      }
      return out;
    }

}

DROMON_NAMESPACE_CLOSE
#endif //DROMON_UTILITY_H
