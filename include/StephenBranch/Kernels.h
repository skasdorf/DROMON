//
// Created by Jake J. Harmon (jake.harmon@ieee.org) on 4/10/22.
//

#ifndef DROMON_KERNELS_H
#define DROMON_KERNELS_H
#include "config.h"
#include <complex>
#include "Point.h"

DROMON_NAMESPACE_OPEN
namespace kernels
{
  template<class Real = double, class MatValueType>
  std::complex<Real> gScaled(const Real& R, const Real& omega, const MatValueType& epsr, const MatValueType& mur, const double wavelength)
  {
    auto tester = constants<Real>::ROOT_EPS0MU0_;
    const Real Rscaled = R*wavelength;
    const std::complex<Real> gamma = std::complex<Real>(0.0,1.0)*omega*sqrt(epsr*mur)*constants<Real>::ROOT_EPS0MU0_;
    const auto test = std::exp(-gamma*R);
    return std::exp(-gamma*Rscaled)/(Real(4.0)*constants<Real>::PI*Rscaled);
  }

  template<class Real = double, class MatValueType>
  std::complex<Real> g(const Real& R, const Real& omega, const MatValueType& epsr, const MatValueType& mur)
  {
    auto tester = constants<Real>::ROOT_EPS0MU0_;

    const std::complex<Real> gamma = std::complex<Real>(0.0,1.0)*omega*sqrt(epsr*mur)*constants<Real>::ROOT_EPS0MU0_;
    const auto test = std::exp(-gamma*R);
    return std::exp(-gamma*R)/(Real(4.0)*constants<Real>::PI*R);
  }

  template<class Real = double, class MatValueType>
  std::complex<Real> gPrint(const Real& R, const Real& omega, const MatValueType& epsr, const MatValueType& mur)
  {
    auto tester = constants<Real>::ROOT_EPS0MU0_;
    const std::complex<Real> gamma = std::complex<Real>(0.0,1.0)*omega*sqrt(epsr*mur)*constants<Real>::ROOT_EPS0MU0_;
    const auto test = std::exp(-gamma*R);

    std::cout << "gamma: " << gamma << "  R: " << R << "   phase: " << std::exp(-gamma*R) << "   denom: " << Real(4.0)*constants<Real>::PI*R << std::endl;

    return std::exp(-gamma*R)/(Real(4.0)*constants<Real>::PI*R);
  }

  // template<class Real = double, class MatValueType>
  // std::complex<Real> dgdw(const Real& R, const Real& omega, const MatValueType& epsr, const MatValueType& mur)
  // {
  //   auto tester = constants<Real>::ROOT_EPS0MU0_;

  //   const std::complex<Real> gamma = std::complex<Real>(0.0,1.0)*omega*sqrt(epsr*mur)*constants<Real>::ROOT_EPS0MU0_;
  //   const auto test = std::exp(-gamma*R);
  //   return -gamma/omega*R;
  // }

  template <unsigned int spacedim, class Real = double, class MatValueType>
  Point<spacedim, std::complex<Real>> grad_gScaled(const Point<spacedim, Real>& Rvec, const Real& omega, const MatValueType& epsr, const MatValueType& mur, const double wavelength)
  {
    std::complex<Real> gamma = std::complex<Real>(0.,1.)*omega*sqrt(epsr*mur)*constants<Real>::ROOT_EPS0MU0_;
    const Real R = Rvec.norm();
    const Real Rscaled = R*wavelength;
    return (-gamma*gamma/(4*constants<Real>::PI*Rscaled)*std::exp(-1.*gamma*Rscaled)*(1./(gamma*Rscaled)+1./(gamma*gamma*Rscaled*Rscaled)))*Rvec*wavelength;
  }

  // template <unsigned int spacedim, class Real = double, class MatValueType>
  // Point<spacedim, std::complex<Real>> grad_g(const Point<spacedim, Real>& Rvec, const Real& omega, const MatValueType& epsr, const MatValueType& mur)
  // {
  //   std::complex<Real> gamma = std::complex<Real>(0.,1.)*omega*sqrt(epsr*mur)*constants<Real>::ROOT_EPS0MU0_;
  //   const Real R = Rvec.norm();
  //   return (-gamma*gamma/(4*constants<Real>::PI*R)*std::exp(-1.*gamma*R)*(1./(gamma*R)+1./(gamma*gamma*R*R)))*Rvec;
  // }

  template <unsigned int spacedim, class Real = double, class MatValueType>
  Point<spacedim, std::complex<Real>> grad_g(const Point<spacedim, Real>& Rvec, const Real& R, const Real& omega, const MatValueType& epsr, const MatValueType& mur)
  {
    std::complex<Real> gamma = std::complex<Real>(0.,1.)*omega*sqrt(epsr*mur)*constants<Real>::ROOT_EPS0MU0_;
    // std::cout << "checking: " <<(-gamma*gamma/(4*constants<Real>::PI*R)*std::exp(-1.*gamma*R)*(1./(gamma*R)+1./(gamma*gamma*R*R)))*Rvec << "  " << (-std::exp()/(4*constants<Real>::PI*R*R*R)*(gamma*R+1)*Rvec)
    return (-gamma*gamma/(4*constants<Real>::PI*R)*std::exp(-1.*gamma*R)*(1./(gamma*R)+1./(gamma*gamma*R*R)))*Rvec;
    // return (-std::exp(-gamma*R)/(4*constants<Real>::PI*R*R*R)*(gamma*R+1.)*Rvec);
  }


  // template <unsigned int spacedim, class Real = double, class MatValueType>
  // Point<spacedim, std::complex<Real>> grad_g(const Point<spacedim, Real>& Rvec, const Real& R, const Real& omega, const MatValueType& epsr, const MatValueType& mur)
  // {
  //   std::complex<Real> gamma = std::complex<Real>(0.,1.)*omega*sqrt(epsr*mur)*constants<Real>::ROOT_EPS0MU0_;
  //   return (-gamma*gamma/(4*constants<Real>::PI*R)*std::exp(-1.*gamma*R)*(1./(gamma*R)+1./(gamma*gamma*R*R)))*Rvec;
  // }
}

DROMON_NAMESPACE_CLOSE

#endif // DROMON_KERNELS_H
