//
// Created by Jake J. Harmon (jake.harmon@ieee.org) on 4/10/22.
// Edited by Christopher A. Erickson (Christopher.Erickson@ieee.org) on 7/2/25
//

#ifndef DROMON_KERNELS_H //checks whether DROMON_KERNELS has not been defined
#define DROMON_KERNELS_H //if it hasnt been defined define it now, this prevents multiple inclusions of this header file
#include "config.h"
#include <complex> //Standard C++ header providing std::complex<T> (complex numbers).
#include "Point.h"

DROMON_NAMESPACE_OPEN
namespace kernels
{
    //Define a template function g that computes a complex-valued electromagnetics Green’s function. The function takes distance, frequency, and material parameters, computes a complex propagation constant gamma, and returns the exponentially decaying field divided by 4piR.
  /**
   * @brief Computes the scalar Green's function for electromagnetic wave propagation. 
   *
   * Computes:
   * \f[
   * g(R) = \frac{\exp(-\gamma R)}{4 \pi R}
   * \f]
   * where
   * \f[
   * \gamma = i \, \omega \, \sqrt{\epsilon_r \mu_r}\sqrt{\epsilon_0 \mu_0}.
   * \f]
   *
   * @tparam Real Floating point type (default: double).
   * @tparam MatValueType Type of material parameters.
   * @param R Distance between source and observation points.
   * @param omega Angular frequency (rad/s).
   * @param epsr Relative permittivity.
   * @param mur Relative permeability.
   * @return Complex-valued Green's function.
   */
  template<class Real = double, class MatValueType> //Template function Real deafauls to double if not specified, MatValueType must be provided
  std::complex<Real> gScaled(
    const Real& R, 
    const Real& omega, 
    const MatValueType& epsr, 
    const MatValueType& mur, 
    const double wavelength)
  {
    auto tester = constants<Real>::ROOT_EPS0MU0_;
    const Real Rscaled = R*wavelength;
    const std::complex<Real> gamma = std::complex<Real>(0.0,1.0)*omega*sqrt(epsr*mur)*constants<Real>::ROOT_EPS0MU0_;
    const auto test = std::exp(-gamma*R);
    return std::exp(-gamma*Rscaled)/(Real(4.0)*constants<Real>::PI*Rscaled);
  }

  template<class Real = double, class MatValueType>
  std::complex<Real> g(
    const Real& R, 
    const Real& omega, 
    const MatValueType& epsr, 
    const MatValueType& mur)
  {
    auto tester = constants<Real>::ROOT_EPS0MU0_; // $$\sqrt{\epsilon_0\mu_0}$$

    const std::complex<Real> gamma = std::complex<Real>(0.0,1.0)*omega*sqrt(epsr*mur)*constants<Real>::ROOT_EPS0MU0_;//(0.0,1.0) is the imaginary unit i, omega*sqrt(epsr*mur) scaling factor, this whole thing is a propagation constant
    const auto test = std::exp(-gamma*R);
    return std::exp(-gamma*R)/(Real(4.0)*constants<Real>::PI*R);
  }

  template<class Real = double, class MatValueType>
  std::complex<Real> gPrint(
    const Real& R, 
    const Real& omega, 
    const MatValueType& epsr, 
    const MatValueType& mur)
  {
    auto tester = constants<Real>::ROOT_EPS0MU0_;
    const std::complex<Real> gamma = std::complex<Real>(0.0,1.0)*omega*sqrt(epsr*mur)*constants<Real>::ROOT_EPS0MU0_;
    const auto test = std::exp(-gamma*R);//$$\exp(-\gamma R)$$

    std::cout << "gamma: " << gamma << "  R: " << R << "   phase: " << std::exp(-gamma*R) << "   denom: " << Real(4.0)*constants<Real>::PI*R << std::endl;

    return std::exp(-gamma*R)/(Real(4.0)*constants<Real>::PI*R);//$$\frac{\exp(-\gamma R)}{4\pi R}$$
  }

  // template<class Real = double, class MatValueType>
  // std::complex<Real> dgdw(const Real& R, const Real& omega, const MatValueType& epsr, const MatValueType& mur)
  // {
  //   auto tester = constants<Real>::ROOT_EPS0MU0_;

  //   const std::complex<Real> gamma = std::complex<Real>(0.0,1.0)*omega*sqrt(epsr*mur)*constants<Real>::ROOT_EPS0MU0_;
  //   const auto test = std::exp(-gamma*R);
  //   return -gamma/omega*R;
  // }
  //This defines two functions grad_g. Each computes the gradient of the scalar Green’s function in the form of a complex-valued vector (Point). One overload computes the distance norm R internally, the other takes it as a parameter (to avoid recomputing if it’s known). Both functions build a complex decay factor and scale Rvec accordingly
  /**
   * @brief Computes the gradient of the scalar Green's function.
   *
   * The gradient is given by:
   * \f[
   * \nabla g(\mathbf{R}) = 
   * -\frac{\gamma^2}{4\pi R}
   * e^{-\gamma R}
   * \left(
   * \frac{1}{\gamma R} + \frac{1}{\gamma^2 R^2}
   * \right) \mathbf{R}.
   * \f]
   *
   * @tparam spacedim Number of spatial dimensions.
   * @tparam Real Floating point type (default: double).
   * @tparam MatValueType Type of material parameters.
   * @param Rvec Vector between source and observation points.
   * @param omega Angular frequency (rad/s).
   * @param epsr Relative permittivity.
   * @param mur Relative permeability.
   * @return Complex-valued gradient vector.
   */
  template <unsigned int spacedim, class Real = double, class MatValueType>
  Point<spacedim, std::complex<Real>> grad_gScaled(
    const Point<spacedim, Real>& Rvec, 
    const Real& omega, 
    const MatValueType& epsr, 
    const MatValueType& mur, 
    const double wavelength)
  {
    std::complex<Real> gamma = std::complex<Real>(0.,1.)*omega*sqrt(epsr*mur)*constants<Real>::ROOT_EPS0MU0_;
    const Real R = Rvec.norm(); //$$R=||R_{vec}||$$ euclidean distance
    const Real Rscaled = R*wavelength; //Attempt at WAvelength rescaling 
    return (-gamma*gamma/(4*constants<Real>::PI*Rscaled)*
    std::exp(-1.*gamma*Rscaled)*
    (1./(gamma*Rscaled)+1./(gamma*gamma*Rscaled*Rscaled)))*
    Rvec*wavelength;//This is the analytical gradient of the Green's function with respect to position
    // Prefactor: $$\frac{\gamma^2}{4\pi R}$$, Decay Factor: $$\exp(-\gamma R)$$, Scaling Term: $$(\frac{1}{\gamma R}+\frac{1}{\gamma^2 R^2})$$, Multiply by Rvec: yields the vector gradient
    //$$\frac{\gamma^2 e^{-\gamma R}}{4\pi R}(\frac{1}{\gamma R}+\frac{1}{\gamma^2 R^2})R_{vec}$$
  }

  // template <unsigned int spacedim, class Real = double, class MatValueType>
  // Point<spacedim, std::complex<Real>> grad_g(const Point<spacedim, Real>& Rvec, const Real& omega, const MatValueType& epsr, const MatValueType& mur)
  // {
  //   std::complex<Real> gamma = std::complex<Real>(0.,1.)*omega*sqrt(epsr*mur)*constants<Real>::ROOT_EPS0MU0_;
  //   const Real R = Rvec.norm();
  //   return (-gamma*gamma/(4*constants<Real>::PI*R)*std::exp(-1.*gamma*R)*(1./(gamma*R)+1./(gamma*gamma*R*R)))*Rvec;
  // }
  /**
   * @brief Computes the gradient of the scalar Green's function, with precomputed distance.
   *
   * This overload is identical to the other gradient version but takes the distance norm R explicitly.
   *
   * The gradient is given by:
   * \f[
   * \nabla g(\mathbf{R}) = 
   * -\frac{\gamma^2}{4\pi R}
   * e^{-\gamma R}
   * \left(
   * \frac{1}{\gamma R} + \frac{1}{\gamma^2 R^2}
   * \right) \mathbf{R}.
   * \f]
   *
   * @tparam spacedim Number of spatial dimensions.
   * @tparam Real Floating point type (default: double).
   * @tparam MatValueType Type of material parameters.
   * @param Rvec Vector between source and observation points.
   * @param R Norm of Rvec.
   * @param omega Angular frequency (rad/s).
   * @param epsr Relative permittivity.
   * @param mur Relative permeability.
   * @return Complex-valued gradient vector.
   */
  template <unsigned int spacedim, class Real = double, class MatValueType>//Same as above but takes a precalculated R instead of the norm of Rvec
  Point<spacedim, std::complex<Real>> grad_g(
    const Point<spacedim, Real>& Rvec, 
    const Real& R, 
    const Real& omega, 
    const MatValueType& epsr, 
    const MatValueType& mur)
  {
    std::complex<Real> gamma = std::complex<Real>(0.,1.)*omega*sqrt(epsr*mur)*constants<Real>::ROOT_EPS0MU0_;
    // std::cout << "checking: " <<(-gamma*gamma/(4*constants<Real>::PI*R)*std::exp(-1.*gamma*R)*(1./(gamma*R)+1./(gamma*gamma*R*R)))*Rvec << "  " << (-std::exp()/(4*constants<Real>::PI*R*R*R)*(gamma*R+1)*Rvec)
    return (-gamma*gamma/(4*constants<Real>::PI*R)*
    std::exp(-1.*gamma*R)*
    (1./(gamma*R)+1./(gamma*gamma*R*R)))*
    Rvec;//This is the analytical gradient of the Green's function with respect to position
    // Prefactor: $$\frac{\gamma^2}{4\pi R}$$, Decay Factor: $$\exp(-\gamma R)$$, Scaling Term: $$(\frac{1}{\gamma R}+\frac{1}{\gamma^2 R^2})$$, Multiply by Rvec: yields the vector gradient
    //$$\frac{\gamma^2 e^{-\gamma R}}{4\pi R}(\frac{1}{\gamma R}+\frac{1}{\gamma^2 R^2})R_{vec}$$
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
