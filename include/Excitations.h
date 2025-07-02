//
// Created by Jake J. Harmon (jake.harmon@ieee.org) on 4/10/22.
// Edited by Christopher A. Erickson (Christopher.Erickson@ieee.org) on 7/2/25
//

#ifndef DROMON_EXCITATIONS_H
#define DROMON_EXCITATIONS_H
#include "config.h"

DROMON_NAMESPACE_OPEN
namespace Excitations
{
/**
 * @brief Abstract base class for excitations (e.g., incident fields) in simulation.
 *
 * @tparam spacedim The embedding dimension (1, 2, or 3).
 * @tparam Real Floating point type (default: double).
 *
 * This class stores frequency information and defines virtual methods to evaluate
 * the excitation field at a point.
 */
template<unsigned int spacedim, class Real = double>
struct Excitation
{
  /**
   * @brief Construct an excitation with a specified frequency.
   *
   * @param frequency Frequency in Hz.
   */
  explicit Excitation(const Real& frequency);
  // const Real freq;
  const std::complex<Real> gamma;
  const Real omega;
  const Real frequency;


  // virtual void set_incident_direction(const Real& theta_inc, const Real& phi_inc);
  /**
   * @brief Evaluate the excitation field at point R.
   *
   * @param R The spatial position.
   * @return Complex vector field at R.
   */
  virtual Point<spacedim, std::complex<Real>> evaluate_excitation(const Point<spacedim, Real>& R) const;
  /**
   * @brief Evaluate the excitation projected in a given direction.
   *
   * @param R The spatial position.
   * @param direction The direction vector.
   * @return Complex scalar projection.
   */
  virtual std::complex<Real> evaluate_excitation_in_direction(const Point<spacedim, Real>& R, const Point<spacedim, Real>& direction) const;
  /**
   * @brief Return the magnitude (norm) of the excitation field.
   *
   * @return Scalar magnitude.
   */
  virtual Real magnitude() const;
};
template < unsigned int spacedim, class Real>
Excitation<spacedim, Real>::Excitation(const Real& frequency) : frequency(frequency), omega(2*constants<Real>::PI*frequency), gamma(2*constants<Real>::PI*frequency*constants<Real>::ROOT_EPS0MU0_*std::complex<Real>(0,1.0)) {}

template <unsigned int spacedim, class Real>
Point<spacedim, std::complex<Real>>
Excitation<spacedim, Real>::evaluate_excitation(
    const Point<spacedim, Real> &R) const {
  assert(false && "Not implemented for base class!");
  return Point<spacedim, std::complex<Real>>();
}
template <unsigned int spacedim, class Real>
std::complex<Real> Excitation<spacedim, Real>::evaluate_excitation_in_direction(
    const Point<spacedim, Real> &R,
    const Point<spacedim, Real> &direction) const {
  assert(false && "Not implemented for base class!");
  return std::complex<Real>();
}
template <unsigned int spacedim, class Real>
Real Excitation<spacedim, Real>::magnitude() const {
  return Real(0);
}

// $$E_{inc}(R) = E_{inc_{mag}}*exp(-j*beta*r\cdot n_{inc}_{hat}$$
// where $$\beta = \omega*sqrt(eps0*mu0)$$
/**
 * @brief Represents a plane wave excitation.
 *
 * @tparam spacedim Embedding dimension.
 * @tparam Real Floating point type.
 */
template<unsigned int spacedim, class Real = double>
struct PlaneWave : public Excitation<spacedim, Real>
{
  /**
   * @brief Construct a plane wave excitation.
   *
   * @param frequency Frequency in Hz.
   * @param E_mag Field amplitude vector.
   */
  explicit PlaneWave(const Real& frequency, const Point<spacedim, std::complex<Real>>& E_mag);
  /**
   * @brief Evaluate the plane wave field at point R.
   *
   * @param R The spatial position.
   * @return Complex vector field.
   */
  virtual Point<spacedim, std::complex<Real>> evaluate_excitation(const Point<spacedim, Real>& R) const override;
  /**
   * @brief Evaluate the plane wave projected in a given direction.
   *
   * @param R The spatial position.
   * @param direction The direction vector.
   * @return Complex scalar projection.
   */
  virtual std::complex<Real> evaluate_excitation_in_direction(const Point<spacedim, Real>& R, const Point<spacedim, Real>& direction) const override;
  /**
   * @brief Return the magnitude of the field.
   *
   * @return Scalar magnitude.
   */
  virtual Real magnitude() const;
};
/**
 * @brief 3D specialization of PlaneWave with incidence angle specification.
 *
 * @tparam Real Floating point type.
 */
template <class Real>
struct PlaneWave<3, Real> : public Excitation<3, Real> {
  /**
   * @brief Construct a 3D plane wave excitation with specified incidence angles.
   *
   * @param frequency Frequency in Hz.
   * @param E_mag Complex vector of field amplitudes in (theta, phi) components.
   * @param theta_inc Incident elevation angle in radians.
   * @param phi_inc Incident azimuth angle in radians.
   */
  explicit PlaneWave(const Real& frequency, const Point<2, std::complex<Real>>& E_mag, const Real& theta_inc, const Real& phi_inc);
  /**
   * @brief Evaluate the 3D plane wave field at point R.
   *
   * @param R The spatial position.
   * @return Complex vector field.
   */
  virtual Point<3, std::complex<Real>> evaluate_excitation(const Point<3, Real>& R) const override;
  /**
   * @brief Evaluate the plane wave projected in a direction.
   *
   * @param R The spatial position.
   * @param direction The projection direction.
   * @return Complex scalar projection.
   */
  virtual std::complex<Real> evaluate_excitation_in_direction(const Point<3, Real>& R, const Point<3, Real>& direction) const override;
  /**
   * @brief Return the magnitude of the field.
   *
   * @return Scalar magnitude.
   */
  virtual Real magnitude() const override;
  /**
   * @brief Set the incident direction using theta and phi angles.
   *
   * @param theta_inc Elevation angle in radians.
   * @param phi_inc Azimuth angle in radians.
   */
  void set_incident_direction(const Real& theta_inc, const Real& phi_inc);

  Point<3, Real> n_hat;
  Point<3, Real> theta_hat;
  Point<3, Real> phi_hat;

private:
  Point<3, std::complex<Real>> E_mag_cart;
};
template <class Real>
PlaneWave<3, Real>::PlaneWave(const Real& frequency, const Point<2, std::complex<Real>>& E_mag, const Real& theta_inc, const Real& phi_inc) : Excitation<3, Real>(frequency)
{
  this->set_incident_direction(theta_inc, phi_inc);
  E_mag_cart = E_mag(0)*theta_hat + E_mag(1)*phi_hat;
}

template <class Real>
void PlaneWave<3, Real>::set_incident_direction(const Real &theta_inc,
                                                       const Real &phi_inc)
{
  const Real sin_theta = sin(theta_inc);
  const Real cos_theta = cos(theta_inc);
  const Real sin_phi = sin(phi_inc);
  const Real cos_phi = cos(phi_inc);

  n_hat = {-sin_theta * cos_phi, -sin_theta * sin_phi, -cos_theta};
  theta_hat = {cos_theta * cos_phi, cos_theta * sin_phi, -sin_theta};
  phi_hat = {-sin_phi, cos_phi, Real(0.0)};
}

template <class Real>
Point<3, std::complex<Real>>
PlaneWave<3, Real>::evaluate_excitation(
    const Point<3, Real> &R) const
{
  const Real R_dot_n_hat = R.dot(n_hat);
  const std::complex<Real> g = exp(-this->gamma*R_dot_n_hat);
  return E_mag_cart*g;
}
template <class Real>
std::complex<Real> PlaneWave<3, Real>::evaluate_excitation_in_direction(
    const Point<3, Real> &R, const Point<3, Real> &direction) const
{
  const Real R_dot_n_hat = R.dot(n_hat);
  const std::complex<Real> g = exp(-this->gamma*R_dot_n_hat);
  return g*(E_mag_cart(0)*direction(0) + E_mag_cart(1)*direction(1) + E_mag_cart(2)*direction(2));
}
template <class Real> Real PlaneWave<3, Real>::magnitude() const {
  return E_mag_cart.norm();
}

}

DROMON_NAMESPACE_CLOSE

#endif // DROMON_EXCITATIONS_H
