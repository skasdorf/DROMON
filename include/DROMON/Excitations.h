//
// Created by Jake J. Harmon (jake.harmon@ieee.org) on 4/10/22.
//

#ifndef DROMON_EXCITATIONS_H
#define DROMON_EXCITATIONS_H
#include "config.h"

DROMON_NAMESPACE_OPEN
namespace Excitations
{
template<unsigned int spacedim, class Real = double>
struct Excitation
{
  explicit Excitation(const Real& frequency);
  // const Real freq;
  std::complex<Real> gamma;
  const Real omega;
  const Real frequency;


  // virtual void set_incident_direction(const Real& theta_inc, const Real& phi_inc);
  virtual Point<spacedim, std::complex<Real>> evaluate_excitation(const Point<spacedim, Real>& R) const;
  virtual std::complex<Real> evaluate_excitation_in_direction(const Point<spacedim, Real>& R, const Point<spacedim, Real>& direction) const;
  virtual std::complex<Real> evaluate_excitation_in_direction_mat(const Point<spacedim, Real>& R, const Point<spacedim, Real>& direction, double epsr) const;
  virtual std::complex<Real> evaluate_excitation_in_direction_perturb(const Point<spacedim, Real>& R, const Point<spacedim, Real>& direction, int perturbVar, double perturbValue) const;
  virtual std::complex<Real> evaluate_hops_excitation_in_direction(const Point<spacedim, Real>& R, const Point<spacedim, Real>& direction) const;
  virtual Real magnitude() const;
};
template < unsigned int spacedim, class Real>
Excitation<spacedim, Real>::Excitation(const Real& frequency) : frequency(frequency), omega(2.0*constants<Real>::PI*frequency), 
gamma(2.0*constants<Real>::PI*frequency*constants<Real>::ROOT_EPS0MU0_*std::complex<Real>(0.0,1.0)) {}

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

// E_inc(R) = E_inc_mag*exp(-j*beta*r\cdot n_inc_hat
// where beta = \omega*sqrt(eps0*mu0)
template<unsigned int spacedim, class Real = double>
struct PlaneWave : public Excitation<spacedim, Real>
{
  explicit PlaneWave(const Real& frequency, const Point<spacedim, std::complex<Real>>& E_mag);

  virtual Point<spacedim, std::complex<Real>> evaluate_excitation(const Point<spacedim, Real>& R) const override;
  virtual std::complex<Real> evaluate_excitation_in_direction(const Point<spacedim, Real>& R, const Point<spacedim, Real>& direction) const override;
  virtual std::complex<Real> evaluate_excitation_in_direction_mat(const Point<spacedim, Real>& R, const Point<spacedim, Real>& direction, double epsr) const override;
  virtual std::complex<Real> evaluate_excitation_in_direction_perturb(const Point<spacedim, Real>& R, const Point<spacedim, Real>& direction, int perturbVar, double perturbValue) const override;
  virtual std::complex<Real> evaluate_hops_excitation_in_direction(const Point<spacedim, Real>& R, const Point<spacedim, Real>& direction) const;
  virtual Real magnitude() const;
};

template <class Real>
struct PlaneWave<3, Real> : public Excitation<3, Real> {
  explicit PlaneWave(const Real& frequency, const Point<2, std::complex<Real>>& E_mag, const Real& theta_inc, const Real& phi_inc, const MaterialData<Real> &materialData);
  virtual Point<3, std::complex<Real>> evaluate_excitation(const Point<3, Real>& R) const override;
  virtual std::complex<Real> evaluate_excitation_in_direction(const Point<3, Real>& R, const Point<3, Real>& direction) const override;
  virtual std::complex<Real> evaluate_excitation_in_direction_mat(const Point<3, Real>& R, const Point<3, Real>& direction, double epsr) const override;
  virtual std::complex<Real> evaluate_excitation_in_direction_perturb(const Point<3, Real>& R, const Point<3, Real>& direction, int perturbVar, double perturbValue) const override;
  virtual std::complex<Real> evaluate_hops_excitation_in_direction(const Point<3, Real>& R, const Point<3, Real>& direction) const;
  virtual Real magnitude() const override;
  void set_incident_direction(const Real& theta_inc, const Real& phi_inc);

  Point<3, Real> n_hat;
  Point<3, Real> theta_hat;
  Point<3, Real> phi_hat;

private:
  Point<3, std::complex<Real>> E_mag_cart;
};
template <class Real>
PlaneWave<3, Real>::PlaneWave(const Real& frequency, const Point<2, std::complex<Real>>& E_mag, const Real& theta_inc, const Real& phi_inc, const MaterialData<Real> &materialData) : Excitation<3, Real>(frequency)
{
  //// this needs to be double checked
  const Material<Real> &exterior_material =
      (materialData.get_material_domain(0))
          .get_exterior();
  double epsr = exterior_material.epsr;
  this->gamma = this->gamma * sqrt(epsr);
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
    // std::cout << "gamma: " << this->gamma << std::endl;

  const std::complex<Real> g = exp(-this->gamma*R_dot_n_hat);
  return g*(E_mag_cart(0)*direction(0) + E_mag_cart(1)*direction(1) + E_mag_cart(2)*direction(2));
}

template <class Real>
std::complex<Real> PlaneWave<3, Real>::evaluate_excitation_in_direction_mat(
    const Point<3, Real> &R, const Point<3, Real> &direction, double epsr) const
{
  const Real R_dot_n_hat = R.dot(n_hat);
    // std::cout << "gamma: " << this->gamma << std::endl;

  const std::complex<Real> g = exp(-this->gamma*R_dot_n_hat*sqrt(epsr));
  return g*(E_mag_cart(0)*direction(0) + E_mag_cart(1)*direction(1) + E_mag_cart(2)*direction(2));
}

template <class Real>
std::complex<Real> PlaneWave<3, Real>::evaluate_excitation_in_direction_perturb(
    const Point<3, Real> &R, const Point<3, Real> &direction, int perturbVar, double perturbValue) const
{

  const Real R_dot_n_hat = R.dot(n_hat);
  std::complex<double> g;
  
  //g==exp(-j*omega*sqrt(eps*mu)*R_dot_n) which is why material is sqrt(perturbValue)
  //perturb material
  if (perturbVar == 1){
    g = exp(-this->gamma*R_dot_n_hat * sqrt(perturbValue));
  }
  //perturb freq
  else if (perturbVar == 2){
    g = exp(-this->gamma*R_dot_n_hat * perturbValue);
  }
  //perturb radius
  else if (perturbVar == 3){
    g = exp(-this->gamma*R_dot_n_hat * perturbValue);
  }

  return g*(E_mag_cart(0)*direction(0) + E_mag_cart(1)*direction(1) + E_mag_cart(2)*direction(2));
}


template <class Real> Real PlaneWave<3, Real>::magnitude() const {
  return E_mag_cart.norm();
}

template <class Real>
std::complex<Real> PlaneWave<3, Real>::evaluate_hops_excitation_in_direction(
    const Point<3, Real> &R, const Point<3, Real> &direction) const
{
  std::complex<Real> complexj(0.0, 1.0);
  const Real R_dot_n_hat = R.dot(n_hat);
  const std::complex<Real> g = exp(-this->gamma*R_dot_n_hat);
  return g*(-complexj*R_dot_n_hat*sqrt(8.8541878128e-12 * 12.566370614359e-7))*(E_mag_cart(0)*direction(0) + E_mag_cart(1)*direction(1) + E_mag_cart(2)*direction(2));
}


}

DROMON_NAMESPACE_CLOSE

#endif // DROMON_EXCITATIONS_H
