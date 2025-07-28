//
// Created by Jake J. Harmon (jake.harmon@ieee.org) on 4/10/22.
// Edited by Christopher A. Erickson (Christopher.Erickson@ieee.org) on 7/3/25
//

#ifndef DROMON_MATERIALS_H
#define DROMON_MATERIALS_H
#include "config.h"
#include <complex>

DROMON_NAMESPACE_OPEN

/**
 * @brief Represents the material properties in an electromagnetic simulation.
 * 
 * @tparam CoefficientType The type used for material coefficients (e.g., double, std::complex).
 * @tparam Real The real number type used for constants.
 */
template <class CoefficientType = double, class Real = double>
struct Material
{
  /**
   * @brief Default constructor. Initializes to vacuum material (epsr = 1, mur = 1).
   */
  Material();
  /**
   * @brief Constructs a material with specified relative permittivity and permeability.
   * 
   * @param epsr Relative permittivity.
   * @param mur Relative permeability.
   * @param is_PEC Flag indicating if the material is a Perfect Electric Conductor (PEC).
   * @param is_PMC Flag indicating if the material is a Perfect Magnetic Conductor (PMC).
   */
  Material(const CoefficientType& epsr, const CoefficientType& mur, const bool& is_PEC, const bool& is_PMC = false);
  /** Relative permittivity (dimensionless). */
  const CoefficientType epsr;
  /** Relative permeability (dimensionless). */
  const CoefficientType mur;
  /** Absolute permittivity (epsr * EPS0). */
  const CoefficientType eps;
  /** Absolute permeability (mur * MU0). */
  const CoefficientType mu;
  /** Indicates whether the material is a Perfect Electric Conductor. */
  const bool is_PEC;
  /** Indicates whether the material is a Perfect Magnetic Conductor. */
  const bool is_PMC;
};
template <class CoefficientType, class Real> Material<CoefficientType, Real>::Material() : epsr(1.0), mur(1.0), eps(constants<Real>::EPS0), mu(constants<Real>::MU0), is_PEC(false), is_PMC(false)
{}
template <class CoefficientType, class Real>
Material<CoefficientType, Real>::Material(const CoefficientType &epsr,
                         const CoefficientType &mur, const bool& is_PEC, const bool& is_PMC) : epsr(epsr), mur(mur), eps(epsr*constants<Real>::EPS0), mu(mur*constants<Real>::MU0), is_PEC(is_PEC), is_PMC(is_PMC) {}

/**
 * @brief Represents a domain consisting of an interior and exterior material.
 * 
 * @tparam CoefficientType The type used for material coefficients.
 */
template <class CoefficientType = double>
class MaterialDomain
{
public:
  /**
   * @brief Default constructor. Initializes both interior and exterior to default materials.
   */
  MaterialDomain();
  /**
   * @brief Constructs a material domain with specified interior and exterior materials.
   * 
   * @param interior The material inside the domain.
   * @param exterior The material outside the domain.
   */
  MaterialDomain(const Material<CoefficientType>& interior, const Material<CoefficientType>& exterior);
  /**
   * @brief Gets the exterior material.
   * 
   * @return Reference to the exterior material.
   */
  const Material<CoefficientType>& get_exterior() const;
  /**
   * @brief Gets the interior material.
   * 
   * @return Reference to the interior material.
   */
  const Material<CoefficientType>& get_interior() const;
private:
  const Material<CoefficientType> interior; /**< Interior material. */
  const Material<CoefficientType> exterior; /**< Exterior material. */
};
template <class CoefficientType> MaterialDomain<CoefficientType>::MaterialDomain() : interior(Material<CoefficientType>()), exterior(Material<CoefficientType>())
{}
template <class CoefficientType>
MaterialDomain<CoefficientType>::MaterialDomain(const Material<CoefficientType> &interior,
                                     const Material<CoefficientType> &exterior) : interior(interior), exterior(exterior)
{}
template <class CoefficientType>
const Material<CoefficientType> &MaterialDomain<CoefficientType>::get_exterior() const {
  return exterior;
}
template <class CoefficientType>
const Material<CoefficientType> &MaterialDomain<CoefficientType>::get_interior() const {
  return interior;
}
/**
 * @brief Container for storing multiple material domains.
 * 
 * @tparam CoefficientType The type used for material coefficients.
 */
template <class CoefficientType = double>
class MaterialData
{
public:
  /**
   * @brief Default constructor.
   */
  MaterialData() = default;
  /**
   * @brief Adds a material domain to the container.
   * 
   * @param domain The material domain to add.
   */
  void push_back(const MaterialDomain<CoefficientType>& domain);
  /**
   * @brief Retrieves a material domain by index.
   * 
   * @param i Index of the material domain.
   * @return Reference to the requested material domain.
   */
  const MaterialDomain<CoefficientType>& get_material_domain(const unsigned int& i) const;
private:
  std::vector<MaterialDomain<CoefficientType>> domains; /**< Stored material domains. */
};

template <class CoefficientType>
void MaterialData<CoefficientType>::push_back(const MaterialDomain<CoefficientType> &domain) {
  domains.push_back(domain);
}
template <class CoefficientType>
const MaterialDomain<CoefficientType> &MaterialData<CoefficientType>::get_material_domain(const unsigned int &i) const {
  return domains[i];
}

DROMON_NAMESPACE_CLOSE

#endif // DROMON_MATERIALS_H
