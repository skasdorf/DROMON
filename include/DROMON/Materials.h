//
// Created by Jake J. Harmon (jake.harmon@ieee.org) on 4/10/22.
//

#ifndef DROMON_MATERIALS_H
#define DROMON_MATERIALS_H
#include "config.h"
#include <complex>

DROMON_NAMESPACE_OPEN

template <class CoefficientType = double, class Real = double>
struct Material
{
  Material();
  Material(const CoefficientType& epsr, const CoefficientType& mur, const bool& is_PEC, const bool& is_PMC = false);
  const CoefficientType epsr;
  const CoefficientType mur;
  const CoefficientType eps;
  const CoefficientType mu;
  const bool is_PEC;
  const bool is_PMC;
};
template <class CoefficientType, class Real> Material<CoefficientType, Real>::Material() : epsr(1.0), mur(1.0), eps(constants<Real>::EPS0), mu(constants<Real>::MU0), is_PEC(false), is_PMC(false)
{}
template <class CoefficientType, class Real>
Material<CoefficientType, Real>::Material(const CoefficientType &epsr,
                         const CoefficientType &mur, const bool& is_PEC, const bool& is_PMC) : epsr(epsr), mur(mur), eps(epsr*constants<Real>::EPS0), mu(mur*constants<Real>::MU0), is_PEC(is_PEC), is_PMC(is_PMC) {}


template <class CoefficientType = double>
class MaterialDomain
{
public:
  MaterialDomain();
  MaterialDomain(const Material<CoefficientType>& interior, const Material<CoefficientType>& exterior);
  const Material<CoefficientType>& get_exterior() const;
  const Material<CoefficientType>& get_interior() const;
private:
  const Material<CoefficientType> interior;
  const Material<CoefficientType> exterior;
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

template <class CoefficientType = double>
class MaterialData
{
public:
  MaterialData() = default;
  void push_back(const MaterialDomain<CoefficientType>& domain);
  const MaterialDomain<CoefficientType>& get_material_domain(const unsigned int& i) const;
private:
  std::vector<MaterialDomain<CoefficientType>> domains;
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
