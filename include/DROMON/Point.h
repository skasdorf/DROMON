//
// Created by Jake J. Harmon (jake.harmon@ieee.org) on 7/18/21.
//

#ifndef DROMON_POINT_H
#define DROMON_POINT_H

#include "config.h"
#include "vector"
#include <cmath>
#include <complex>

//template <typename T> struct add_ref_const { typedef T const type; };
//template <typename T> struct add_ref_const<T&> { typedef T const& type; };

template <typename T>
struct norm_type
{
  typedef T type;
};

template <>
struct norm_type<std::complex<double>>
{
  typedef double type;
};

template <>
struct norm_type<std::complex<float>>
{
  typedef float type;
};


template<unsigned int dim, class Real>
    class Point {
    public:
        Point()  {
            for (auto &coord : coords)
                coord = (Real) 0;
        }

        Point(std::vector<Real> coords_in);
        Point(Real x);
        Point (Real x, Real y);
        Point (Real x, Real y, Real z);

        const Real& operator()(const unsigned int& pos) const;
        Real& operator[](const unsigned int& pos);

        Real dot(const Point<dim, Real>& p) const;
        std::complex<Real> dot(const Point<dim, std::complex<Real>>& p) const;

        Point<dim, Real> operator+(const Point<dim, Real>& p) const;
        Point<dim, Real> operator-(const Point<dim, Real>& p) const;
        Point<dim, Real> operator*(const Real& scalar) const;
        void multiply_scalar(const Real& scalar, Point<dim, Real>* output) const;
        void add_multiple_of(const Point<dim, Real>& p, const Real& scalar);
        // Point<dim, Real> cross(Point<dim, Real>& p) const;
        void operator/=(const Real& scalar);
        void operator*=(const Real& scalar);
        void operator+=(const Point<dim, Real>& p);
        void operator^=(const Real& scalar);
       // Real norm() const;
       auto norm() const -> typename norm_type<Real>::type;
       Point<dim, Real> conj() const;
       friend bool operator==(const Point<dim, Real> &p1,
                              const Point<dim, Real> &p2)
       {
         for (unsigned int i =0; i < dim; ++i)
           if (p1.coords[i] != p2.coords[i])
             return false;
         return true;
       }

       friend bool operator!=(const Point<dim, Real> &p1,
                              const Point<dim, Real> &p2)
       {
         return !(p1 == p2);
       }
    private:
        Real coords [dim];
    };

    template<unsigned int dim, class Real>
    Point<dim, Real>::Point(std::vector<Real> coords_in)  {
#ifdef DEBUG
        assert(coords_in.size() == dim && "Mismatch in the dimension of the Point and the input coordinates!");
#endif
        for (unsigned int i = 0; i < dim; ++i)
        coords[i]=coords_in[i];
    }

    template<unsigned int dim, class Real>
    const Real& Point<dim, Real>::operator()(const unsigned int& pos) const {
#ifdef DEBUG
        assert(pos <  dim && "Input position out of range for coords!");
#endif
        return this->coords[pos];
    }

template<unsigned int dim, class Real>
Point<dim, Real> Point<dim, Real>::operator+(const Point<dim, Real>& p) const {
    Point<dim, Real> out;
    for (unsigned int i = 0; i < dim; ++i)
        out[i] = this->coords[i] + p(i);
    return out;
}

template<unsigned int dim, class Real>
Real &Point<dim, Real>::operator[](const unsigned int& pos) {
    return coords[pos];
}

template<unsigned int dim, class Real>
Point<dim, Real> Point<dim, Real>::operator-(const Point<dim, Real>& p) const {
    Point<dim, Real> out;
    for (unsigned int i = 0; i < dim; ++i)
        out[i] = this->coords[i] - p(i);
    return out;
}

template<unsigned int dim, class Real>
void Point<dim, Real>::operator/=(const Real& scalar) {
    for (auto& coord : coords)
        coord /= scalar;
}

template<unsigned int dim, class Real>
void Point<dim, Real>::operator+=(const Point<dim, Real> &p) {
    for (unsigned int i = 0; i < dim; ++i)
        this->coords[i] += p(i);
}

template<unsigned int dim, class Real>
void Point<dim, Real>::operator*=(const Real& scalar) {
    for (auto& coord : coords)
        coord *= scalar;
}

template<unsigned int dim, class Real>
void Point<dim, Real>::operator^=(const Real& scalar) {
    for (auto& coord : coords)
        coord = std::pow(coord, scalar);
}

template<unsigned int dim, class Real>
Point<dim, Real> Point<dim, Real>::operator*(const Real& scalar) const {
    Point<dim, Real> out = *this;
    out *= scalar;
    return out;
}
template <unsigned int dim, class Real> Point<dim, Real>::Point(Real x)  {
#ifdef DEBUG
  assert(dim == 1 && "This constructor only works for dim == 1!");
#endif
  coords[0] = x;
}

template <unsigned int dim, class Real> Point<dim, Real>::Point(Real x, Real y)  {
#ifdef DEBUG
  assert(dim == 2 && "This constructor only works for dim == 2!");
#endif
  coords[0] = x;
  coords[1] = y;
}

template <unsigned int dim, class Real> Point<dim, Real>::Point(Real x, Real y, Real z)  {
  #ifdef DEBUG
  assert(dim == 3 && "This constructor only works for dim == 3!");
  #endif

  coords[0] = x;
  coords[1] = y;
  coords[2] = z;
}
template <unsigned int dim, class Real>
Real Point<dim, Real>::dot(const Point<dim, Real> &p) const {
  Real out = Real(0);
  for (unsigned int i = 0; i < dim; ++i)
    out += this->coords[i]*p(i);

  return out;
}
template <unsigned int dim, class Real>
std::complex<Real> Point<dim, Real>::dot(const Point<dim, std::complex<Real>> &p) const {
  std::complex<Real> out = std::complex<Real>(0);
  for (unsigned int i = 0; i < dim; ++i)
    out += this->coords[i]*p(i);

  return out;
}

template <unsigned int dim, class Real>
void Point<dim, Real>::add_multiple_of(const Point<dim, Real> &p,
                                       const Real &scalar)
{
  for (unsigned int i = 0; i < dim; ++i)
  {
    this->coords[i] += p(i)*scalar;
  }
}

template <>
inline double Point<3, double>::dot(const Point<3, double> &p) const{
  return this->coords[0]*p.coords[0]+this->coords[1]*p.coords[1]+this->coords[2]*p.coords[2];
}

template <>
inline double Point<2, double>::dot(const Point<2, double> &p) const{
  return this->coords[0]*p.coords[0]+this->coords[1]*p.coords[1];
}

template <>
inline Point<2,double> Point<2, double>::operator*(const double& scalar) const{
  return {this->coords[0]*scalar, this->coords[1]*scalar};
}

template <>
inline void Point<2, double>::operator*=(const double& scalar){
  this->coords[0]*=scalar;
  this->coords[1]*=scalar;
}
//template <>
//auto Point<3, double>::norm() -> decltype(abs(double(0))) const {
//  double sum = 0.0;
//  for (const auto& coord : coords) {
//    sum += coord * coord;
//  }
//  return std::sqrt(sum);
//}
//template <class Real>
//Real Point<2, Real>::dot(const Point<2, Real> &p) const
//{
//  return coords[0];
//}
//template <unsigned int dim, class Real>
//Point<dim, Real> Point<dim, Real>::cross(Point<dim, Real> &p) const {
//  Point<dim, Real> out;
//  for (unsigned int i = 0; i < dim; ++i)
//  {
//
//  }
template <class Real>
inline Point<3, Real> cross(const Point<3, Real>& p1, const Point<3, Real>& p2)
{
  Point<3, Real> out;
  out[0] = p1(1)*p2(2) - p1(2)*p2(1);
  out[1] = p1(2)*p2(0) - p1(0)*p2(2);
  out[2] = p1(0)*p2(1) - p1(1)*p2(0);
  return out;
}

template <class Real>
inline Point<3, Real> cross(const Point<2, Real>& p1, const Point<2, Real>& p2)
{
  Point<3, Real> out;
  out[2] = p1(0)*p2(1) - p1(1)*p2(0);
  return out;
}

template <class Real>
inline Point<3, Real> cross(const Point<1, Real>& p1, const Point<1, Real>& p2)
{
  Point<3, Real> out;
  return out;
}

template <class Real>
inline Point<3, std::complex<Real>> operator*(const std::complex<Real>& scalar, const Point<3, Real>& p)
{
  return {scalar*p(0), scalar*p(1), scalar*p(2)};
}

template <class Real>
inline Point<3, std::complex<Real>> operator*(const Real& scalar, const Point<3, std::complex<Real>>& p)
{
  return {scalar*p(0), scalar*p(1), scalar*p(2)};
}


//template<unsigned int dim, class Real>
//Real Point<dim, Real>::norm() const {
//  Real sum = (Real)0;
//  for (const auto& coord : coords)
//    sum += coord*coord;
//  return std::sqrt(sum);
//}

template <unsigned int dim, class Real>
auto Point<dim, Real>::norm() const -> typename norm_type<Real>::type
{
    typename norm_type<Real>::type sum(0);
    for (const auto& coord : coords) {
      const auto abs_val = std::abs(coord);
      sum += abs_val*abs_val;
    }
    return std::sqrt(sum);
}
template <unsigned int dim, class Real>
Point<dim, Real> Point<dim, Real>::conj() const {
  Point<dim, Real> out;
  for (unsigned int i = 0; i < dim; ++i)
    out[i] = std::conj(coords[i]);
  return out;
}


#endif //DROMON_POINT_H
