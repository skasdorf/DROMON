//
// Created by Jake J. Harmon (jake.harmon@ieee.org) on 7/16/21.
//

#ifndef DROMON_CONFIG_H
#define DROMON_CONFIG_H

#include <assert.h>
#include <iostream>
#include <vector>
#include <utility>
#include <limits>
#include <complex>

#define DROMON_NAMESPACE_OPEN namespace  dromon {
#define DROMON_NAMESPACE_CLOSE }


#define LINEARP 0
#define QUADRATICP 1
#define CUBICP 2


DROMON_NAMESPACE_OPEN
enum DoFDirection
{
    u_dir, v_dir, w_dir
};

enum CurrentType
{
    Electric,
    Magnetic,
    Nothing
};
enum UpdateFlags
{
  update_default = 0,
  output_active = 0x000001,
  verbose_output = 0x000040,
  suppress_comments = 0x000002
};

enum RegularityType
{
  SelfSingular, // A cell overlap singularity
  FaceSingular, // A face singularity
  EdgeSingular, // An edge singularity
  PointSingular, // A point singularity
  FaceNearSingular,
  EdgeNearSingular,
  PointNearSingular,
  Regular // No singular behavior
};

enum AdjacencyType
{
  Disjoint,
  Self,
  Face,
  Edge,
  Vertex,
};

inline UpdateFlags
operator|(const UpdateFlags f1, const UpdateFlags f2)
{
  return static_cast<UpdateFlags>(static_cast<unsigned int>(f1) |
                                  static_cast<unsigned int>(f2));
}

inline UpdateFlags operator&(const UpdateFlags f1, const UpdateFlags f2)
{
  return static_cast<UpdateFlags>(static_cast<unsigned int>(f1) &
                                  static_cast<unsigned int>(f2));
}

template <unsigned int dim, unsigned int spacedim, unsigned int surf_type>
struct GeometryInfo
{
  static constexpr unsigned int nodes_per_cell = dim == 1 ? surf_type + 2 : dim == 2 ? (surf_type+2)*(surf_type+2) : (surf_type+2)*(surf_type+2)*(surf_type+2);
  static constexpr unsigned int vertices_per_cell = 1 << dim;
  static constexpr unsigned int faces_per_cell = 2*dim;
  static constexpr unsigned int vertices_per_face = 2 << (dim-1);
};

inline std::ostream& operator<<(std::ostream& os, const DoFDirection& dir)
{
  if (dir == u_dir)
    return os << "u_dir";
  else if (dir == v_dir)
    return os << "v_dir";
  else
   return os << "w_dir";

}

inline std::ostream& operator<<(std::ostream& os, const CurrentType& cur_type)
{
  if (cur_type == Electric)
    return os << "Electric";
  else if (cur_type == Magnetic)
    return os << "Magnetic";
  else
    return os << "Nothing";
}

//namespace constants
//{
//  static constexpr double PI = 3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117067982148086513282306647;
//  static constexpr double C0 = 299792458.;
//  static constexpr double EPS0 = 8.8541878128e-12;
//  static constexpr double MU0 = PI*4.0e-7;
//  static constexpr double ROOT_EPS0MU0_ = 3.3356409510735270059318962774541534256017524049995961e-9;
//
//}

template <class Real>
struct constants
{
  static constexpr Real PI = 3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117067982148086513282306647;
  static constexpr Real C0 = 299792458.;
  static constexpr Real EPS0 = 8.8541878128e-12;
  static constexpr Real MU0 = PI*4.0e-7;
  static constexpr Real ROOT_EPS0MU0_ = 3.3356409510735270059318962774541534256017524049995961e-9;
  static constexpr std::complex<Real> complexj = {Real(0), Real(1)};
};

//template <>
//struct constants<double>
//{
//    static constexpr double PI = 3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117067982148086513282306647;
//    static constexpr double C0 = 299792458.;
//    static constexpr double EPS0 = 8.8541878128e-12;
//    static constexpr double MU0 = PI*4.0e-7;
//    static constexpr double ROOT_EPS0MU0_ = 3.3356409510735270059318962774541534256017524049995961e-9;
//    static constexpr std::complex<double> complexj = {double(0), double(1)};
//};




DROMON_NAMESPACE_CLOSE

#endif //DROMON_CONFIG_H
