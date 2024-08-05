//
// Created by Jake J. Harmon (jake.harmon@ieee.org) on 3/31/22.
//
#ifndef DROMON_DOFMASK_H
#define DROMON_DOFMASK_H

#include "config.h"

DROMON_NAMESPACE_OPEN
// DoFMask is essentialy a vector of booleans with some additional
// functionality to allow "turning off" DoFs when integrating contributions
// for assembling Galerkin systems from previous iterations.
class DoFMask
{
public:
  DoFMask() = default;
  explicit DoFMask(const unsigned int& n);
  void mask_on();
  void mask_off();
  void mask_on(const unsigned int& i);
  void mask_off(const unsigned int& i);

  const unsigned int& size() const;

  void resize(const unsigned int& new_n);

  const char& operator()(const unsigned int& i) const { return mask[i]; }
  const char& operator[](const unsigned int& i) const { return mask[i]; }

private:
  unsigned int length;
  std::vector<char> mask;
};
DoFMask::DoFMask(const unsigned int &n)  {
  mask = std::vector<char>(n, true);
  length = n;
}
void DoFMask::mask_on() {
  mask = std::vector<char>(mask.size(), false);
}
void DoFMask::mask_off() {
  mask = std::vector<char>(mask.size(), true);
}
void DoFMask::mask_on(const unsigned int &i) {
  assert(i < mask.size() && "Index i out of range!");
  mask[i] = false;
}
void DoFMask::mask_off(const unsigned int &i) {
  assert(i < mask.size() && "Index i out of range!");
  mask[i] =  true;
}
void DoFMask::resize(const unsigned int &new_n)
{
  mask.resize(new_n, true);
  length = new_n;
}

const unsigned int& DoFMask::size() const
{
  return length;
}

//// DoFMask is essentialy a vector of booleans with some additional
//// functionality to allow "turning off" DoFs when integrating contributions
//// for assembling Galerkin systems from previous iterations.
//class DoFMask
//{
//public:
//  DoFMask() = default;
//  explicit DoFMask(const unsigned int& n);
//  void mask_on();
//  void mask_off();
//  void mask_on(const unsigned int& i);
//  void mask_off(const unsigned int& i);
//
//  const unsigned int& size() const;
//
//  void resize(const unsigned int& new_n);
//
//  const auto& operator()(const unsigned int& i) const { return mask[i]; }
//  const auto& operator[](const unsigned int& i) const { return mask[i]; }
//
//private:
//  unsigned int length;
//  //std::vector<bool> mask;
//  std::unique_ptr<bool[]> mask;
//};
//DoFMask::DoFMask(const unsigned int &n)  {
//  //mask = std::vector<bool>(n, true);
//  mask = std::unique_ptr<bool[]>(new bool[n]);
//  length = n;
//  std::fill_n(mask.get(), length, true);
//}
//void DoFMask::mask_on() {
//  std::fill_n(mask.get(), length, false);
//}
//void DoFMask::mask_off() {
//  std::fill_n(mask.get(), length, true);
//}
//void DoFMask::mask_on(const unsigned int &i) {
//  assert(i < length && "Index i out of range!");
//  (mask)[i] = false;
//}
//void DoFMask::mask_off(const unsigned int &i) {
//  assert(i < length && "Index i out of range!");
//  mask[i] =  true;
//}
//void DoFMask::resize(const unsigned int &new_n)
//{
//  mask.resize(new_n, true);
//  length = new_n;
//}
//
//const unsigned int& DoFMask::size() const
//{
//  return length;
//}

DROMON_NAMESPACE_CLOSE
#endif // DROMON_DOFMASK_H
