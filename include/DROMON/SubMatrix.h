//
// Created by Jake J. Harmon (jake.harmon@ieee.org) on 3/30/22.
// Edited by Christopher A. Erickson (Christopher.Erickson@ieee.org) on 7/3/25
//
/**
 * @file
 * @brief Data structures for dense and sparse matrix and vector assembly with numerical stability.
 *
 * This header defines classes for:
 * - Contiguous matrices for efficient storage
 * - Kahan summation vectors and matrices for compensated summation
 * - High-order compensated matrices
 * - Sparse-dense hybrid submatrices and subvectors
 *
 * These are used in assembling system matrices and excitation vectors in boundary element formulations.
 *
 */
/**
 * @defgroup AssemblyContainers
 * @brief Matrix and vector data structures for numerical assembly.
 */

#ifndef DROMON_SUBMATRIX_H
#define DROMON_SUBMATRIX_H

DROMON_NAMESPACE_OPEN
/**
 * @brief A contiguous 2D matrix with basic accumulation functionality.
 *
 * Stores a flat m x n array of coefficients with simple access and update methods.
 *
 * @tparam CoefficientType Type of the matrix entries.
 *
 * @ingroup AssemblyContainers
 */
template <class CoefficientType> class ContiguousMatrix {
public:
  ContiguousMatrix(const unsigned int &m, const unsigned int &n);

  CoefficientType &operator()(const unsigned int &i, const unsigned int &j);
  void accumulate(const unsigned int &i, const unsigned int &j,
                  const CoefficientType &value_to_add);
  void apply_carry() { return; };

private:
  const unsigned int m;
  const unsigned int n;
  std::vector<CoefficientType> matrix_data;
};
template <class CoefficientType>
ContiguousMatrix<CoefficientType>::ContiguousMatrix(const unsigned int &m,
                                                    const unsigned int &n)
    : m(m), n(n) {
  matrix_data = std::vector<CoefficientType>(m * n, CoefficientType(0));
}
template <class CoefficientType>
CoefficientType &
ContiguousMatrix<CoefficientType>::operator()(const unsigned int &i,
                                              const unsigned int &j) {
  return matrix_data[this->n * i + j];
}
template <class CoefficientType>
void ContiguousMatrix<CoefficientType>::accumulate(
    const unsigned int &i, const unsigned int &j,
    const CoefficientType &value_to_add) {
  this->operator()(i, j) += value_to_add;
}

/**
 * @brief A vector with Neumaier-compensated summation for improved numerical stability.
 *
 * Accumulation operations track correction terms to reduce floating-point error.
 *
 * @tparam CoefficientType Type of the vector entries.
 *
 * @ingroup AssemblyContainers
 */
template <class CoefficientType> class KahanVector
{
public:
  explicit KahanVector(const unsigned int &m);
  // CoefficientType &operator()(const unsigned int& i);
  CoefficientType &operator[](const unsigned int& i);
  void accumulate(const unsigned int& i, const CoefficientType& value_to_add);
  unsigned int size();
  void apply_carry();
  std::vector<CoefficientType>& data();
  /**
   * @brief Compute the total sum of all entries (excluding any un-applied carry).
   *
   * @return The sum of the vector entries.
   */
  CoefficientType sum() const;

private:
  const unsigned int m;
  std::vector<CoefficientType> vector;
  std::vector<CoefficientType> carry;
};
template <class CoefficientType>
KahanVector<CoefficientType>::KahanVector(const unsigned int &m) : m(m){
  vector = std::vector<CoefficientType>(m, CoefficientType(0));
  carry = std::vector<CoefficientType>(m, CoefficientType(0));
}

template <class CoefficientType>
CoefficientType &
KahanVector<CoefficientType>::operator[](const unsigned int &i) {
  return vector[i];
}
template <class CoefficientType>
void KahanVector<CoefficientType>::accumulate(const unsigned int &i,
                                              const CoefficientType &value_to_add)
{
  auto &current_value = vector[i];
  // Neumaier Version
  const auto t = current_value + value_to_add;
  if (std::abs(current_value) >= std::abs(value_to_add))
    carry[i] += (current_value - t) + value_to_add;
  else
    carry[i] += (value_to_add - t) + current_value;

  current_value = t;
}
template <class CoefficientType>
void KahanVector<CoefficientType>::apply_carry() {
  for (unsigned int i = 0; i < vector.size(); ++i)
    vector[i] += carry[i];

  carry = std::vector<CoefficientType>(vector.size(), 0);
}
template <class CoefficientType>
unsigned int KahanVector<CoefficientType>::size() {
  return vector.size();
}
template <class CoefficientType>
std::vector<CoefficientType> &KahanVector<CoefficientType>::data() {
  return vector;
}
template <class CoefficientType>
CoefficientType KahanVector<CoefficientType>::sum() const {
  CoefficientType out = CoefficientType(0);
  for (const auto& entry : this->vector)
    out += entry;

  return out;
}
/**
 * @brief A contiguous matrix with Neumaier-compensated summation for each entry.
 *
 * Provides stable accumulation of values with per-element carry tracking.
 *
 * @tparam CoefficientType Type of the matrix entries.
 *
 * @ingroup AssemblyContainers
 */

template <class CoefficientType> class KahanContiguousMatrix {
public:
  KahanContiguousMatrix(const unsigned int &m, const unsigned int &n);
  CoefficientType &operator()(const unsigned int &i, const unsigned int &j);
  void accumulate(const unsigned int &i, const unsigned int &j,
                  const CoefficientType &value_to_add);
  void apply_carry();


private:
  CoefficientType &get_carry(const unsigned int &i, const unsigned int &j);
  const unsigned int m;
  const unsigned int n;
  std::vector<CoefficientType> matrix_data;
  std::vector<CoefficientType> carry;
};
template <class CoefficientType>
KahanContiguousMatrix<CoefficientType>::KahanContiguousMatrix(
    const unsigned int &m, const unsigned int &n)
    : m(m), n(n) {
  matrix_data = std::vector<CoefficientType>(m * n, CoefficientType(0));
  carry = std::vector<CoefficientType>(m * n, CoefficientType(0));
}


template <class CoefficientType>
CoefficientType &
KahanContiguousMatrix<CoefficientType>::operator()(const unsigned int &i,
                                                   const unsigned int &j) {
  return matrix_data[this->n * i + j];
}
template <class CoefficientType>
void KahanContiguousMatrix<CoefficientType>::accumulate(
    const unsigned int &i, const unsigned int &j,
    const CoefficientType &value_to_add) {
  auto &current_value = this->operator()(i, j);
  // Neumaier Version
  const auto t = current_value + value_to_add;
  if (std::abs(current_value) >= std::abs(value_to_add))
    get_carry(i, j) += (current_value - t) + value_to_add;
  else
    get_carry(i, j) += (value_to_add - t) + current_value;

  current_value = t;

  //  // Original Kahan
  //    const auto y = value_to_add - get_carry(i,j);
  //    const auto t = current_value + y;
  //    get_carry(i,j) = (t-current_value) - y;
  //    current_value = t;
}
template <class CoefficientType>
CoefficientType &
KahanContiguousMatrix<CoefficientType>::get_carry(const unsigned int &i,
                                                  const unsigned int &j) {
  return carry[this->n * i + j];
}
template <class CoefficientType>
void KahanContiguousMatrix<CoefficientType>::apply_carry() {
  for (unsigned int i = 0; i < matrix_data.size(); ++i)
    matrix_data[i] += carry[i];

  carry = std::vector<CoefficientType>(matrix_data.size(), 0);
}
/**
 * @brief A contiguous matrix with higher-order Kahan summation for increased accuracy.
 *
 * Uses Neumaier compensation logic to accumulate values robustly.
 *
 * @tparam CoefficientType Type of the matrix entries.
 *
 * @ingroup AssemblyContainers
 */

template <class CoefficientType> class HighOrderKahanContiguousMatrix {
public:
  HighOrderKahanContiguousMatrix(const unsigned int &m, const unsigned int &n);
  CoefficientType &operator()(const unsigned int &i, const unsigned int &j);
  void accumulate(const unsigned int &i, const unsigned int &j,
                  const CoefficientType &value_to_add);
  void apply_carry();

private:
  CoefficientType &get_carry(const unsigned int &i, const unsigned int &j);
  const unsigned int m;
  const unsigned int n;
  std::vector<CoefficientType> matrix_data;
  std::vector<CoefficientType> carry;
};
template <class CoefficientType>
HighOrderKahanContiguousMatrix<CoefficientType>::HighOrderKahanContiguousMatrix(
    const unsigned int &m, const unsigned int &n)
    : m(m), n(n) {
  matrix_data = std::vector<CoefficientType>(m * n, CoefficientType(0));
  carry = std::vector<CoefficientType>(m * n, CoefficientType(0));
}

template <class CoefficientType>
CoefficientType &HighOrderKahanContiguousMatrix<CoefficientType>::operator()(
    const unsigned int &i, const unsigned int &j) {
  return matrix_data[this->n * i + j];
}
template <class CoefficientType>
void HighOrderKahanContiguousMatrix<CoefficientType>::accumulate(
    const unsigned int &i, const unsigned int &j,
    const CoefficientType &value_to_add) {
  auto &current_value = this->operator()(i, j);
  // Neumaier Version
  const auto t = current_value + value_to_add;
  if (std::abs(current_value) >= std::abs(value_to_add))
    get_carry(i, j) += (current_value - t) + value_to_add;
  else
    get_carry(i, j) += (value_to_add - t) + current_value;

  current_value = t;

  //  // Original Kahan
  //    const auto y = value_to_add - get_carry(i,j);
  //    const auto t = current_value + y;
  //    get_carry(i,j) = (t-current_value) - y;
  //    current_value = t;
}
template <class CoefficientType>
CoefficientType &HighOrderKahanContiguousMatrix<CoefficientType>::get_carry(
    const unsigned int &i, const unsigned int &j) {
  return carry[this->n * i + j];
}
template <class CoefficientType>
void HighOrderKahanContiguousMatrix<CoefficientType>::apply_carry() {
  for (unsigned int i = 0; i < matrix_data.size(); ++i)
    matrix_data[i] += carry[i];

  carry = std::vector<CoefficientType>(matrix_data.size(), 0);
}
/**
 * @brief A hybrid sparse-dense submatrix structure with dynamic DoF mapping.
 *
 * Allows accumulation into a matrix that is dense in one dimension and sparse in the other.
 * Also tracks which entries have been filled.
 *
 * @tparam CoefficientType Type of the matrix entries.
 *
 * @ingroup AssemblyContainers
 */

template <class CoefficientType> class DenseSubMatrix {
public:
  /**
   *
   * @param m - The sparse dimension
   * @param n - The dense dimension
   */
  DenseSubMatrix(const unsigned int &m, const unsigned int &n);

  // CoefficientType& operator[](const unsigned int& i, const unsigned int& j);
  CoefficientType &operator()(const unsigned int &i, const unsigned int &j);
  CoefficientType &at(const unsigned int &i, const unsigned int &j);
  char is_filled_at(const unsigned int &i, const unsigned int &j) const;
  char &get_is_filled_at(const unsigned int &i, const unsigned int &j);

  const CoefficientType &const_at(const unsigned int &i,
                                  const unsigned int &j) const;
  /**
   * @brief Generates a mask of which entries in the submatrix have been filled.
   *
   * For each (test DoF, trial DoF) pair, marks whether an entry exists in this submatrix.
   *
   * @tparam DoFCellType Type representing a cell with DoFs.
   * @param cell_test Test cell.
   * @param cell_trial Trial cell.
   * @param fill_data Output mask matrix.
   */
  template <class DoFCellType>
  void generate_filled_data(const DoFCellType &cell_test,
                            const DoFCellType &cell_trial,
                            ContiguousMatrix<char> *fill_data) const;
  /**
   * @brief Generates a mask of filled entries considering DoF masks.
   *
   * Same as the other overload, but respects additional masks indicating active DoFs.
   *
   * @tparam DoFCellType Type representing a cell with DoFs.
   * @param cell_test Test cell.
   * @param cell_trial Trial cell.
   * @param test_mask Active DoFs in the test cell.
   * @param trial_mask Active DoFs in the trial cell.
   * @param fill_data Output mask matrix.
   */
  template <class DoFCellType>
  void generate_filled_data(const DoFCellType &cell_test,
                            const DoFCellType &cell_trial,
                            const DoFMask &test_mask, const DoFMask &trial_mask,
                            ContiguousMatrix<char> *fill_data) const;

  void resize(const unsigned int &new_m, const unsigned int &new_n);
  unsigned int n_rows() const;
  unsigned int n_cols() const;

private:
  // While perhaps a nested vector is not ideal in terms of access efficiency
  // it is the simplest
  std::vector<std::vector<CoefficientType>> matrix;
  std::vector<std::vector<char>> is_filled;

  // As the matrix is dense in one dimension and sparse in the other,
  // this map allows us to take global indices to the correct local index
  // since the indices might not be contiguous...
  std::unordered_map<unsigned int, unsigned int> dof_mapper;
  unsigned int m = 0;
  unsigned int n = 0;
};

template <class CoefficientType>
DenseSubMatrix<CoefficientType>::DenseSubMatrix(const unsigned int &m,
                                                const unsigned int &n) {
  this->matrix = std::vector<std::vector<CoefficientType>>(
      m, std::vector<CoefficientType>(n));
  this->is_filled = std::vector<std::vector<char>>(
      m, std::vector<char>(n, false));
  this->m = m;
  this->n = n;
}
template <class CoefficientType>
CoefficientType &
DenseSubMatrix<CoefficientType>::operator()(const unsigned int &i,
                                            const unsigned int &j) {
  auto row_location = dof_mapper.find(i);

  if (row_location == dof_mapper.end()) {
    const auto res = dof_mapper.insert(
        std::pair<unsigned int, unsigned int>(i, dof_mapper.size()));
    row_location = res.first;
  }

  return matrix[row_location->second][j];
}

template <class CoefficientType>
void DenseSubMatrix<CoefficientType>::resize(const unsigned int &new_m,
                                             const unsigned int &new_n) {
  this->matrix.resize(new_m);
  this->is_filled.resize(new_m);
  for (auto &entry : this->matrix)
    entry.resize(new_n);
  for (auto &entry : this->is_filled)
    entry.resize(new_n, false);
  this->m = new_m;
  this->n = new_n;
}
template <class CoefficientType>
CoefficientType &DenseSubMatrix<CoefficientType>::at(const unsigned int &i,
                                                     const unsigned int &j) {
  auto row_location = dof_mapper.find(i);

  if (row_location == dof_mapper.end()) {
    const auto res = dof_mapper.insert(
        std::pair<unsigned int, unsigned int>(i, dof_mapper.size()));
    row_location = res.first;
  }
  auto temp_index = row_location->second;
  return (matrix[row_location->second])[j];
}
template <class CoefficientType>
const CoefficientType &
DenseSubMatrix<CoefficientType>::const_at(const unsigned int &i,
                                          const unsigned int &j) const {
  auto row_location = dof_mapper.find(i);
#ifdef DEBUG
  assert(row_location != dof_mapper.end() && "The key does not exist!");
#endif
  const auto temp_index = row_location->second;
  return (matrix[row_location->second])[j];
}
template <class CoefficientType>
char DenseSubMatrix<CoefficientType>::is_filled_at(
    const unsigned int &i, const unsigned int &j) const {
  auto row_location = dof_mapper.find(i);
  if (row_location == dof_mapper.end())
    return false;
  else
    return (is_filled[row_location->second])[j];
}
template <class CoefficientType>
char &DenseSubMatrix<CoefficientType>::get_is_filled_at(const unsigned int &i,
                                                        const unsigned int &j) {
  auto row_location = dof_mapper.find(i);
#ifdef DEBUG
  assert(row_location != dof_mapper.end() && "The key does not exist!");
#endif

  return (is_filled[row_location->second])[j];
}
template <class CoefficientType>
template <class DoFCellType>
void DenseSubMatrix<CoefficientType>::generate_filled_data(
    const DoFCellType &cell_test, const DoFCellType &cell_trial,
    ContiguousMatrix<char> *fill_data) const {
  for (unsigned int i = 0; i < cell_test.n_dofs(); ++i) {
    const auto &dof_test = cell_test.get_dof(i);
    // get the possible index here...
    auto row_location = dof_mapper.find(dof_test.global_index);

    for (unsigned int j = 0; j < cell_trial.n_dofs(); ++j) {
      const auto &dof_trial = cell_test.get_dof(i);
      if (row_location == dof_mapper.end())
        fill_data->operator()(i,j) = false;
      else
        fill_data->operator()(i,j) = is_filled[row_location->second][dof_trial.global_index];
    }
  }
}

// The same function as above, but involves DoFMasks in order to
// make the decision of whether we need to do any work  more clear...
template <class CoefficientType>
template <class DoFCellType>
void DenseSubMatrix<CoefficientType>::generate_filled_data(
    const DoFCellType &cell_test, const DoFCellType &cell_trial,
    const DoFMask &test_mask, const DoFMask &trial_mask,
    ContiguousMatrix<char> *fill_data) const
{
  for (unsigned int i = 0; i < cell_test.n_dofs(); ++i) {
    const auto &dof_test = cell_test.get_dof(i);
    // get the possible index here...
    auto row_location = dof_mapper.find(dof_test.global_index);

    for (unsigned int j = 0; j < cell_trial.n_dofs(); ++j) {
      const auto &dof_trial = cell_trial.get_dof(j);
      if (test_mask[i]==false || trial_mask[j] == false)
        fill_data->operator()(i,j) = true;
      else if (row_location == dof_mapper.end())
        fill_data->operator()(i,j) = false;
      else
        fill_data->operator()(i,j) = is_filled[row_location->second][dof_trial.global_index];
    }
  }
}
template <class CoefficientType>
unsigned int DenseSubMatrix<CoefficientType>::n_rows() const {
  return m;
}
template <class CoefficientType>
unsigned int DenseSubMatrix<CoefficientType>::n_cols() const {
  return n;
}
/**
 * @brief A dense vector supporting dynamic DoF indexing.
 *
 * Provides accumulation, access, and reset capabilities for excitation vectors or sub-assemblies.
 *
 * @tparam CoefficientType Type of the vector entries.
 *
 * @ingroup AssemblyContainers
 */

template <class CoefficientType> class DenseSubVector {
public:
  /**
   *
   * @param n - The dense dimension
   */
  DenseSubVector(const unsigned int &n);

  // CoefficientType& operator[](const unsigned int& i, const unsigned int& j);
  CoefficientType &operator()(const unsigned int &i);
  CoefficientType &at(const unsigned int &i);
  const CoefficientType &const_at(const unsigned int &i) const;

  void resize(const unsigned int &new_n);
  void zero_out();

private:
  // While perhaps a nested vector is not ideal in terms of access efficiency
  // it is the simplest
  std::vector<CoefficientType> vector;

  // As the matrix is dense in one dimension and sparse in the other,
  // this map allows us to take global indices to the correct local index
  // since the indices might not be contiguous...
  std::map<unsigned int, unsigned int> dof_mapper;
};

template <class CoefficientType>
DenseSubVector<CoefficientType>::DenseSubVector(const unsigned int &n) {
  this->vector = std::vector<CoefficientType>(n);
}
template <class CoefficientType>
CoefficientType &
DenseSubVector<CoefficientType>::operator()(const unsigned int &i) {
  auto row_location = dof_mapper.find(i);

  if (row_location == dof_mapper.end()) {
    const auto res = dof_mapper.insert(
        std::pair<unsigned int, unsigned int>(i, dof_mapper.size()));
    row_location = res.first;
  }

  return vector[row_location->second];
}

template <class CoefficientType>
void DenseSubVector<CoefficientType>::resize(const unsigned int &new_n) {
  this->vector.resize(new_n);
}
template <class CoefficientType>
CoefficientType &DenseSubVector<CoefficientType>::at(const unsigned int &i) {
  auto row_location = dof_mapper.find(i);

  if (row_location == dof_mapper.end()) {
    const auto res = dof_mapper.insert(
        std::pair<unsigned int, unsigned int>(i, dof_mapper.size()));
    row_location = res.first;
  }

  return (vector[row_location->second]);
}
template <class CoefficientType>
const CoefficientType &
DenseSubVector<CoefficientType>::const_at(const unsigned int &i) const {
  auto row_location = dof_mapper.find(i);

#ifdef DEBUG
  assert(row_location != dof_mapper.end() && "The key does not exist!");
#endif

  return (vector[row_location->second]);
}

/**
 * @brief Zeros out the vector, resetting all entries to zero.
 *
 * @tparam CoefficientType Type of the vector entries.
 */
template <class CoefficientType>
void DenseSubVector<CoefficientType>::zero_out() {
  std::fill(this->vector.begin(), this->vector.end(), CoefficientType(0));
}
DROMON_NAMESPACE_CLOSE
#endif // DROMON_SUBMATRIX_H
