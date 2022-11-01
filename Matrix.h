#pragma once
#include <iostream>
#include <vector>

using std::vector;
using std::cout;

template<typename T>
class Matrix {
private:
  vector<vector<T>> data;
  size_t row_size, col_size;
public:
  Matrix() {};
  ~Matrix() {};
  Matrix(size_t s): row_size(s), col_size(s) {
    this->data.resize(s, vector<T>(s));
  };
  Matrix(size_t rs, size_t cs): row_size(rs), col_size(cs) {
    this->data.resize(rs, vector<T>(cs));
  };
  Matrix(const vector<vector<T>> &data) : row_size(data.size()), col_size(data[0].size()) {
    this->data = data;
  }
  Matrix(size_t s, const vector<vector<T>> &data) : row_size(s), col_size(s) {
    if (data.size() != s * s) throw std::runtime_error("uncorrect vector size");
    this->data = data;
  }
  Matrix(size_t rows, size_t cols, const vector<vector<T>> &data) : row_size(rows), col_size(cols) {
    if (data.size() != rows * cols) throw std::runtime_error("uncorrect vector size");
    this->data = data;
  }

  size_t n_rows() { return row_size; }
  size_t n_cols() { return col_size; }

  friend std::ostream &operator<<(std::ostream &os, const Matrix &matrix) {
    for (const auto& row: matrix.data) {
      for (const auto& item: row)
        os << item << "\t";
      os << "\n";
    }
    return os;
  }

  Matrix fill(T value) {
    for (size_t i = 0; i < row_size; i++)
      for (size_t j = 0; j < col_size; j++)
        data[i][j] = value;
    return *this;
  }

  T &operator()(size_t i, size_t j) {
    if (j >= col_size || i >= row_size) throw std::runtime_error("Exiting the boundaries of the vectors");
    return data[i][j];
  }

  T operator()(size_t i, size_t j) const {
    if (j >= col_size || i >= row_size) throw std::runtime_error("Exiting the boundaries of the vectors");
    return data[i][j];
  }

  Matrix operator*(const Matrix &b) const {
    if (col_size != b.row_size) throw std::runtime_error("The mismatch of the sizes of matrices");
    Matrix result(row_size, b.col_size);
    result.fill(0);

    for (size_t i = 0; i < result.n_rows(); i++) {
      for (size_t k = 0; k < col_size; k++) {
        double tmp = operator()(i, k);
        for (size_t j = 0; j < result.n_cols(); j++)
          result.data[i][j] += tmp * b(k, j);
      }
    }
    return result;
  }
};