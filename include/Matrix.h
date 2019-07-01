#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <memory>

#include "Dimension.h"
#include "ThreadPool.h"
#include "Utility.h"

template<class T = float>
class Matrix {
  public:
    static unsigned int n_threads;
    static Matrix<T> Strassen_multiply(Matrix<T>, Matrix<T>, int = 0);

  private:
    std::vector<std::vector<T>> value;
    unsigned int n_rows, n_columns;

    bool is_square() const { return n_rows == n_columns; }
    bool is_triangular() const;
    bool is_upper_triangular() const;
    bool is_lower_triangular() const;

  public:
    Matrix(unsigned int _n_rows, unsigned int _n_columns)
      : n_rows(_n_rows), n_columns(_n_columns) {
        auto tmp = std::vector<T>(n_columns);
        value = std::vector<std::vector<T>>(n_rows, tmp);
      }

    Matrix(const Dimension2D dimension)
      : Matrix(dimension.getRow(), dimension.getColumn()) {}

    Matrix(const Matrix<T> &o) : n_rows(o.n_rows), n_columns(o.n_columns) {
      value = o.value;
    }

    Matrix(const std::vector<std::vector<T>> &_value) {
      n_rows = _value.size();
      n_columns = _value[0].size();
      value = _value;
    }

    ~Matrix() {}

    Matrix<T> operator+(const Matrix &) const;
    Matrix<T> operator-(const Matrix<T> &) const;
    Matrix<T> operator*(const Matrix<T> &)const;
    bool operator==(const Matrix<T> &) const;
    bool operator!=(const Matrix<T> &) const;

    template<class U>
      friend std::ostream& operator<< (std::ostream & os, const Matrix<U> &);

    Dimension2D getDimension() const { return Dimension2D({n_rows, n_columns}); }
    unsigned int getNRows() const { return n_rows; }
    unsigned int getNColumns() const { return n_columns; }
    std::vector<std::vector<T>> getValues() const { return value; }
    void set(unsigned int i, unsigned int j, T x) { value[i][j] = x; }

    bool same_dimension(const Matrix<T> &o) const {
      return n_rows == o.n_rows && n_columns == o.n_columns;
    }

    Matrix<T> sub_matrix(unsigned int, unsigned int, unsigned int,
        unsigned int) const;

    T det(int = 0) const;

    std::pair<Matrix, Matrix> lu_decomposition() const;

  private:
    void _add_thread(const Matrix<T> &, unsigned int, unsigned int, Matrix<T> &) const;
    void _subtract_thread(const Matrix<T> &, unsigned int, unsigned int,
        Matrix<T> &) const;
    void _multiply_thread(const Matrix<T> &, unsigned int, unsigned int,
        Matrix<T> &) const;
    T _laplace_determinant() const;
    T _division_free() const;
    T _det_triangular_matrix() const;
};

template<class T>
std::ostream& operator<<(std::ostream &os, const Matrix<T> &matrix) {
  for (unsigned int i = 0; i < matrix.n_rows; i++) {
    os << ((i == 0) ? "[" : " ");
    for (unsigned int j = 0; j < matrix.n_columns; j++) {
      os << matrix.value[i][j];
      if (j + 1 == matrix.n_columns)
        if (i + 1 == matrix.n_rows)
          os << "]";
        else
          os << "\n";
      else
        os << ", ";
    }
  }
  return os;
}

template<class T>
Matrix<T> Matrix<T>::operator+(const Matrix<T> &o) const {
  assert(n_rows == o.getNRows());
  assert(n_columns == o.getNColumns());

  Matrix<T> ans(n_rows, n_columns);
  unsigned int total_load = 1ll * n_rows * n_columns;
  unsigned int start = 0;
  unsigned int load_per_thread = total_load / (n_threads + 1);
  for (unsigned int i = 0; i < n_threads; i++) {
    unsigned int load_for_thread =
      load_per_thread + (i < total_load % (n_threads + 1) ? 1 : 0);
    std::thread(&Matrix::_add_thread, *this, std::ref(o), start,
        start + load_for_thread, std::ref(ans))
      .join();
    start += load_for_thread;
  }
  _add_thread(o, start, total_load, ans);
  return ans;
}

template<class T>
Matrix<T> Matrix<T>::operator-(const Matrix<T> &o) const {
  assert(n_rows == o.getNRows());
  assert(n_columns == o.getNColumns());

  Matrix<T> ans(n_rows, n_columns);
  unsigned int total_load = 1ll * n_rows * n_columns;
  unsigned int start = 0;
  unsigned int load_per_thread = total_load / (n_threads + 1);
  for (unsigned int i = 0; i < n_threads; i++) {
    unsigned int load_for_thread =
      load_per_thread + (i < total_load % (n_threads + 1) ? 1 : 0);
    std::thread(&Matrix::_subtract_thread, *this, std::ref(o), start,
        start + load_for_thread, std::ref(ans))
      .join();
    start += load_for_thread;
  }
  _subtract_thread(o, start, total_load, ans);
  return ans;
}

template<class T>
Matrix<T> Matrix<T>::operator*(const Matrix<T> &o) const {
  assert(n_columns == o.getNRows());

  Matrix<T> ans(n_rows, o.getNColumns());
  unsigned int total_load = 1ll * n_rows * o.getNColumns();
  unsigned int start = 0;
  unsigned int load_per_thread = total_load / (n_threads + 1);
  for (unsigned int i = 0; i < n_threads; i++) {
    unsigned int load_for_thread =
      load_per_thread + (i < total_load % (n_threads + 1) ? 1 : 0);
    std::thread(&Matrix::_multiply_thread, *this, std::ref(o), start,
        start + load_for_thread, std::ref(ans))
      .join();
    start += load_for_thread;
  }
  _multiply_thread(o, start, total_load, ans);
  return ans;
}

template<class T>
bool Matrix<T>::operator==(const Matrix<T> & o) const {
  if (n_columns != o.n_columns) 
    return false;
  if (n_rows != o.n_rows)
    return false;
  for(unsigned int i=0; i<n_rows; i++)
    for(unsigned int j=0; j<n_columns; j++)
      if (value[i][j] != o.value[i][j])
        return false;
  return true;
}

template<class T>
void Matrix<T>::_add_thread(const Matrix<T> &o, unsigned int start, unsigned int end,
    Matrix<T> &ans) const {
  for (unsigned int cell = start; cell < end; cell++) {
    unsigned int row = cell / n_columns;
    unsigned int col = cell % n_columns;
    ans.value[row][col] = value[row][col] + o.value[row][col];
  }
}

template<class T>
void Matrix<T>::_subtract_thread(const Matrix<T> &o, unsigned int start,
    unsigned int end, Matrix<T> &ans) const {
  for (unsigned int cell = start; cell < end; cell++) {
    unsigned int row = cell / n_columns;
    unsigned int col = cell % n_columns;
    ans.value[row][col] = value[row][col] - o.value[row][col];
  }
}

template<class T>
void Matrix<T>::_multiply_thread(const Matrix<T> &o, unsigned int start,
    unsigned int end, Matrix<T> &ans) const {
  for (unsigned int cell = start; cell < end; cell++) {
    unsigned int row = cell / n_columns;
    unsigned int col = cell % n_columns;
    for (unsigned int i = 0; i < n_columns; i++)
      ans.value[row][col] += value[row][i] * o.value[i][col];
  }
}

template<class T>
Matrix<T> Matrix<T>::Strassen_multiply(Matrix<T> a, Matrix<T> b, int max_recursive) {
  assert(a.is_square());
  assert(a.same_dimension(b));
  if (a.n_rows % 2) return a * b;
  if (max_recursive <= 0) return a * b;
  unsigned int half = a.n_rows / 2;
  Matrix<T> a00 = a.sub_matrix(0, half, 0, half);
  Matrix<T> a01 = a.sub_matrix(0, half, half, a.n_columns);
  Matrix<T> a10 = a.sub_matrix(half, a.n_rows, 0, half);
  Matrix<T> a11 = a.sub_matrix(half, a.n_rows, half, a.n_columns);
  Matrix<T> b00 = b.sub_matrix(0, half, 0, half);
  Matrix<T> b01 = b.sub_matrix(0, half, half, b.n_columns);
  Matrix<T> b10 = b.sub_matrix(half, b.n_rows, 0, half);
  Matrix<T> b11 = b.sub_matrix(half, b.n_rows, half, b.n_columns);
  Matrix<T> m0 = Strassen_multiply(a00 + a11, b00 + b11, max_recursive - 1);
  Matrix<T> m1 = Strassen_multiply(a10 + a11, b00, max_recursive - 1);
  Matrix<T> m2 = Strassen_multiply(a00, b01 - b11, max_recursive - 1);
  Matrix<T> m3 = Strassen_multiply(a11, b10 - b00, max_recursive - 1);
  Matrix<T> m4 = Strassen_multiply(a00 + a01, b11, max_recursive - 1);
  Matrix<T> m5 = Strassen_multiply(a10 - a00, b00 + b01, max_recursive - 1);
  Matrix<T> m6 = Strassen_multiply(a01 - a11, b10 + b11, max_recursive - 1);
  Matrix<T> c00 = m0 + m3 - m4 + m6;
  Matrix<T> c01 = m2 + m4;
  Matrix<T> c10 = m1 + m3;
  Matrix<T> c11 = m0 - m1 + m2 + m5;
  Matrix<T> ans(a.n_rows, a.n_columns);
  for (unsigned int i = 0; i < half; i++)
    for (unsigned int j = 0; j < half; j++) {
      ans.value[i][j] = c00.value[i][j];
      ans.value[i + half][j] = c10.value[i][j];
      ans.value[i][j + half] = c01.value[i][j];
      ans.value[i + half][j + half] = c11.value[i][j];
    }
  return ans;
}

template<class T>
Matrix<T> Matrix<T>::sub_matrix(unsigned int from_row, unsigned int to_row,
    unsigned int from_col, unsigned int to_col) const {
  Matrix<T> ans(to_row - from_row, to_col - from_col);
  for (unsigned int i = from_row; i < to_row; i++)
    for (unsigned int j = from_col; j < to_col; j++)
      ans.value[i - from_row][j - from_col] = value[i][j];
  return ans;
}

template<class T>
T Matrix<T>::det(int option) const {
  assert(is_square());
  if (is_triangular())
    return _det_triangular_matrix();
  switch (option) {
    case 0:
      return _laplace_determinant();
    case 1:
      return _division_free();
    case 2: 
      {
        auto lu = lu_decomposition();
        return lu.first.det() * lu.second.det();
      }
    default: 
      throw "Option invalid";
  }
}

template<class T>
T Matrix<T>::_laplace_determinant() const {
  if (n_rows == 2)
    return value[0][0] * value[1][1] - value[0][1] * value[1][0];
  Matrix<T> tmp(n_rows-1, n_columns-1);
  for(unsigned int i=1; i<n_rows; i++)
    for(unsigned int j=1; j<n_columns; j++)
      tmp.value[i-1][j-1] = value[i][j];
  int sign = 1;
  T ans = tmp._laplace_determinant() * value[0][0];;
  for(unsigned int i=1; i<n_rows; i++) {
    for(unsigned int j=1; j<n_columns; j++)
      tmp.value[i-1][j-1] = value[i-1][j];
    sign = -sign;
    ans += sign * tmp._laplace_determinant() * value[i][0];
  }
  return ans;
}

template<class T>
T Matrix<T>::_division_free() const {
  const unsigned int n = n_columns;
  std::vector<std::vector<std::vector<std::vector<T>>>> f = vector_4d<T>(n, n+1, n+1, 2, 0);
  for(unsigned int i=0; i<n; i++) 
    f[0][i][i][1] = 1;
  for(unsigned int i=0; i<n-1; i++) 
    for(unsigned int j=0; j<n; j++) 
      for(unsigned int k=j; k<n; k++)
        for(unsigned int s=0; s<2; s++) {
          for(unsigned int h=j+1; h<n; h++)
            f[i+1][h][h][1-s] += f[i][j][k][s] * value[k][j];
          for(unsigned int h=j+1; h<n; h++)
            f[i+1][j][h][s] += f[i][j][k][s] * value[k][h];
        }
  T ans = 0;
  for(unsigned int i=0; i<n; i++) 
    for(unsigned int j=i; j<n; j++)
      for(unsigned int s=0; s<2; s++)
        ans += f[n-1][i][j][s] * value[j][i] * (s == 0 ? -1 : 1);
  return (n % 2 == 0) ? -ans : ans;
}

template<class T>
std::pair<Matrix<T>, Matrix<T>> Matrix<T>::lu_decomposition() const {
  assert(is_square());
  unsigned int n = n_rows;
  Matrix<T> l(n, n);
  Matrix<T> u(n, n);
  for(unsigned int i=0; i<n; i++) {
    u.value[i][i] = 1;
    T tmp = 0;
    for(unsigned int j=0; j<i; j++) {
      tmp = 0;
      for(unsigned int k=0; k<j; k++) 
        tmp += l.value[i][k] * u.value[k][j];
      l.value[i][j] = value[i][j] - tmp;
      tmp = 0;
      for(unsigned int k=0; k<j; k++) 
        tmp += l.value[j][k] * u.value[k][i];
      u.value[j][i] = (value[j][i] - tmp) / l.value[j][j];
    }
    tmp = 0;
    for(unsigned int j=0; j<i; j++) 
      tmp += l.value[i][j] * u.value[j][i];
    l.value[i][i] = value[i][i] - tmp;
  }
  return std::make_pair(l, u);
}

template<class T>
bool Matrix<T>::is_triangular() const {
  return is_upper_triangular() || is_lower_triangular();
}

template<class T>
bool Matrix<T>::is_upper_triangular() const {
  assert(is_square());
  for(unsigned long i=0; i<n_rows; i++)
    for(unsigned long j=0; j<i; j++) 
      if (value[i][j] != 0)
        return false;
  return true;
}

template<class T>
bool Matrix<T>::is_lower_triangular() const {
  assert(is_square());
  for(unsigned long i=0; i<n_rows; i++)
    for(unsigned long j=i+1; j<n_columns; j++) 
      if (value[i][j] != 0)
        return false;
  return true;
}

template<class T>
T Matrix<T>::_det_triangular_matrix() const {
  assert(is_triangular());
  T det = 1;
  for(unsigned long i=0; i<n_rows; i++) 
    det *= value[i][i];
  return det;
}
