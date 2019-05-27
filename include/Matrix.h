#pragma once

#include <algorithm>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <memory>

#include "Dimension.h"
#include "ThreadPool.h"

class Matrix {

public:
  static unsigned int n_threads;

private:
  float **value;
  unsigned int n_rows, n_columns;

public:
  Matrix(unsigned int _n_rows, unsigned int _n_columns)
      : n_rows(_n_rows), n_columns(_n_columns) {
    value = new float *[n_rows];
    for (unsigned int i = 0; i < n_rows; i++)
      value[i] = new float[n_columns];
  }

  Matrix(const Dimension2D dimension)
      : Matrix(dimension.getRow(), dimension.getColumn()) {}

  Matrix(const Matrix &o) : n_rows(o.n_rows), n_columns(o.n_columns) {
    value = new float *[n_rows];
    for (unsigned int i = 0; i < n_rows; i++) {
      value[i] = new float[n_columns];
      std::copy(o.value[i], o.value[i] + n_columns, value[i]);
    }
  }

  Matrix(const std::vector<std::vector<float>> &_value) {
    n_rows = _value.size();
    n_columns = 0;
    for (auto row : _value)
      n_columns = std::max((unsigned int)row.size(), n_columns);
    value = new float *[n_rows];
    for (unsigned int i = 0; i < n_rows; i++) {
      value[i] = new float[n_columns];
      std::copy(_value[i].begin(), _value[i].end(), value[i]);
    }
  }

  Matrix(unsigned int _n_rows, unsigned int _n_columns, float **_value)
      : n_rows(_n_rows), n_columns(_n_columns) {
    value = new float *[n_rows];
    for (unsigned int i = 0; i < n_rows; i++) {
      value[i] = new float[n_columns];
      std::copy(_value[i], _value[i] + n_columns, value[i]);
    }
  }

  ~Matrix() {}

  Matrix operator+(const Matrix &) const;
  Matrix operator-(const Matrix &) const;
  Matrix operator*(const Matrix &)const;
  bool operator==(const Matrix &) const;
  bool operator!=(const Matrix &) const;
  friend std::ostream &operator<<(std::ostream &, const Matrix &);

  Dimension2D getDimension() const { return Dimension2D({n_rows, n_columns}); }
  unsigned int getNRows() const { return n_rows; }
  unsigned int getNColumns() const { return n_columns; }
  float **getValues() const { return value; }

private:
  void _add_thread(const Matrix &, unsigned int, unsigned int, Matrix &) const;
  void _subtract_thread(const Matrix &, unsigned int, unsigned int,
                        Matrix &) const;
  void _multiply_thread(const Matrix &, unsigned int, unsigned int,
                        Matrix &) const;
};
