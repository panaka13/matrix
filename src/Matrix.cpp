#include <thread>

#include "Matrix.h"

std::ostream &operator<<(std::ostream &os, const Matrix &matrix) {
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

Matrix Matrix::operator+(const Matrix &o) const {
  assert(n_rows == o.getNRows());
  assert(n_columns == o.getNColumns());

  Matrix ans(n_rows, n_columns);
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

Matrix Matrix::operator-(const Matrix &o) const {
  assert(n_rows == o.getNRows());
  assert(n_columns == o.getNColumns());

  Matrix ans(n_rows, n_columns);
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

Matrix Matrix::operator*(const Matrix &o) const {
  assert(n_columns == o.getNRows());

  Matrix ans(n_rows, o.getNColumns());
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

void Matrix::_add_thread(const Matrix &o, unsigned int start, unsigned int end,
                         Matrix &ans) const {
  for (unsigned int cell = start; cell < end; cell++) {
    unsigned int row = cell / n_columns;
    unsigned int col = cell % n_columns;
    ans.value[row][col] = value[row][col] + o.value[row][col];
  }
}

void Matrix::_subtract_thread(const Matrix &o, unsigned int start,
                              unsigned int end, Matrix &ans) const {
  for (unsigned int cell = start; cell < end; cell++) {
    unsigned int row = cell / n_columns;
    unsigned int col = cell % n_columns;
    ans.value[row][col] = value[row][col] - o.value[row][col];
  }
}

void Matrix::_multiply_thread(const Matrix &o, unsigned int start,
                              unsigned int end, Matrix &ans) const {
  for (unsigned int cell = start; cell < end; cell++) {
    unsigned int row = cell / n_columns;
    unsigned int col = cell % n_columns;
    for (unsigned int i = 0; i < n_columns; i++)
      ans.value[row][col] += value[row][i] * o.value[i][col];
  }
}
