#include "Matrix.h"

#include <thread>

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

Matrix Matrix::Strassen_multiply(Matrix a, Matrix b, int max_recursive) {
  assert(a.is_square());
  assert(a.same_dimension(b));
  if (a.n_rows % 2) return a * b;
  if (max_recursive <= 0) return a * b;
  unsigned int half = a.n_rows / 2;
  Matrix a00 = a.sub_matrix(0, half, 0, half);
  Matrix a01 = a.sub_matrix(0, half, half, a.n_columns);
  Matrix a10 = a.sub_matrix(half, a.n_rows, 0, half);
  Matrix a11 = a.sub_matrix(half, a.n_rows, half, a.n_columns);
  Matrix b00 = b.sub_matrix(0, half, 0, half);
  Matrix b01 = b.sub_matrix(0, half, half, b.n_columns);
  Matrix b10 = b.sub_matrix(half, b.n_rows, 0, half);
  Matrix b11 = b.sub_matrix(half, b.n_rows, half, b.n_columns);
  Matrix m0 = Strassen_multiply(a00 + a11, b00 + b11, max_recursive - 1);
  Matrix m1 = Strassen_multiply(a10 + a11, b00, max_recursive - 1);
  Matrix m2 = Strassen_multiply(a00, b01 - b11, max_recursive - 1);
  Matrix m3 = Strassen_multiply(a11, b10 - b00, max_recursive - 1);
  Matrix m4 = Strassen_multiply(a00 + a01, b11, max_recursive - 1);
  Matrix m5 = Strassen_multiply(a10 - a00, b00 + b01, max_recursive - 1);
  Matrix m6 = Strassen_multiply(a01 - a11, b10 + b11, max_recursive - 1);
  Matrix c00 = m0 + m3 - m4 + m6;
  Matrix c01 = m2 + m4;
  Matrix c10 = m1 + m3;
  Matrix c11 = m0 - m1 + m2 + m5;
  Matrix ans(a.n_rows, a.n_columns);
  for (unsigned int i = 0; i < half; i++)
    for (unsigned int j = 0; j < half; j++) {
      ans.value[i][j] = c00.value[i][j];
      ans.value[i + half][j] = c10.value[i][j];
      ans.value[i][j + half] = c01.value[i][j];
      ans.value[i + half][j + half] = c11.value[i][j];
    }
  return ans;
}

Matrix Matrix::sub_matrix(unsigned int from_row, unsigned int to_row,
                          unsigned int from_col, unsigned int to_col) const {
  Matrix ans(to_row - from_row, to_col - from_col);
  for (unsigned int i = from_row; i < to_row; i++)
    for (unsigned int j = from_col; j < to_col; j++)
      ans.value[i - from_row][j - from_col] = value[i][j];
  return ans;
}
