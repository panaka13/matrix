#include "Matrix.h"

float Matrix::_laplace_determinant() const {
  if (n_rows == 2)
    return value[0][0] * value[1][1] - value[0][1] * value[1][0];
  Matrix tmp(n_rows-1, n_columns-1);
  for(unsigned int i=1; i<n_rows; i++)
    for(unsigned int j=1; j<n_columns; j++)
      tmp.value[i-1][j-1] = value[i][j];
  int sign = 1;
  float ans = tmp._laplace_determinant() * value[0][0];;
  for(unsigned int i=1; i<n_rows; i++) {
    for(unsigned int j=1; j<n_columns; j++)
      tmp.value[i-1][j-1] = value[i-1][j];
    sign = -sign;
    ans += sign * tmp._laplace_determinant() * value[i][0];
  }
  return ans;
}


