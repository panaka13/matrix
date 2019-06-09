#include "Matrix.h"
#include "Utility.h"

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

float Matrix::_division_free() const {
  const unsigned long n = n_columns;
  std::vector<std::vector<std::vector<std::vector<float>>>> f;
  f = vector_4d<float>(n, n+1, n+1, 2, 0);
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
  float ans = 0;
  for(unsigned int i=0; i<n; i++) 
    for(unsigned int j=i; j<n; j++)
      for(unsigned int s=0; s<2; s++)
        ans += f[n-1][i][j][s] * value[j][i] * (s == 0 ? -1 : 1);
  return (n % 2 == 0) ? -ans : ans;
}

