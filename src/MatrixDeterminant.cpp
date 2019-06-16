#include "Matrix.h"
#include "Utility.h"

template class Matrix<int>;

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
  const unsigned long n = n_columns;
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

