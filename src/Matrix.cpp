#include "Matrix.h"

std::ostream& operator<<(std::ostream& os, const Matrix& matrix) {
  for(unsigned int i = 0; i < matrix.n_rows; i++) {
    os << ((i == 0) ? "[" : " ");
    for(unsigned int j = 0; j < matrix.n_columns; j++) {
      os << matrix.value[i][j];
      if (j+1 == matrix.n_columns) 
        if (i+1 == matrix.n_rows) 
          os << "]\n";
        else 
          os << "\n";
      else
        os << ", ";
    }
  }
  return os;
}
