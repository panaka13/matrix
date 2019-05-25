#include <cstdio>
#include <iostream>
#include <vector>

#include "Matrix.h"

using namespace std;

int main() {
  vector<vector<float>> a{{1, 2}, {3, 4}};
  Matrix ma(a);
  Matrix mb(ma);
  cout << ma << endl;
  cout << mb << endl;
}
