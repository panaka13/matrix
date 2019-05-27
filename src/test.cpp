#include <cstdio>
#include <iostream>
#include <memory>
#include <vector>

#include "Matrix.h"
#include "ThreadPool.h"

using namespace std;

unsigned int Matrix::n_threads;

int main() {
  srand(time(0));
  cout << "Number of threads: ";
  cin >> Matrix::n_threads;
  cout << "Size of the matrix: ";
  int n;
  cin >> n;
  printf("Init value \n");
  printf("Create matrix a\n");
  Matrix ma(n, n);
  printf("Create matrix b\n");
  Matrix mb(n, n);
  auto start = std::chrono::system_clock::now();
  printf("Create matrix c\n");
  Matrix mc = ma + mb;
  Matrix md = ma - mb;
  Matrix me = ma * mb;
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end - start;
  cout << elapsed_seconds.count();
}
