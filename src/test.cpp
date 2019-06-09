#include <cmath>
#include <cstdio>
#include <iostream>
#include <memory>
#include <vector>

#include "Matrix.h"
#include "ThreadPool.h"

using namespace std;

unsigned int Matrix::n_threads;

Matrix random_matrix(int n, int m) {
  srand(time(0));
  vector<vector<float>> ans;
  for (int i = 0; i < n; i++) {
    vector<float> tmp;
    for (int j = 0; j < m; j++) tmp.push_back(rand() % 10);
    ans.push_back(tmp);
  }
  return Matrix(ans);
}

void run() {
  srand(time(0));
  cout << "Number of threads: ";
  cin >> Matrix::n_threads;
  cout << "Size of the matrix: ";
  int n;
  cin >> n;
  printf("Init value \n");
  printf("Create matrix a\n");
  Matrix ma = random_matrix(n, n);
  printf("Create matrix b\n");
  Matrix mb = random_matrix(n, n);
  auto start = std::chrono::system_clock::now();
  Matrix me = ma * mb;
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end - start;
  cout << "Multithread: " << elapsed_seconds.count() << endl;

  start = std::chrono::system_clock::now();
  Matrix mg = Matrix::Strassen_multiply(ma, mb, 1);
  end = std::chrono::system_clock::now();
  elapsed_seconds = end - start;
  cout << "Strassen once: " << elapsed_seconds.count() << endl;

  start = std::chrono::system_clock::now();
  Matrix mf = Matrix::Strassen_multiply(ma, mb, 100);
  end = std::chrono::system_clock::now();
  elapsed_seconds = end - start;
  cout << "Strassen til can't: " << elapsed_seconds.count() << endl;
}

void test_det() {
  srand(time(0));
  cout << "Number of threads: ";
  cin >> Matrix::n_threads;
  cout << "Size of the matrix: ";
  int n;
  cin >> n;
  printf("Init value \n");
  printf("Create matrix a\n");
  Matrix ma = random_matrix(n, n);
  cout << ma << endl;
  cout << int(ma.det()) << endl;
  cout << int(ma.det(1)) << endl;
}

int main() {
  test_det();
}
