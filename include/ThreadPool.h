#pragma once

#include <thread>

struct ThreadPool {
  int n_threads;

  ThreadPool(int n) : n_threads(n) {}
};
