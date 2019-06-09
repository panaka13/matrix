#include <vector>

template<class T>
std::vector<T> initialize_vector(unsigned int size, T& t) {
  return std::vector<T>(size, t);
}

template<class T>
std::vector<std::vector<std::vector<std::vector<T>>>> vector_4d(
    unsigned int d1,
    unsigned int d2,
    unsigned int d3,
    unsigned int d4,
    T t) {
  auto v4 = initialize_vector<T>(d4, t);
  auto v3 = initialize_vector<std::vector<T>>(d3, v4);
  auto v2 = initialize_vector<std::vector<std::vector<T>>>(d2, v3);
  return initialize_vector<std::vector<std::vector<std::vector<T>>>>(d1, v2);
}
