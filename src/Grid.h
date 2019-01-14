#ifndef GRID_H_
#define GRID_H_

#include <cstring>
#include <cstdlib>
#include <iostream>

template<typename T>
class Grid {
public:
  Grid(int X, int Y);
  ~Grid();
  
  /// Implements periodic boundary conditions
  inline int periodicIndex(int i, int N) {
    int pi;
    if (i >= 0) {
      pi = i % N;
    } else {
      pi = N-1 - (-i-1)%N;
    }
    return pi;
  }
  
  inline T& get(int x, int y) {
    return m_data[periodicIndex(y, m_Y) * m_X + periodicIndex(x, m_X)];
  }
  
  inline int X() const {
    return m_X;
  }
  
  inline int Y() const {
    return m_Y;
  }

private:
  int m_X;
  int m_Y;
  T* m_data;
};

template<typename T>
Grid<T>::Grid(int X, int Y)
  : m_X(X), m_Y(Y)
{
  int err = posix_memalign(reinterpret_cast<void**>(&m_data), ALIGNMENT, X*Y*sizeof(T));
  if (err) {
    std::cerr << "Failed to allocate " << X*Y*sizeof(T) << " bytes in " << __FILE__ << std::endl;
    exit(EXIT_FAILURE);
  }
  #pragma omp parallel for collapse(2)
  for (int y = 0; y < m_Y; ++y) {
    for (int x = 0; x < m_X; ++x) {
      T& data = get(x, y);
      memset(&data, 0, sizeof(T));
    }
  }
}

template<typename T>
Grid<T>::~Grid() {
  free(m_data);
}

#endif // GRID_H_
