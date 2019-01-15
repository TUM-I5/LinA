#ifndef GRID_H_
#define GRID_H_

#include <cstring>
#include <cstdlib>
#include <iostream>

template<typename T>
class Grid {
public:
  Grid(int X, int Y, int Z);
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
  
  inline T& get(int x, int y, int z) {
    return m_data[periodicIndex(z, m_Z) * m_X * m_Z + periodicIndex(y, m_Y) * m_X + periodicIndex(x, m_X)];
  }
  
  inline int X() const {
    return m_X;
  }
  
  inline int Y() const {
    return m_Y;
  }
  
  inline int Z() const {
    return m_Z;
  }

private:
  int m_X;
  int m_Y;
  int m_Z;
  T* m_data;
};

template<typename T>
Grid<T>::Grid(int X, int Y, int Z)
  : m_X(X), m_Y(Y), m_Z(Z)
{
  int err = posix_memalign(reinterpret_cast<void**>(&m_data), ALIGNMENT, X*Y*Z*sizeof(T));
  if (err) {
    std::cerr << "Failed to allocate " << X*Y*Z*sizeof(T) << " bytes in " << __FILE__ << std::endl;
    exit(EXIT_FAILURE);
  }
  #pragma omp parallel for collapse(3)
  for (int z = 0; z < m_Z; ++z) {
    for (int y = 0; y < m_Y; ++y) {
      for (int x = 0; x < m_X; ++x) {
        T& data = get(x, y, z);
        memset(&data, 0, sizeof(T));
      }
    }
  }
}

template<typename T>
Grid<T>::~Grid() {
  free(m_data);
}

#endif // GRID_H_
