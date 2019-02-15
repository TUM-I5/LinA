#ifndef TYPEDEFS_H_
#define TYPEDEFS_H_

#include <cmath>
#include "constants.h"
#include <generated_code/tensor.h>

#if REAL_SIZE == 8
typedef double real;
#elif REAL_SIZE == 4
typedef float real;
#else
#  error REAL_SIZE not supported.
#endif

typedef real DegreesOfFreedom[lina::tensor::Q::size()];
typedef real EdgeDOFs[lina::tensor::Q1::size()];

struct GlobalConstants {
  double hx;
  double hy;
  double hz;
  int X;
  int Y;
  int Z;
  double maxTimestep;
  double endTime;
};

struct Material {
  double K0;
  double rho0;
  
  inline double wavespeed() const { return sqrt(K0/rho0); }
};

struct SourceTerm {
  SourceTerm() : x(-1), y(-1), z(-1) {} // -1 == invalid
  int x;
  int y;
  int z;
  real phi[NUMBER_OF_BASIS_FUNCTIONS];
  double (*antiderivative)(double);
  int quantity;
};

#endif // TYPEDEFS_H_
