#ifndef GEMM_H_
#define GEMM_H_

#include "typedefs.h"

/// Generalized matrix-matrix multiplication for column-major storage
void DGEMM(unsigned M, unsigned N, unsigned K, real alpha, real const* A, unsigned ldA, real const* B, unsigned ldB, real beta, double* C, unsigned ldC);

#endif // GEMM_H_
