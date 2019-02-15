#include "GEMM.h"

void DGEMM(unsigned M, unsigned N, unsigned K, real alpha, real const* A, unsigned ldA, real const* B, unsigned ldB, real beta, real* C, unsigned ldC) {
  for (unsigned j = 0; j < N; ++j) {
    for (unsigned i = 0; i < M; ++i) {
      real cij = 0.0;
      for (unsigned k = 0; k < K; ++k) {
        cij += A[k*ldA + i] * B[j*ldB + k];
      }
      C[j*ldC + i] = alpha * cij + beta * C[j*ldC + i];
    }
  }
}
