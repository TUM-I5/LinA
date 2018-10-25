#ifndef DGMATRICES_H_
#define DGMATRICES_H_

#include <generated_code/tensor.h>

class GlobalMatrices {
public:
  double* kDivM;
  double* kDivMT;
  lina::tensor::FDivM::Container<double const*> FDivM;
  
  GlobalMatrices();
  ~GlobalMatrices();

private:
  double* m_matrixMem;
};

struct LocalMatrices {
  double Astar[lina::tensor::star::size(0)];
  double Bstar[lina::tensor::star::size(1)];
  double fluxSolver[2][2][2][lina::tensor::fluxSolver::size()];
};

#endif // DGMATRICES_H_
