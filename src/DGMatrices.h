#ifndef DGMATRICES_H_
#define DGMATRICES_H_

#include <generated_code/tensor.h>

class GlobalMatrices {
public:
  double* kDivM;
  double* kTDivM;
  double* kDivMT;
  double* kTDivMT;
  lina::tensor::FDivM::Container<double const*> FDivM;
  lina::tensor::FDivMT::Container<double const*> FDivMT;
  double* quadrature;
  
  GlobalMatrices();
  ~GlobalMatrices();

private:
  double* m_matrixMem;
};

struct LocalMatrices {
  double Astar[lina::tensor::star::size(0)];
  double Bstar[lina::tensor::star::size(1)];
  double Cstar[lina::tensor::star::size(2)];
  double fluxSolver[3][2][2][lina::tensor::fluxSolver::size()];
};

#endif // DGMATRICES_H_
