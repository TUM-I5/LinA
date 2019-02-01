#ifndef DGMATRICES_H_
#define DGMATRICES_H_

#include "typedefs.h"
#include <generated_code/tensor.h>

class GlobalMatrices {
public:
  real* kDivM;
  real* kTDivM;
  real* kDivMT;
  real* kTDivMT;
  lina::tensor::FDivM::Container<real const*> FDivM;
  lina::tensor::FDivMT::Container<real const*> FDivMT;
  real* quadrature;
  
  GlobalMatrices();
  ~GlobalMatrices();

private:
  real* m_matrixMem;
};

struct LocalMatrices {
  real Astar[lina::tensor::star::size(0)];
  real Bstar[lina::tensor::star::size(1)];
  real fluxSolver[2][2][2][lina::tensor::fluxSolver::size()];
};

#endif // DGMATRICES_H_
