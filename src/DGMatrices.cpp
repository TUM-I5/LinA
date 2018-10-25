#include "DGMatrices.h"
#include <yateto.h>
#include <generated_code/init.h>

GlobalMatrices::GlobalMatrices() {
  size_t alignedReals = yateto::alignedReals<double>(ALIGNMENT);
  unsigned globalMatrixMemSize = 0;
  globalMatrixMemSize += yateto::alignedUpper(lina::init::kDivM::size(), alignedReals);
  globalMatrixMemSize += yateto::alignedUpper(lina::init::kDivMT::size(), alignedReals);
  globalMatrixMemSize += yateto::computeFamilySize<lina::init::FDivM>(alignedReals);

  posix_memalign(reinterpret_cast<void**>(&m_matrixMem), ALIGNMENT, globalMatrixMemSize * sizeof(double));

  double* mem = m_matrixMem;
  yateto::copyTensorToMemAndSetPtr<lina::init::kDivM,  double>(mem, kDivM, ALIGNMENT);
  yateto::copyTensorToMemAndSetPtr<lina::init::kDivMT, double>(mem, kDivMT, ALIGNMENT);
  yateto::copyFamilyToMemAndSetPtr<lina::init::FDivM,  double>(mem, FDivM, ALIGNMENT);

  for (unsigned i = 0; i < lina::init::kDivMT::size(); ++i) {
    kDivMT[i] *= -1.0;
  }
}

GlobalMatrices::~GlobalMatrices() {
  free(m_matrixMem);
}
