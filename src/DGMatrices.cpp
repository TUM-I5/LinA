#include "DGMatrices.h"
#include <yateto.h>
#include <generated_code/init.h>

GlobalMatrices::GlobalMatrices() {
  size_t alignedReals = yateto::alignedReals<double>(ALIGNMENT);
  unsigned globalMatrixMemSize = 0;
  globalMatrixMemSize += yateto::alignedUpper(lina::init::kDivM::size(), alignedReals);
  globalMatrixMemSize += yateto::alignedUpper(lina::init::kDivMT::size(), alignedReals);
  globalMatrixMemSize += yateto::alignedUpper(lina::init::kTDivM::size(), alignedReals);
  globalMatrixMemSize += yateto::alignedUpper(lina::init::kTDivMT::size(), alignedReals);
  globalMatrixMemSize += yateto::computeFamilySize<lina::init::FDivM>(alignedReals);
  globalMatrixMemSize += yateto::computeFamilySize<lina::init::FDivMT>(alignedReals);
  globalMatrixMemSize += yateto::alignedUpper(lina::init::quadrature::size(), alignedReals);

  posix_memalign(reinterpret_cast<void**>(&m_matrixMem), ALIGNMENT, globalMatrixMemSize * sizeof(double));

  double* mem = m_matrixMem;
  yateto::copyTensorToMemAndSetPtr<lina::init::kDivM,   double>(mem, kDivM, ALIGNMENT);
  yateto::copyTensorToMemAndSetPtr<lina::init::kDivMT,  double>(mem, kDivMT, ALIGNMENT);
  yateto::copyTensorToMemAndSetPtr<lina::init::kTDivM,  double>(mem, kTDivM, ALIGNMENT);
  yateto::copyTensorToMemAndSetPtr<lina::init::kTDivMT, double>(mem, kTDivMT, ALIGNMENT);
  yateto::copyFamilyToMemAndSetPtr<lina::init::FDivM,   double>(mem, FDivM, ALIGNMENT);
  yateto::copyFamilyToMemAndSetPtr<lina::init::FDivMT,  double>(mem, FDivMT, ALIGNMENT);
  yateto::copyTensorToMemAndSetPtr<lina::init::quadrature,   double>(mem, quadrature, ALIGNMENT);

  for (unsigned i = 0; i < lina::init::kTDivM::size(); ++i) {
    kTDivM[i] *= -1.0;
  }
  for (unsigned i = 0; i < lina::init::kTDivMT::size(); ++i) {
    kTDivMT[i] *= -1.0;
  }
}

GlobalMatrices::~GlobalMatrices() {
  free(m_matrixMem);
}
