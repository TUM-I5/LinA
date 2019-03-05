#include "DGMatrices.h"
#include <yateto.h>
#include <generated_code/init.h>
#include <iostream>

GlobalMatrices::GlobalMatrices() {
  size_t alignedReals = yateto::alignedReals<real>(ALIGNMENT);
  unsigned globalMatrixMemSize = 0;
  globalMatrixMemSize += yateto::alignedUpper(lina::init::kDivM::size(), alignedReals);
  globalMatrixMemSize += yateto::alignedUpper(lina::init::kDivMT::size(), alignedReals);
  globalMatrixMemSize += yateto::alignedUpper(lina::init::kTDivM::size(), alignedReals);
  globalMatrixMemSize += yateto::alignedUpper(lina::init::kTDivMT::size(), alignedReals);
  globalMatrixMemSize += yateto::computeFamilySize<lina::init::FDivM>(alignedReals);
  globalMatrixMemSize += yateto::computeFamilySize<lina::init::F>(alignedReals);
  globalMatrixMemSize += yateto::alignedUpper(lina::init::quadrature::size(), alignedReals);
  globalMatrixMemSize += yateto::alignedUpper(lina::init::mInv::size(), alignedReals);

  int err = posix_memalign(reinterpret_cast<void**>(&m_matrixMem), ALIGNMENT, globalMatrixMemSize * sizeof(real));
  if (err) {
    std::cerr << "Failed to allocate " << globalMatrixMemSize * sizeof(real) << " bytes in " << __FILE__ << std::endl;
    exit(EXIT_FAILURE);
  }

  real* mem = m_matrixMem;
  yateto::copyTensorToMemAndSetPtr<lina::init::kDivM,   real>(mem, kDivM, ALIGNMENT);
  yateto::copyTensorToMemAndSetPtr<lina::init::kDivMT,  real>(mem, kDivMT, ALIGNMENT);
  yateto::copyTensorToMemAndSetPtr<lina::init::kTDivM,  real>(mem, kTDivM, ALIGNMENT);
  yateto::copyTensorToMemAndSetPtr<lina::init::kTDivMT, real>(mem, kTDivMT, ALIGNMENT);
  yateto::copyFamilyToMemAndSetPtr<lina::init::FDivM,   real>(mem, FDivM, ALIGNMENT);
  yateto::copyFamilyToMemAndSetPtr<lina::init::F,       real>(mem, F, ALIGNMENT);
  yateto::copyTensorToMemAndSetPtr<lina::init::quadrature,   real>(mem, quadrature, ALIGNMENT);
  yateto::copyTensorToMemAndSetPtr<lina::init::mInv,    real>(mem, mInv, ALIGNMENT);

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
