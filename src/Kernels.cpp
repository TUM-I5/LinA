#include "Kernels.h"
#include "GEMM.h"
#include "Model.h"
#include <generated_code/init.h>
#include <generated_code/kernel.h>

using namespace lina;

void computeAder( real                  timestep,
                  GlobalMatrices const&   globalMatrices,
                  LocalMatrices const&    localMatrices,
                  DegreesOfFreedom const& degreesOfFreedom,
                  DegreesOfFreedom&       timeIntegrated )
{
  real dQcur[tensor::dQcur::size()] __attribute__((aligned(ALIGNMENT)));
  real dQnext[tensor::dQnext::size()] __attribute__((aligned(ALIGNMENT)));

  kernel::derivative krnl;
  krnl.kTDivM = globalMatrices.kTDivM;
  krnl.kTDivMT = globalMatrices.kTDivMT;
  krnl.star(0) = localMatrices.Astar;
  krnl.star(1) = localMatrices.Bstar;
  krnl.dQcur = degreesOfFreedom;
  krnl.dQnext = dQnext;

  kernel::derivativeTaylorExpansion intKrnl;
  intKrnl.I = timeIntegrated;
  intKrnl.Q = degreesOfFreedom;

  intKrnl.power = timestep;
  intKrnl.execute0();

  krnl.execute();
  krnl.dQcur = dQcur;
  intKrnl.dQnext = krnl.dQnext;
  intKrnl.power *= timestep / 2;
  intKrnl.execute1();

  for (unsigned der = 2; der < CONVERGENCE_ORDER; ++der) {
    real const* tmp = krnl.dQnext;
    krnl.dQnext = const_cast<real*>(krnl.dQcur);
    krnl.dQcur = tmp;

    krnl.execute();

    intKrnl.power *= timestep / (der+1);
    intKrnl.dQnext = krnl.dQnext;
    intKrnl.execute1();
  }
}

void flopsAder( unsigned int        &nonZeroFlops,
                unsigned int        &hardwareFlops ) {
  // initialization
  nonZeroFlops  += kernel::derivativeTaylorExpansion::nonZeroFlops(0);
  hardwareFlops += kernel::derivativeTaylorExpansion::hardwareFlops(0);

  // interate over derivatives
  for( unsigned der = 1; der < CONVERGENCE_ORDER; ++der ) {
    nonZeroFlops  += kernel::derivative::NonZeroFlops;
    hardwareFlops += kernel::derivative::HardwareFlops;

    // update of time integrated DOFs
    nonZeroFlops  += kernel::derivativeTaylorExpansion::nonZeroFlops(1);
    hardwareFlops += kernel::derivativeTaylorExpansion::hardwareFlops(1);
  }
}


void computeVolumeIntegral( GlobalMatrices const&   globalMatrices,
                            LocalMatrices const&    localMatrices,
                            DegreesOfFreedom const& timeIntegrated,
                            DegreesOfFreedom&       degreesOfFreedom )
{
  kernel::volume krnl;
  krnl.I = timeIntegrated;
  krnl.Q = degreesOfFreedom;
  krnl.kDivM = globalMatrices.kDivM;
  krnl.kDivMT = globalMatrices.kDivMT;
  krnl.star(0) = localMatrices.Astar;
  krnl.star(1) = localMatrices.Bstar;

  krnl.execute();
}

void flopsVolume( unsigned int        &nonZeroFlops,
                  unsigned int        &hardwareFlops ) {
  nonZeroFlops  += kernel::volume::NonZeroFlops;
  hardwareFlops += kernel::volume::HardwareFlops;
}

void computeLocalFlux(  GlobalMatrices const&   globalMatrices,
                        LocalMatrices const&    localMatrices,
                        DegreesOfFreedom        timeIntegrated,
                        real*                   timeIntegratedEdge[2][2],
                        DegreesOfFreedom        degreesOfFreedom )
{
  kernel::flux krnl;
  krnl.FDivM = globalMatrices.FDivM;
  krnl.Q = degreesOfFreedom;
  
  kernel::evaluateSide evalKrnl;
  evalKrnl.F = globalMatrices.F;
  evalKrnl.I = timeIntegrated;
  
  for (unsigned dim = 0; dim < 2; ++dim) {
    for (unsigned side = 0; side < 2; ++side) {
      evalKrnl.Q1 = timeIntegratedEdge[dim][side];
      evalKrnl.execute(dim, side);
      krnl.Q1 = timeIntegratedEdge[dim][side];
      krnl.fluxSolver = localMatrices.fluxSolver[dim][side][side];
      krnl.execute(dim, side);
    }
  }
}

void flopsLocalFlux( unsigned int        &nonZeroFlops,
                     unsigned int        &hardwareFlops ) {
  for (unsigned dim = 0; dim < 2; ++dim) {
    for (unsigned side = 0; side < 2; ++side) {
      nonZeroFlops  += kernel::evaluateSide::nonZeroFlops(dim, side);
      hardwareFlops += kernel::evaluateSide::hardwareFlops(dim, side);
      nonZeroFlops  += kernel::flux::nonZeroFlops(dim, side);
      hardwareFlops += kernel::flux::hardwareFlops(dim, side);
    }
  }
}


void computeNeighbourFlux(  GlobalMatrices const&   globalMatrices,
                            LocalMatrices const&    localMatrices,
                            real*                   timeIntegratedEdge[2][2],
                            DegreesOfFreedom        degreesOfFreedom )
{
  kernel::flux krnl;
  krnl.FDivM = globalMatrices.FDivM;
  krnl.Q = degreesOfFreedom;
  
  for (unsigned dim = 0; dim < 2; ++dim) {
    for (unsigned side1 = 0; side1 < 2; ++side1) {
      unsigned side2 = 1-side1;
      krnl.Q1 = timeIntegratedEdge[dim][side1];
      krnl.fluxSolver = localMatrices.fluxSolver[dim][side1][side2];
      krnl.execute(dim, side1);
    }
  }
}

void flopsNeighbourFlux( unsigned int        &nonZeroFlops,
                         unsigned int        &hardwareFlops ) {
  for (unsigned dim = 0; dim < 2; ++dim) {
    for (unsigned side1 = 0; side1 < 2; ++side1) {
      nonZeroFlops  += kernel::flux::nonZeroFlops(dim, side1);
      hardwareFlops += kernel::flux::hardwareFlops(dim, side1);
    }
  }
}

