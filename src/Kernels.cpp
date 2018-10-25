#include "Kernels.h"
#include "GEMM.h"
#include "Model.h"
#include <generated_code/init.h>
#include <generated_code/kernel.h>

using namespace lina;

void computeAder( double                  timestep,
                  GlobalMatrices const&   globalMatrices,
                  LocalMatrices const&    localMatrices,
                  DegreesOfFreedom const& degreesOfFreedom,
                  DegreesOfFreedom&       timeIntegrated )
{
  double derivativesBuffer[yateto::computeFamilySize<tensor::dQ>()] __attribute__((aligned(ALIGNMENT)));

  kernel::derivative krnl;
  krnl.kTDivM = globalMatrices.kTDivM;
  krnl.kTDivMT = globalMatrices.kTDivMT;
  krnl.star(0) = localMatrices.Astar;
  krnl.star(1) = localMatrices.Bstar;

  kernel::derivativeTaylorExpansion intKrnl;
  intKrnl.I = timeIntegrated;

  krnl.dQ(0) = const_cast<double*>(degreesOfFreedom);
  intKrnl.dQ(0) = degreesOfFreedom;
  unsigned offset = 0;
  for (unsigned i = 1; i < yateto::numFamilyMembers<tensor::dQ>(); ++i) {
    krnl.dQ(i) = derivativesBuffer + offset;
    intKrnl.dQ(i) = derivativesBuffer + offset;
    offset += tensor::dQ::size(i);
  }
  
  intKrnl.power = timestep;
  intKrnl.execute0();
  
  for (unsigned der = 1; der < CONVERGENCE_ORDER; ++der) {  
    krnl.execute(der);

    intKrnl.power *= timestep / (der+1);
    intKrnl.execute(der);
  }
}

void flopsAder( unsigned int        &nonZeroFlops,
                unsigned int        &hardwareFlops ) {
  // initialization
  nonZeroFlops  += kernel::derivativeTaylorExpansion::nonZeroFlops(0);
  hardwareFlops += kernel::derivativeTaylorExpansion::hardwareFlops(0);

  // interate over derivatives
  for( unsigned der = 1; der < CONVERGENCE_ORDER; ++der ) {
    nonZeroFlops  += kernel::derivative::nonZeroFlops(der);
    hardwareFlops += kernel::derivative::hardwareFlops(der);

    // update of time integrated DOFs
    nonZeroFlops  += kernel::derivativeTaylorExpansion::nonZeroFlops(der);
    hardwareFlops += kernel::derivativeTaylorExpansion::hardwareFlops(der);
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
                        DegreesOfFreedom const& timeIntegrated,
                        DegreesOfFreedom        degreesOfFreedom )
{
  kernel::flux krnl;
  krnl.FDivM = globalMatrices.FDivM;
  krnl.FDivMT = globalMatrices.FDivMT;
  krnl.I = timeIntegrated;
  krnl.Q = degreesOfFreedom;
  krnl._prefetch.I = timeIntegrated + tensor::I::size();
  krnl._prefetch.Q = degreesOfFreedom + tensor::Q::size();
  
  for (unsigned dim = 0; dim < 2; ++dim) {
    for (unsigned side1 = 0; side1 < 2; ++side1) {
      krnl.fluxSolver = localMatrices.fluxSolver[dim][side1][side1];
      krnl.execute(dim, side1, side1);
    }
  }
}

void flopsLocalFlux( unsigned int        &nonZeroFlops,
                     unsigned int        &hardwareFlops ) {
  for (unsigned dim = 0; dim < 2; ++dim) {
    for (unsigned side1 = 0; side1 < 2; ++side1) {
      nonZeroFlops  += kernel::flux::nonZeroFlops(dim, side1, side1);
      hardwareFlops += kernel::flux::hardwareFlops(dim, side1, side1);
    }
  }
}


void computeNeighbourFlux(  GlobalMatrices const&   globalMatrices,
                            LocalMatrices const&    localMatrices,
                            double*                 timeIntegrated[2][2],
                            DegreesOfFreedom        degreesOfFreedom )
{
  kernel::flux krnl;
  krnl.FDivM = globalMatrices.FDivM;
  krnl.FDivMT = globalMatrices.FDivMT;
  krnl.Q = degreesOfFreedom;
  
  for (unsigned dim = 0; dim < 2; ++dim) {
    for (unsigned side1 = 0; side1 < 2; ++side1) {
      unsigned side2 = 1-side1;
      krnl.I = timeIntegrated[dim][side1];
      if (dim != 1 || side1 != 1) {
        krnl._prefetch.I = timeIntegrated[dim+(side1+1)/2][(side1+1)%2];
      }
      krnl.fluxSolver = localMatrices.fluxSolver[dim][side1][side2];
      krnl.execute(dim, side1, side2);
    }
  }
}

void flopsNeighbourFlux( unsigned int        &nonZeroFlops,
                         unsigned int        &hardwareFlops ) {
  for (unsigned dim = 0; dim < 2; ++dim) {
    for (unsigned side1 = 0; side1 < 2; ++side1) {
      unsigned side2 = 1-side1;
      nonZeroFlops  += kernel::flux::nonZeroFlops(dim, side1, side2);
      hardwareFlops += kernel::flux::hardwareFlops(dim, side1, side2);
    }
  }
}

