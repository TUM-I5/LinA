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
  krnl.kDivMT = globalMatrices.kDivMT;
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

void computeVolumeIntegral( GlobalMatrices const&   globalMatrices,
                            LocalMatrices const&    localMatrices,
                            DegreesOfFreedom const& timeIntegrated,
                            DegreesOfFreedom&       degreesOfFreedom )
{
  kernel::volume krnl;
  krnl.I = timeIntegrated;
  krnl.Q = degreesOfFreedom;
  krnl.kDivM = globalMatrices.kDivM;
  krnl.star(0) = localMatrices.Astar;
  krnl.star(1) = localMatrices.Bstar;

  krnl.execute();
}

void computeLocalFlux(  GlobalMatrices const&   globalMatrices,
                        LocalMatrices const&    localMatrices,
                        DegreesOfFreedom const& timeIntegrated,
                        DegreesOfFreedom        degreesOfFreedom )
{
  kernel::flux krnl;
  krnl.FDivM = globalMatrices.FDivM;
  krnl.I = timeIntegrated;
  krnl.Q = degreesOfFreedom;
  
  for (unsigned dim = 0; dim < 2; ++dim) {
    for (unsigned side1 = 0; side1 < 2; ++side1) {
      krnl.fluxSolver = localMatrices.fluxSolver[dim][side1][side1];
      krnl.execute(dim, side1, side1);
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
  krnl.Q = degreesOfFreedom;
  
  for (unsigned dim = 0; dim < 2; ++dim) {
    for (unsigned side1 = 0; side1 < 2; ++side1) {
      unsigned side2 = 1-side1;
      krnl.I = timeIntegrated[dim][side1];
      krnl.fluxSolver = localMatrices.fluxSolver[dim][side1][side2];
      krnl.execute(dim, side1, side2);
    }
  }
}
