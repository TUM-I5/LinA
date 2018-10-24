#include "Kernels.h"
#include "GEMM.h"
#include "Model.h"
#include "generated_code/init.h"
#include "generated_code/kernel.h"

using namespace lina;

void computeAder( double                  timestep,
                  GlobalConstants const&  globals,
                  Material const&         material,
                  DegreesOfFreedom const& degreesOfFreedom,
                  DegreesOfFreedom&       timeIntegrated )
{
  double derivativesBuffer[yateto::computeFamilySize<tensor::dQ>()] __attribute__((aligned(ALIGNMENT)));
  double A[NUMBER_OF_QUANTITIES * NUMBER_OF_QUANTITIES];
  double B[NUMBER_OF_QUANTITIES * NUMBER_OF_QUANTITIES];

  computeA(material, A);
  computeB(material, B);
  
  for (unsigned i = 0; i < NUMBER_OF_QUANTITIES*NUMBER_OF_QUANTITIES; ++i) {
    A[i] *= -1.0 / globals.hx;
    B[i] *= -1.0 / globals.hy;
  }

  kernel::derivative krnl;
  krnl.kDivMT = init::kDivMT::Values;
  krnl.star(0) = A;
  krnl.star(1) = B;

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

void computeVolumeIntegral( GlobalConstants const&  globals,
                            Material const&         material,
                            DegreesOfFreedom const& timeIntegrated,
                            DegreesOfFreedom&       degreesOfFreedom )
{
  double A[NUMBER_OF_QUANTITIES * NUMBER_OF_QUANTITIES];
  double B[NUMBER_OF_QUANTITIES * NUMBER_OF_QUANTITIES];
  
  computeA(material, A);
  computeB(material, B);
  
  for (unsigned i = 0; i < NUMBER_OF_QUANTITIES*NUMBER_OF_QUANTITIES; ++i) {
    A[i] *= 1.0 / globals.hx;
    B[i] *= 1.0 / globals.hy;
  }

  kernel::volume krnl;
  krnl.I = timeIntegrated;
  krnl.Q = degreesOfFreedom;
  krnl.kDivM = lina::init::kDivM::Values;
  krnl.star(0) = A;
  krnl.star(1) = B;

  krnl.execute();
}

void computeFlux( double                  factor,
                  unsigned                dim,
                  unsigned                side1,
                  unsigned                side2,
                  double                  rotatedFluxSolver[NUMBER_OF_QUANTITIES*NUMBER_OF_QUANTITIES],
                  DegreesOfFreedom const& timeIntegrated,
                  DegreesOfFreedom        degreesOfFreedom )
{
  for (unsigned i = 0; i < NUMBER_OF_QUANTITIES*NUMBER_OF_QUANTITIES; ++i) {
    rotatedFluxSolver[i] *= factor;
  }
  
  kernel::flux krnl;
  krnl.FDivM(side1,side2) = init::FDivM::Values[side2*2+side1];
  krnl.I = timeIntegrated;
  krnl.Q = degreesOfFreedom;
  krnl.fluxSolver = rotatedFluxSolver;
  
  krnl.execute(dim, side1, side2);
}
