#include "Model.h"

#include <cstring>
#include <generated_code/kernel.h>

void computeA(Material const& material, real A[lina::tensor::star::size(0)], double scale)
{
  auto Av = lina::init::star::view<0>::create(A);
  Av.setZero();
  Av(1,0) = scale * material.K0;
  Av(0,1) = scale / material.rho0;
}

void computeB(Material const& material, real B[lina::tensor::star::size(1)], double scale)
{
  auto Bv = lina::init::star::view<1>::create(B);
  Bv.setZero();
  Bv(2,0) = scale * material.K0;
  Bv(0,2) = scale / material.rho0;
}

void rotateFluxSolver(  double        nx,
                        double        ny,
                        real const    Apm[lina::tensor::Apm::size()],
                        real          fluxSolver[lina::tensor::fluxSolver::size()],
                        double        scale )
{
  real T[lina::tensor::T::size()] __attribute__((aligned(ALIGNMENT))) = {}; // zero initialisation
  real TT[lina::tensor::TT::size()] __attribute__((aligned(ALIGNMENT))) = {}; // zero initialisation
  
  auto Tv = lina::init::T::view::create(T);
  auto TTv = lina::init::TT::view::create(TT);
  
  Tv(0,0) = 1.0;
  Tv(1,1) = nx;
  Tv(2,1) = ny;
  Tv(1,2) = -ny;
  Tv(2,2) = nx;
  
  TTv(0,0) = 1.0;
  TTv(1,1) = nx;
  TTv(2,1) = -ny;
  TTv(1,2) = ny;
  TTv(2,2) = nx;

  lina::kernel::computeFluxSolver krnl;
  krnl.fluxScale = scale;
  krnl.fluxSolver = fluxSolver;
  krnl.T = T;
  krnl.TT = TT;
  krnl.Apm = Apm;
  krnl.execute();
}

void computeAplus( Material const&  local,
                   Material const&  neighbour,
                   real             Aplus[lina::tensor::Apm::size()] )
{
  auto Ap = lina::init::Apm::view::create(Aplus);
  Ap.setZero();
  
  double cm = local.wavespeed();
  double cp = neighbour.wavespeed();
  double div1 = 1.0 / (local.K0 * cp + neighbour.K0 * cm);
  double div2 = div1 / local.rho0;
  Ap(0,0) = local.K0 * cm * cp * div1;
  Ap(1,0) = local.K0 * local.K0 * cp * div1;
  Ap(0,1) = neighbour.K0 * cm * div2;
  Ap(1,1) = local.K0 * neighbour.K0 * div2;
}

void computeAminus( Material const& local,
                    Material const& neighbour,
                    real            Aminus[lina::tensor::Apm::size()] )
{  
  auto Am = lina::init::Apm::view::create(Aminus);
  Am.setZero();
  
  double cm = local.wavespeed();
  double cp = neighbour.wavespeed();
  double div1 = 1.0 / (local.K0 * cp + neighbour.K0 * cm);
  double div2 = div1 / local.rho0;

  Am(0,0) = -local.K0 * cm * cp * div1;
  Am(1,0) = local.K0 * neighbour.K0 * cm * div1;
  Am(0,1) = local.K0 * cp * div2;
  Am(1,1) = -local.K0 * neighbour.K0 * div2;
}
