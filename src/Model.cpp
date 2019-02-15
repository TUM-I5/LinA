#include "Model.h"

#include <cstring>
#include <generated_code/kernel.h>

void rotateFluxSolver(  double        nx,
                        double        ny,
                        double        nz,
                        real const    Apm[lina::tensor::Apm::size()],
                        real          fluxSolver[lina::tensor::fluxSolver::size()],
                        double        scale )
{
  real T[lina::tensor::T::size()] __attribute__((aligned(ALIGNMENT))) = {}; // zero initialisation
  real TT[lina::tensor::TT::size()] __attribute__((aligned(ALIGNMENT))) = {}; // zero initialisation
  
  auto Tv = lina::init::T::view::create(T);
  auto TTv = lina::init::TT::view::create(TT);
  
  double sx = 0.0, sy = 0.0, sz = 0.0, tx = 0.0, ty = 0.0, tz = 0.0;
  if (fabs(nx) > fabs(ny)) {
    sx = nz;
    sz = -nx;
  } else {
    sy = -nz;
    sz = ny;
  }
  double norm = sqrt(sx*sx + sy*sy + sz*sz);
  sx /= norm;
  sy /= norm;
  sz /= norm;
  
  // t = n x s
  tx = ny*sz - nz*sy;
  ty = nz*sx - nx*sz;
  tz = nx*sy - ny*sx;
  
  Tv(0,0) = 1.0;
  Tv(1,1) = nx;
  Tv(2,1) = ny;
  Tv(3,1) = nz;
  Tv(1,2) = sx;
  Tv(2,2) = sy;
  Tv(3,2) = sz;
  Tv(1,3) = tx;
  Tv(2,3) = ty;
  Tv(3,3) = tz;
  
  TTv(0,0) = 1.0;
  TTv(1,1) = nx;
  TTv(1,2) = ny;
  TTv(1,3) = nz;
  TTv(2,1) = sx;
  TTv(2,2) = sy;
  TTv(2,3) = sz;
  TTv(3,1) = tx;
  TTv(3,2) = ty;
  TTv(3,3) = tz;

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
