#ifndef MODEL_H_
#define MODEL_H_

#include "typedefs.h"
#include <generated_code/init.h>

template<int Dim>
void computeJacobian(Material const& material, real A[lina::tensor::star::size(Dim)], double scale)
{
  auto Av = lina::init::star::view<Dim>::create(A);
  Av.setZero();
  Av(1+Dim,0) = scale * material.K0;
  Av(0,1+Dim) = scale / material.rho0;
}

/** Returns rotated flux solver in column-major storage for face-aligned coordinate system in direction (nx, ny, nz).
  * (nx,ny,nz) must be a unit vector, i.e. nx^2 + ny^2 + nz^2 = 1. */
void rotateFluxSolver(  double        nx,
                        double        ny,
                        double        nz,
                        real const    Apm[lina::tensor::Apm::size()],
                        real          fluxSolver[lina::tensor::fluxSolver::size()],
                        double        scale );

/** Returns A^+ in column-major storage */
void computeAplus( Material const&  local,
                   Material const&  neighbour,
                   real             Aplus[lina::tensor::Apm::size()] );

/** Returns A^- in column-major storage */
void computeAminus( Material const& local,
                    Material const& neighbour,
                    real            Aminus[lina::tensor::Apm::size()] );

#endif // MODEL_H_
