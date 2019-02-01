#ifndef MODEL_H_
#define MODEL_H_

#include "typedefs.h"
#include <generated_code/init.h>

/** Returns A in column-major storage */
void computeA(Material const& material, real A[lina::tensor::star::size(0)], double scale);

/** Returns B in column-major storage */
void computeB(Material const& material, real B[lina::tensor::star::size(1)], double scale);

/** Returns rotated flux solver in column-major storage for face-aligned coordinate system in direction (nx, ny).
  * (nx,ny) must be a unit vector, i.e. nx^2 + ny^2 = 1. */
void rotateFluxSolver(  double        nx,
                        double        ny,
                        real   const  Apm[lina::tensor::Apm::size()],
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
