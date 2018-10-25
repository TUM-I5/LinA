#ifndef KERNELS_H_
#define KERNELS_H_

#include "typedefs.h"
#include "DGMatrices.h"

void computeAder( double                  timestep,
                  GlobalMatrices const&   globalMatrices,
                  LocalMatrices const&    localMatrices,
                  DegreesOfFreedom const& degreesOfFreedom,
                  DegreesOfFreedom&       timeIntegrated );
                  
void computeVolumeIntegral( GlobalMatrices const&   globalMatrices,
                            LocalMatrices const&    localMatrices,
                            DegreesOfFreedom const& timeIntegrated,
                            DegreesOfFreedom&       degreesOfFreedom );

void computeLocalFlux(  GlobalMatrices const&   globalMatrices,
                        LocalMatrices const&    localMatrices,
                        DegreesOfFreedom const& timeIntegrated,
                        DegreesOfFreedom        degreesOfFreedom );

void computeNeighbourFlux(  GlobalMatrices const&   globalMatrices,
                            LocalMatrices const&    localMatrices,
                            double*                 timeIntegrated[2][2],
                            DegreesOfFreedom        degreesOfFreedom );

#endif // KERNELS_H_
