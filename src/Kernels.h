#ifndef KERNELS_H_
#define KERNELS_H_

#include "typedefs.h"
#include "DGMatrices.h"

void computeAder( real                  timestep,
                  GlobalMatrices const&   globalMatrices,
                  LocalMatrices const&    localMatrices,
                  DegreesOfFreedom const& degreesOfFreedom,
                  DegreesOfFreedom&       timeIntegrated );

void flopsAder( unsigned &nonZeroFlops, unsigned &hardwareFlops );
                  
void computeVolumeIntegral( GlobalMatrices const&   globalMatrices,
                            LocalMatrices const&    localMatrices,
                            DegreesOfFreedom const& timeIntegrated,
                            DegreesOfFreedom&       degreesOfFreedom );

void flopsVolume( unsigned &nonZeroFlops, unsigned &hardwareFlops );

void computeLocalFlux(  GlobalMatrices const&   globalMatrices,
                        LocalMatrices const&    localMatrices,
                        DegreesOfFreedom        timeIntegrated,
                        real*                   timeIntegratedEdge[2][2],
                        DegreesOfFreedom        degreesOfFreedom );

void flopsLocalFlux( unsigned &nonZeroFlops, unsigned &hardwareFlops );

void computeNeighbourFlux(  GlobalMatrices const&   globalMatrices,
                            LocalMatrices const&    localMatrices,
                            real*                 timeIntegratedEdge[2][2],
                            DegreesOfFreedom        degreesOfFreedom );

void flopsNeighbourFlux( unsigned &nonZeroFlops, unsigned &hardwareFlops );

#endif // KERNELS_H_
