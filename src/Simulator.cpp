#include "Simulator.h"

#include <algorithm>
#include <cmath>
#include <iostream>

#include "constants.h"
#include "Kernels.h"
#include "Model.h"
#include "Stopwatch.h"

double determineTimestep(double hx, double hy, double hz, double cfl, Grid<Material>& materialGrid)
{
  double maxWaveSpeed = 0.0;
  for (int z = 0; z < materialGrid.Z(); ++z) {
    for (int y = 0; y < materialGrid.Y(); ++y) {
      for (int x = 0; x < materialGrid.X(); ++x) {
        maxWaveSpeed = std::max(maxWaveSpeed, materialGrid.get(x, y, z).wavespeed());
      }
    }
  }

  double PNPM[10]  = {1.0, 0.33, 0.17, 0.1, 0.069, 0.045, 0.038, 0.03, 0.02, 0.015};
  double factor = (CONVERGENCE_ORDER < 10) ? PNPM[CONVERGENCE_ORDER] : 1.0/(2*CONVERGENCE_ORDER-1);

  if (CONVERGENCE_ORDER >= 10 && cfl == 1.0) {
    std::cerr << "Warning: Stability factor unknown for order " << CONVERGENCE_ORDER << ". Defaulting to guess 1/(2N+1). Please adjust the CFL number with option -c." << std::endl;
  }

  return cfl * factor / maxWaveSpeed / (1.0/hx + 1.0/hy + 1.0/hz);
}

int simulate( GlobalConstants const&  globals,
              GlobalMatrices const&   globalMatrices,
              Grid<Material>&         materialGrid,
              Grid<DegreesOfFreedom>& degreesOfFreedomGrid,
              WaveFieldWriter&        waveFieldWriter,
              SourceTerm&             sourceterm  )
{
  unsigned nonZeroFlops = 0;
  unsigned hardwareFlops = 0;
  flopsAder(nonZeroFlops, hardwareFlops);
  flopsVolume(nonZeroFlops, hardwareFlops);
  flopsLocalFlux(nonZeroFlops, hardwareFlops);
  flopsNeighbourFlux(nonZeroFlops, hardwareFlops);

  std::cout << "Non-zero flops / element: " << nonZeroFlops << std::endl;
  std::cout << "Hardware flops / element: " << hardwareFlops << std::endl;

  Grid<EdgeDOFs>* timeIntegratedEdges[3][2];
  for (int dim = 0; dim < 3; ++dim) {
    for (int side = 0; side < 2; ++side) {
      timeIntegratedEdges[dim][side] = new Grid<EdgeDOFs>(globals.X, globals.Y, globals.Z);
    }
  }

  Grid<LocalMatrices> localMatricesGrid(globals.X, globals.Y, globals.Z);
  #pragma omp parallel for collapse(3)
  for (int z = 0; z < globals.Z; ++z) {
    for (int y = 0; y < globals.Y; ++y) {
      for (int x = 0; x < globals.X; ++x) {
        Material& material = materialGrid.get(x, y, z);
        LocalMatrices& localMatrices = localMatricesGrid.get(x, y, z);
        computeJacobian<0>(material, localMatrices.Astar, 1.0 / globals.hx);
        computeJacobian<1>(material, localMatrices.Bstar, 1.0 / globals.hy);
        computeJacobian<2>(material, localMatrices.Cstar, 1.0 / globals.hz);

        real Apm[lina::tensor::Apm::size()]; 
        double h[3] = {globals.hx, globals.hy, globals.hz};
        for (int dim = 0; dim < 3; ++dim) {
          for (int side1 = 0; side1 < 2; ++side1) {
            int off = 2*side1-1;
            int xn = x + ((dim == 0) ? off : 0);
            int yn = y + ((dim == 1) ? off : 0);
            int zn = z + ((dim == 2) ? off : 0);
            double nx = (dim == 0) ? (2.0*side1-1.0) : 0.0;
            double ny = (dim == 1) ? (2.0*side1-1.0) : 0.0;
            double nz = (dim == 2) ? (2.0*side1-1.0) : 0.0;
            for (int side2 = 0; side2 < 2; ++side2) {
              if (side1 != side2) {
                computeAminus(material, materialGrid.get(xn, yn, zn), Apm);
              } else {
                computeAplus(material, materialGrid.get(xn, yn, zn), Apm);
              }
              double fluxScale = -1.0 / h[dim];
              rotateFluxSolver(nx, ny, nz, Apm, localMatrices.fluxSolver[dim][side1][side2], fluxScale);
            }
          }
        }
      }
    }
  }

  Stopwatch sw;
  sw.start();
  
  double time;
  int step = 0;
  for (time = 0.0; time < globals.endTime; time += globals.maxTimestep) {
    waveFieldWriter.writeTimestep(time, degreesOfFreedomGrid);

    double timestep = std::min(globals.maxTimestep, globals.endTime - time);
    #pragma omp parallel
    {
      DegreesOfFreedom timeIntegrated __attribute((aligned(ALIGNMENT)));
      #pragma omp for collapse(3)
      for (int z = 0; z < globals.Z; ++z) {
        for (int y = 0; y < globals.Y; ++y) {
          for (int x = 0; x < globals.X; ++x) {
            LocalMatrices& localMatrices = localMatricesGrid.get(x, y, z);
            DegreesOfFreedom& degreesOfFreedom = degreesOfFreedomGrid.get(x, y, z);

            computeAder(timestep, globalMatrices, localMatrices, degreesOfFreedom, timeIntegrated);

            computeVolumeIntegral(globalMatrices, localMatrices, timeIntegrated, degreesOfFreedom);

            real* timeIntegratedEdge[3][2];
            for (int dim = 0; dim < 3; ++dim) {
              for (int side = 0; side < 2; ++side) {
                timeIntegratedEdge[dim][side] = timeIntegratedEdges[dim][side]->get(x, y, z);
              }
            }
            computeLocalFlux(globalMatrices, localMatrices, timeIntegrated, timeIntegratedEdge);
          }
        }
      }
    }

    #pragma omp parallel for collapse(3)
    for (int z = 0; z < globals.Z; ++z) {
      for (int y = 0; y < globals.Y; ++y) {
        for (int x = 0; x < globals.X; ++x) {
          LocalMatrices& localMatrices = localMatricesGrid.get(x, y, z);
          DegreesOfFreedom& degreesOfFreedom = degreesOfFreedomGrid.get(x, y, z);
        
          real* timeIntegratedEdge[3][2][2];
          for (int dim = 0; dim < 3; ++dim) {
            for (int side = 0; side < 2; ++side) {
              int off = 2*side-1;
              int xn = x + ((dim == 0) ? off : 0);
              int yn = y + ((dim == 1) ? off : 0);
              int zn = z + ((dim == 2) ? off : 0);
              timeIntegratedEdge[dim][side][0] = timeIntegratedEdges[dim][side]->get(x, y, z);
              timeIntegratedEdge[dim][side][1] = timeIntegratedEdges[dim][1-side]->get(xn, yn, zn);
            }
          }

          computeNeighbourFlux(globalMatrices, localMatrices, timeIntegratedEdge, degreesOfFreedom);
        }
      }
    }
    
    /*if (sourceterm.x >= 0 && sourceterm.x < globals.X && sourceterm.y >= 0 && sourceterm.y < globals.Y) {
      double areaInv = 1. / (globals.hx*globals.hy);
      DegreesOfFreedom& degreesOfFreedom = degreesOfFreedomGrid.get(sourceterm.x, sourceterm.y);
      double timeIntegral = (*sourceterm.antiderivative)(time + timestep) - (*sourceterm.antiderivative)(time);
      for (unsigned b = 0; b < NUMBER_OF_BASIS_FUNCTIONS; ++b) {
        degreesOfFreedom[sourceterm.quantity * NUMBER_OF_BASIS_FUNCTIONS + b] += areaInv * timeIntegral * sourceterm.phi[b];
      }
    }*/
    
    ++step;
    if (step % 100 == 0) {
      std::cout << "At time / timestep: " << time << " / " << step << std::endl;
    }
  }

  auto wallTime = sw.stop();
  
  waveFieldWriter.writeTimestep(globals.endTime, degreesOfFreedomGrid, true);

  uint64_t totalNonZeroFlops = nonZeroFlops;
  totalNonZeroFlops *= step * globals.X * globals.Y * globals.Z;
  uint64_t totalHardwareFlops = hardwareFlops;
  totalHardwareFlops *= step * globals.X * globals.Y * globals.Z;

  std::cout << "Time (s): " << wallTime << std::endl;
  std::cout << "Performance (NZ-GFLOPS): " << totalNonZeroFlops / wallTime * 1.0e-9 << std::endl;
  std::cout << "Performance (HW-GFLOPS): " << totalHardwareFlops / wallTime * 1.0e-9 << std::endl;
  
  for (int dim = 0; dim < 3; ++dim) {
    for (int side = 0; side < 2; ++side) {
      delete timeIntegratedEdges[dim][side];
    }
  }

  return step;
}
