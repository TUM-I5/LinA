#include "Simulator.h"

#include <algorithm>
#include <cmath>
#include <iostream>

#include "constants.h"
#include "Kernels.h"
#include "Model.h"
#include "Stopwatch.h"
#include "generated_code/kernel.h" 

double determineTimestep(double hx, double hy, double cfl, Grid<Material>& materialGrid)
{
  double maxWaveSpeed = 0.0;
  for (int y = 0; y < materialGrid.Y(); ++y) {
    for (int x = 0; x < materialGrid.X(); ++x) {
      maxWaveSpeed = std::max(maxWaveSpeed, materialGrid.get(x, y).wavespeed());
    }
  }

  double PNPM[10]  = {1.0, 0.33, 0.17, 0.1, 0.069, 0.045, 0.038, 0.03, 0.02, 0.015};
  double factor = (CONVERGENCE_ORDER < 10) ? PNPM[CONVERGENCE_ORDER] : 1.0/(2*CONVERGENCE_ORDER-1);

  if (CONVERGENCE_ORDER >= 10 && cfl == 1.0) {
    std::cerr << "Warning: Stability factor unknown for order " << CONVERGENCE_ORDER << ". Defaulting to guess 1/(2N+1). Please adjust the CFL number with option -c." << std::endl;
  }

  return cfl * factor * hx * hy / (maxWaveSpeed * (hx + hy));
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

  Grid<EdgeDOFs>* timeIntegratedEdges[2][2];
  for (int dim = 0; dim < 2; ++dim) {
    for (int side = 0; side < 2; ++side) {
      timeIntegratedEdges[dim][side] = new Grid<EdgeDOFs>(globals.X, globals.Y);
    }
  }

  Grid<LocalMatrices> localMatricesGrid(globals.X, globals.Y);
  #pragma omp parallel for collapse(2)
  for (int y = 0; y < globals.Y; ++y) {
    for (int x = 0; x < globals.X; ++x) {
      Material& material = materialGrid.get(x, y);
      LocalMatrices& localMatrices = localMatricesGrid.get(x, y);
      computeA(material, localMatrices.Astar, 1.0 / globals.hx);
      computeB(material, localMatrices.Bstar, 1.0 / globals.hx);

      real Apm[lina::tensor::Apm::size()];
      real fluxScale = -globals.hy / (globals.hx * globals.hy);
      for (int dim = 0; dim < 2; ++dim) {
        for (int side1 = 0; side1 < 2; ++side1) {
          int xn = x + (1-dim)*(2*side1-1);
          int yn = y +    dim *(2*side1-1);
          double nx = (1-dim)*(2.0*side1-1.0);
          double ny =    dim *(2.0*side1-1.0);
          for (unsigned side2 = 0; side2 < 2; ++side2) {
            if (side1 != side2) {
              computeAminus(material, materialGrid.get(xn, yn), Apm);
            } else {
              computeAplus(material, materialGrid.get(xn, yn), Apm);
            }
            rotateFluxSolver(nx, ny, Apm, localMatrices.fluxSolver[dim][side1][side2], fluxScale);
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
      #pragma omp for collapse(2)
      for (int y = 0; y < globals.Y; ++y) {
        for (int x = 0; x < globals.X; ++x) {
          LocalMatrices& localMatrices = localMatricesGrid.get(x, y);
          DegreesOfFreedom& degreesOfFreedom = degreesOfFreedomGrid.get(x, y);

          computeAder(timestep, globalMatrices, localMatrices, degreesOfFreedom, timeIntegrated);

          computeVolumeIntegral(globalMatrices, localMatrices, timeIntegrated, degreesOfFreedom);

          real* timeIntegratedEdge[2][2];
          for (int dim = 0; dim < 2; ++dim) {
            for (int side = 0; side < 2; ++side) {
              timeIntegratedEdge[dim][side] = timeIntegratedEdges[dim][side]->get(x, y);
            }
          }
          computeLocalFlux(globalMatrices, localMatrices, timeIntegrated, timeIntegratedEdge);
        }
      }
    }

    #pragma omp parallel for collapse(2)
    for (int y = 0; y < globals.Y; ++y) {
      for (int x = 0; x < globals.X; ++x) {
        LocalMatrices& localMatrices = localMatricesGrid.get(x, y);
        DegreesOfFreedom& degreesOfFreedom = degreesOfFreedomGrid.get(x, y);
        
        real* timeIntegratedEdge[2][2][2];
        for (int dim = 0; dim < 2; ++dim) {
          for (int side = 0; side < 2; ++side) {
            int xn = x + (1-dim)*(2*side-1);
            int yn = y +    dim *(2*side-1);
            timeIntegratedEdge[dim][side][0] = timeIntegratedEdges[dim][side]->get(x, y);
            timeIntegratedEdge[dim][side][1] = timeIntegratedEdges[dim][1-side]->get(xn, yn);
          }
        }

        computeNeighbourFlux(globalMatrices, localMatrices, timeIntegratedEdge, degreesOfFreedom);
      }
    }
    
    if (sourceterm.x >= 0 && sourceterm.x < globals.X && sourceterm.y >= 0 && sourceterm.y < globals.Y) {
      double areaInv = 1. / (globals.hx*globals.hy);
      DegreesOfFreedom& degreesOfFreedom = degreesOfFreedomGrid.get(sourceterm.x, sourceterm.y);
      double timeIntegral = (*sourceterm.antiderivative)(time + timestep) - (*sourceterm.antiderivative)(time);
      real source[NUMBER_OF_QUANTITIES] __attribute__((aligned(ALIGNMENT))) = {};
      source[sourceterm.quantity] = areaInv * timeIntegral;
      lina::kernel::applySource krnl;
      krnl.phi(0) = sourceterm.phiDivM[0];
      krnl.phi(1) = sourceterm.phiDivM[1];
      krnl.source = source;
      krnl.Q = degreesOfFreedom;
      krnl.execute();
    }
    
    ++step;
    if (step % 100 == 0) {
      std::cout << "At time / timestep: " << time << " / " << step << std::endl;
    }
  }

  auto wallTime = sw.stop();
  
  waveFieldWriter.writeTimestep(globals.endTime, degreesOfFreedomGrid, true);

  uint64_t totalNonZeroFlops = nonZeroFlops;
  totalNonZeroFlops *= step * globals.X * globals.Y;
  uint64_t totalHardwareFlops = hardwareFlops;
  totalHardwareFlops *= step * globals.X * globals.Y;

  std::cout << "Time (s): " << wallTime << std::endl;
  std::cout << "Performance (NZ-GFLOPS): " << totalNonZeroFlops / wallTime * 1.0e-9 << std::endl;
  std::cout << "Performance (HW-GFLOPS): " << totalHardwareFlops / wallTime * 1.0e-9 << std::endl;
  
  for (int dim = 0; dim < 2; ++dim) {
    for (int side = 0; side < 2; ++side) {
      delete timeIntegratedEdges[dim][side];
    }
  }

  return step;
}
