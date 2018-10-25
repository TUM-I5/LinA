#include "Simulator.h"

#include <algorithm>
#include <cmath>
#include <iostream>

#include "constants.h"
#include "Kernels.h"
#include "Model.h"
#include "Stopwatch.h"

double determineTimestep(double hx, double hy, Grid<Material>& materialGrid)
{
  double maxWaveSpeed = 0.0;
  for (int y = 0; y < materialGrid.Y(); ++y) {
    for (int x = 0; x < materialGrid.X(); ++x) {
      maxWaveSpeed = std::max(maxWaveSpeed, materialGrid.get(x, y).wavespeed());
    }
  }
  
  double PNPM[10]  = {1.0, 0.33, 0.17, 0.1, 0.069, 0.045, 0.038, 0.03, 0.02, 0.015};
  double factor = (CONVERGENCE_ORDER < 10) ? PNPM[CONVERGENCE_ORDER] : 0.25/(2*CONVERGENCE_ORDER-1);

  return factor * hx * hy / (maxWaveSpeed * (hx + hy));
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

  Grid<DegreesOfFreedom> timeIntegratedGrid(globals.X, globals.Y);

  Grid<LocalMatrices> localMatricesGrid(globals.X, globals.Y);
  #pragma omp parallel for collapse(2)
  for (int y = 0; y < globals.Y; ++y) {
    for (int x = 0; x < globals.X; ++x) {
      Material& material = materialGrid.get(x, y);
      LocalMatrices& localMatrices = localMatricesGrid.get(x, y);
      computeA(material, localMatrices.Astar, 1.0 / globals.hx);
      computeB(material, localMatrices.Bstar, 1.0 / globals.hx);

      double Apm[lina::tensor::Apm::size()];
      double fluxScale = -globals.hy / (globals.hx * globals.hy);
      for (unsigned dim = 0; dim < 2; ++dim) {
        for (unsigned side1 = 0; side1 < 2; ++side1) {
          unsigned xn = x + (1-dim)*(2*side1-1);
          unsigned yn = y +    dim *(2*side1-1);
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
    #pragma omp parallel for collapse(2)
    for (int y = 0; y < globals.Y; ++y) {
      for (int x = 0; x < globals.X; ++x) {        
        LocalMatrices& localMatrices = localMatricesGrid.get(x, y);
        DegreesOfFreedom& degreesOfFreedom = degreesOfFreedomGrid.get(x, y);
        DegreesOfFreedom& timeIntegrated = timeIntegratedGrid.get(x, y);
        
        computeAder(timestep, globalMatrices, localMatrices, degreesOfFreedom, timeIntegrated);
        
        computeVolumeIntegral(globalMatrices, localMatrices, timeIntegrated, degreesOfFreedom);

        computeLocalFlux(globalMatrices, localMatrices, timeIntegrated, degreesOfFreedom);
      }
    }

    #pragma omp parallel for collapse(2)
    for (int y = 0; y < globals.Y; ++y) {
      for (int x = 0; x < globals.X; ++x) {
        LocalMatrices& localMatrices = localMatricesGrid.get(x, y);
        DegreesOfFreedom& degreesOfFreedom = degreesOfFreedomGrid.get(x, y);
        
        double* timeIntegrated[2][2];
        for (unsigned dim = 0; dim < 2; ++dim) {
          for (unsigned side1 = 0; side1 < 2; ++side1) {
            unsigned xn = x + (1-dim)*(2*side1-1);
            unsigned yn = y +    dim *(2*side1-1);
            timeIntegrated[dim][side1] = timeIntegratedGrid.get(xn, yn);
          }
        }

        computeNeighbourFlux(globalMatrices, localMatrices, timeIntegrated, degreesOfFreedom);
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
  totalNonZeroFlops *= step * globals.X * globals.Y;
  uint64_t totalHardwareFlops = hardwareFlops;
  totalHardwareFlops *= step * globals.X * globals.Y;

  std::cout << "Time (s): " << wallTime << std::endl;
  std::cout << "Performance (NZ-GFLOPS): " << totalNonZeroFlops / wallTime * 1.0e-9 << std::endl;
  std::cout << "Performance (HW-GFLOPS): " << totalHardwareFlops / wallTime * 1.0e-9 << std::endl;
  
  return step;
}
