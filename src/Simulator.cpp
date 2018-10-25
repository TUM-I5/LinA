#include "Simulator.h"

#include <algorithm>
#include <cmath>
#include <iostream>

#include "constants.h"
#include "Kernels.h"
#include "Model.h"
#include "DGMatrices.h"

double determineTimestep(double hx, double hy, Grid<Material>& materialGrid)
{
  double maxWaveSpeed = 0.0;
  for (int y = 0; y < materialGrid.Y(); ++y) {
    for (int x = 0; x < materialGrid.X(); ++x) {
      maxWaveSpeed = std::max(maxWaveSpeed, materialGrid.get(x, y).wavespeed());
    }
  }
  
  return 0.25 * std::min(hx, hy)/((2*CONVERGENCE_ORDER-1) * maxWaveSpeed);
}

int simulate( GlobalConstants const&  globals,
              Grid<Material>&         materialGrid,
              Grid<DegreesOfFreedom>& degreesOfFreedomGrid,
              WaveFieldWriter&        waveFieldWriter,
              SourceTerm&             sourceterm  )
{
  Grid<DegreesOfFreedom> timeIntegratedGrid(globals.X, globals.Y);
  
  GlobalMatrices globalMatrices;
  Grid<LocalMatrices> localMatricesGrid(globals.X, globals.Y);
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
  
  double time;
  int step = 0;
  for (time = 0.0; time < globals.endTime; time += globals.maxTimestep) {
    waveFieldWriter.writeTimestep(time, degreesOfFreedomGrid);
  
    double timestep = std::min(globals.maxTimestep, globals.endTime - time);
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
  
  waveFieldWriter.writeTimestep(globals.endTime, degreesOfFreedomGrid, true);
  
  return step;
}
