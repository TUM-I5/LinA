#ifndef SIMULATOR_H_
#define SIMULATOR_H_

#include "typedefs.h"
#include "Grid.h"
#include "WaveFieldWriter.h"
#include "DGMatrices.h"

double determineTimestep(double hx, double hy, double hz, Grid<Material>& materialGrid);

int simulate( GlobalConstants const&  globals,
              GlobalMatrices const&   globalMatrices,
              Grid<Material>&         materialGrid,
              Grid<DegreesOfFreedom>& degreesOfFreedomGrid,
              WaveFieldWriter&        waveFieldWriter,
              SourceTerm&             sourceterm  );

#endif // SIMULATOR_H_
