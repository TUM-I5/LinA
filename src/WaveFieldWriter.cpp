#include "WaveFieldWriter.h"

#include <iostream>
#include <iomanip>
#include <sstream>
#include <cstdio>
#include <limits>

#include "basisfunctions.h"
#include "GEMM.h"
#include "nodes.h"
#include "generated_code/init.h"

WaveFieldWriter::WaveFieldWriter(std::string const& baseName, GlobalConstants const& globals, double interval)
  : m_step(0), m_interval(interval), m_lastTime(-std::numeric_limits<double>::max())
{
  if (!baseName.empty()) {
    std::size_t lastFound = 0;
    std::size_t found;
    while ((found = baseName.find("/", lastFound+1)) != std::string::npos) {
      lastFound = found;
    }
    if (lastFound > 0) {
      ++lastFound;
    }
    m_dirName = baseName.substr(0, lastFound);
    m_baseName = baseName.substr(lastFound);

    std::stringstream xGridName, yGridName;
    xGridName << m_baseName << "_x.bin";
    yGridName << m_baseName << "_y.bin";

    m_xdmf.open((baseName + ".xdmf").c_str());
    m_xdmf  << "<?xml version=\"1.0\" ?>" << std::endl
            << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\">" << std::endl
            << "<Xdmf Version=\"2.0\">" << std::endl
            << "  <Domain>" << std::endl
            << "    <Topology TopologyType=\"3DRectMesh\" Dimensions=\"1 " << NUMBER_OF_BASIS_FUNCTIONS * globals.Y << " " << NUMBER_OF_BASIS_FUNCTIONS * globals.X << "\"/>" << std::endl
            << "    <Geometry GeometryType=\"VXVYVZ\">" << std::endl
            << "      <DataItem Format=\"Binary\" NumberType=\"Float\" Precision=\"4\" Dimensions=\"" << NUMBER_OF_BASIS_FUNCTIONS * globals.X << "\">" << std::endl
            << "        " << xGridName.str() << std::endl
            << "      </DataItem>" << std::endl
            << "      <DataItem Format=\"Binary\" NumberType=\"Float\" Precision=\"4\" Dimensions=\"" << NUMBER_OF_BASIS_FUNCTIONS * globals.Y << "\">" << std::endl
            << "        " << yGridName.str() << std::endl
            << "      </DataItem>" << std::endl
            << "      <DataItem Format=\"XML\" Dimensions=\"1\">0.0</DataItem>" << std::endl
            << "    </Geometry>" << std::endl
            << "    <Grid Name=\"TimeSeries\" GridType=\"Collection\" CollectionType=\"Temporal\">" << std::endl;

    int gridSize = NUMBER_OF_BASIS_FUNCTIONS * NUMBER_OF_BASIS_FUNCTIONS * globals.X * globals.Y;
    m_pressure = new float[gridSize];
    m_uvel = new float[gridSize];
    m_vvel = new float[gridSize];

    float* xGrid = new float[NUMBER_OF_BASIS_FUNCTIONS*globals.X];
    float* yGrid = new float[NUMBER_OF_BASIS_FUNCTIONS*globals.Y];

    for (unsigned x = 0; x < globals.X; ++x) {
      for (unsigned i = 0; i < NUMBER_OF_BASIS_FUNCTIONS; ++i) {
        xGrid[x*NUMBER_OF_BASIS_FUNCTIONS + i] = (x + lina::LGLNodes[i])*globals.hx;
      }
    }

    for (unsigned y = 0; y < globals.Y; ++y) {
      for (unsigned i = 0; i < NUMBER_OF_BASIS_FUNCTIONS; ++i) {
        yGrid[y*NUMBER_OF_BASIS_FUNCTIONS + i] = (y + lina::LGLNodes[i])*globals.hy;
      }
    }

    FILE* xGridGile = fopen((m_dirName + xGridName.str()).c_str(), "wb");
    fwrite(xGrid, sizeof(float), NUMBER_OF_BASIS_FUNCTIONS*globals.X, xGridGile);
    fclose(xGridGile);

    FILE* yGridGile = fopen((m_dirName + yGridName.str()).c_str(), "wb");
    fwrite(yGrid, sizeof(float), NUMBER_OF_BASIS_FUNCTIONS*globals.Y, yGridGile);
    fclose(yGridGile);

    delete[] xGrid;
    delete[] yGrid;
  }
}

WaveFieldWriter::~WaveFieldWriter()
{
  if (!m_baseName.empty()) {
    m_xdmf  << "    </Grid>" << std::endl
            << "  </Domain>" << std::endl
            << "</Xdmf>" << std::endl;
    m_xdmf.close();

    delete[] m_pressure;
    delete[] m_uvel;
    delete[] m_vvel;
  }
}

void WaveFieldWriter::writeTimestep(double time, Grid<DegreesOfFreedom>& degreesOfFreedomGrid, bool forceWrite)
{
  if (!m_baseName.empty() && (time >= m_lastTime + m_interval || forceWrite)) {
    m_lastTime = time;
    
    std::stringstream pressureFileName, uvelFileName, vvelFileName;
    pressureFileName << m_baseName << "_pressure" << m_step << ".bin";
    uvelFileName << m_baseName << "_u" << m_step << ".bin";
    vvelFileName << m_baseName << "_v" << m_step << ".bin";
    
    m_xdmf  << "      <Grid Name=\"step_" << m_step << "\" GridType=\"Uniform\">" << std::setw(0) << std::endl
            << "        <Topology Reference=\"/Xdmf/Domain/Topology[1]\"/>" << std::endl
            << "        <Geometry Reference=\"/Xdmf/Domain/Geometry[1]\"/>" << std::endl
            << "        <Time Value=\"" << time << "\"/>" << std::endl
            << "        <Attribute Name=\"pressure\" Center=\"Node\">" << std::endl
            << "          <DataItem Format=\"Binary\" NumberType=\"Float\" Precision=\"4\" Dimensions=\"1 " << NUMBER_OF_BASIS_FUNCTIONS * degreesOfFreedomGrid.Y() << " " << NUMBER_OF_BASIS_FUNCTIONS * degreesOfFreedomGrid.X() << "\">" << std::endl
            << "            " << pressureFileName.str() << std::endl
            << "          </DataItem>" << std::endl
            << "        </Attribute>" << std::endl
            << "        <Attribute Name=\"u\" Center=\"Node\">" << std::endl
            << "          <DataItem Format=\"Binary\" NumberType=\"Float\" Precision=\"4\" Dimensions=\"1 " << NUMBER_OF_BASIS_FUNCTIONS * degreesOfFreedomGrid.Y() << " " << NUMBER_OF_BASIS_FUNCTIONS * degreesOfFreedomGrid.X() << "\">" << std::endl
            << "            " << uvelFileName.str() << std::endl
            << "          </DataItem>" << std::endl
            << "       </Attribute>" << std::endl
            << "        <Attribute Name=\"v\" Center=\"Node\">" << std::endl
            << "          <DataItem Format=\"Binary\" NumberType=\"Float\" Precision=\"4\" Dimensions=\"1 " << NUMBER_OF_BASIS_FUNCTIONS * degreesOfFreedomGrid.Y() << " " << NUMBER_OF_BASIS_FUNCTIONS * degreesOfFreedomGrid.X() << "\">" << std::endl
            << "            " << vvelFileName.str() << std::endl
            << "          </DataItem>" << std::endl
            << "        </Attribute>" << std::endl
            << "      </Grid>" << std::endl;

    for (int y = 0; y < degreesOfFreedomGrid.Y(); ++y) {
      for (int x = 0; x < degreesOfFreedomGrid.X(); ++x) {
        DegreesOfFreedom& degreesOfFreedom  = degreesOfFreedomGrid.get(x, y);
        auto dofs = lina::init::Q::view::create(degreesOfFreedom);
        for (int j = 0; j < NUMBER_OF_BASIS_FUNCTIONS; ++j) {
          for (int i = 0; i < NUMBER_OF_BASIS_FUNCTIONS; ++i) {
            unsigned targetIndex = (y*NUMBER_OF_BASIS_FUNCTIONS+j)*NUMBER_OF_BASIS_FUNCTIONS*degreesOfFreedomGrid.X() + (x*NUMBER_OF_BASIS_FUNCTIONS+i);
            m_pressure[targetIndex] = dofs(i,j,0);
            m_uvel[targetIndex] = dofs(i,j,1);
            m_vvel[targetIndex] = dofs(i,j,2);
          }
        }
      }
    }

    unsigned subGridSize = NUMBER_OF_BASIS_FUNCTIONS * NUMBER_OF_BASIS_FUNCTIONS;
    FILE* pressureFile = fopen((m_dirName + pressureFileName.str()).c_str(), "wb");
    fwrite(m_pressure, sizeof(float), subGridSize*degreesOfFreedomGrid.X()*degreesOfFreedomGrid.Y(), pressureFile);
    fclose(pressureFile);
    
    FILE* uFile = fopen((m_dirName + uvelFileName.str()).c_str(), "wb");
    fwrite(m_uvel, sizeof(float), subGridSize*degreesOfFreedomGrid.X()*degreesOfFreedomGrid.Y(), uFile);
    fclose(uFile);
    
    FILE* vFile = fopen((m_dirName + vvelFileName.str()).c_str(), "wb");
    fwrite(m_vvel, sizeof(float), subGridSize*degreesOfFreedomGrid.X()*degreesOfFreedomGrid.Y(), vFile);
    fclose(vFile);
    
    ++m_step;
  }
}
