#include "InitialCondition.h"
#include "Quadrature.h"
#include "basisfunctions.h"
#include "generated_code/init.h"
#include "generated_code/kernel.h"

void planeWave(Material const& material, double t, double xp, double yp, double zp, double f[NUMBER_OF_QUANTITIES]) {
  double k = 10.0 * M_PI;
  double omega = sqrt(3.) * k * material.wavespeed();
  double sn = sin(omega*t - k*xp - k*yp - k*zp);
  double scaledWavespeed = sqrt(3.) * material.wavespeed() / 3.;
  f[0] = material.K0*sn;
  f[1] = scaledWavespeed*sn;
  f[2] = scaledWavespeed*sn;
  f[3] = scaledWavespeed*sn;
}

void initialCondition(  GlobalConstants const& globals,
                        GlobalMatrices const& globalMatrices,
                        Grid<Material>& materialGrid,
                        Grid<DegreesOfFreedom>& degreesOfFreedomGrid  )
{
  int const npoints = CONVERGENCE_ORDER+1;
  double points[npoints];
  double weights[npoints];

  seissol::quadrature::GaussLegendre(points, weights, npoints);

#pragma omp parallel
{
  real icBuffer[lina::tensor::initialCond::Size] __attribute__((aligned(ALIGNMENT))) = {};
  auto ic = lina::init::initialCond::view::create(icBuffer);
  double f[NUMBER_OF_QUANTITIES];

  lina::kernel::quadrature krnl;
  krnl.initialCond = icBuffer;
  krnl.quadrature = globalMatrices.quadrature;

  #pragma omp for collapse(3)
  for (int z = 0; z < globals.Z; ++z) {
    for (int y = 0; y < globals.Y; ++y) {
      for (int x = 0; x < globals.X; ++x) {
        DegreesOfFreedom& degreesOfFreedom = degreesOfFreedomGrid.get(x, y, z);
        Material& material = materialGrid.get(x, y, z);
        
        for (int i = 0; i < npoints; ++i) {
          double xi = (points[i]+1.)/2.;
          for (unsigned j = 0; j < npoints; ++j) {
            double eta = (points[j]+1.)/2.;
            for (unsigned k = 0; k < npoints; ++k) {
              double zeta = (points[k]+1.)/2.;

              double xp = xi*globals.hx + x*globals.hx;
              double yp = eta*globals.hy + y*globals.hy;
              double zp = zeta*globals.hz + z*globals.hz;

              planeWave(material, 0.0, xp, yp, zp, f);

              for (int p = 0; p < NUMBER_OF_QUANTITIES; ++p) {
                ic(i,j,k,p) = f[p];
              }
            }
          }
        }
        
        krnl.Q = degreesOfFreedom;
        krnl.execute();
      }
    }
  }
}
}

void L2error( double time,
              GlobalConstants const& globals,
              Grid<Material>& materialGrid,
              Grid<DegreesOfFreedom>& degreesOfFreedomGrid,
              double l2error[NUMBER_OF_QUANTITIES]  )
{
  int const npoints = CONVERGENCE_ORDER+1;
  double points[npoints];
  double weights[npoints];
  
  seissol::quadrature::GaussLegendre(points, weights, npoints);

  memset(l2error, 0, NUMBER_OF_QUANTITIES * sizeof(double));
  
  double area = globals.hx*globals.hy*globals.hz;

  #pragma omp parallel for collapse(3) reduction(+:l2error[:NUMBER_OF_QUANTITIES])
  for (int z = 0; z < globals.Z; ++z) {
    for (int y = 0; y < globals.Y; ++y) {
      for (int x = 0; x < globals.X; ++x) {
        DegreesOfFreedom& degreesOfFreedom = degreesOfFreedomGrid.get(x, y, z);
        Material& material = materialGrid.get(x, y, z);
        
        auto dofs = lina::init::Q::view::create(degreesOfFreedom);
        
        for (int i = 0; i < npoints; ++i) {
          double xi = (points[i]+1.)/2.;
          for (unsigned j = 0; j < npoints; ++j) {
            double eta = (points[j]+1.)/2.;
            for (unsigned k = 0; k < npoints; ++k) {
              double zeta = (points[k]+1.)/2.;
              double weight = weights[i] * weights[j] * weights[k] / 8.;

              double xp = xi*globals.hx + x*globals.hx;
              double yp = eta*globals.hy + y*globals.hy;
              double zp = zeta*globals.hz + z*globals.hz;

              double f[NUMBER_OF_QUANTITIES];
              planeWave(material, time, xp, yp, zp, f);

              double Q[NUMBER_OF_QUANTITIES] = {};
              for (unsigned n = 0; n < NUMBER_OF_BASIS_FUNCTIONS; ++n) {
                for (unsigned m = 0; m < NUMBER_OF_BASIS_FUNCTIONS; ++m) {
                  for (unsigned l = 0; l < NUMBER_OF_BASIS_FUNCTIONS; ++l) {
                    double bf = lina::basisFunction(xi, l) * lina::basisFunction(eta, m) * lina::basisFunction(zeta, n);
                    for (unsigned p = 0; p < NUMBER_OF_QUANTITIES; ++p) {
                      Q[p] += dofs(l,m,n,p) * bf;
                    }
                  }
                }
              }
              for (unsigned q = 0; q < NUMBER_OF_QUANTITIES; ++q) {
                double diff = Q[q] - f[q];
                l2error[q] += diff * diff * weight * area;
              }
            }
          }
        }
      }
    }
  }
  for (unsigned q = 0; q < NUMBER_OF_QUANTITIES; ++q) {
    l2error[q] = sqrt(l2error[q]);
  }
}


void initSourcetermPhi(double xi, double eta, SourceTerm& sourceterm) {
  /*for (unsigned k = 0; k < NUMBER_OF_BASIS_FUNCTIONS; ++k) {
    sourceterm.phi[k] = basisFunctions[k](xi, eta) * GlobalMatrices::Minv[k*NUMBER_OF_BASIS_FUNCTIONS + k];
  }*/
}
