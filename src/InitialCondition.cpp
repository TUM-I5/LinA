#include "InitialCondition.h"
#include "Quadrature.h"
#include "basisfunctions.h"
#include "generated_code/init.h"
#include "generated_code/kernel.h"

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

  lina::kernel::quadrature krnl;
  krnl.initialCond = icBuffer;
  krnl.quadrature = globalMatrices.quadrature;

  #pragma omp for collapse(2)  
  for (int y = 0; y < globals.Y; ++y) {
    for (int x = 0; x < globals.X; ++x) {
      DegreesOfFreedom& degreesOfFreedom = degreesOfFreedomGrid.get(x, y);
      Material& material = materialGrid.get(x, y);
      
      double scaledWavespeed = sqrt(2.) * material.wavespeed() / 2.;
      
      for (int i = 0; i < npoints; ++i) {
        double xi = (points[i]+1.)/2.;
        for (unsigned j = 0; j < npoints; ++j) {
          double eta = (points[j]+1.)/2.;

          double xp = xi*globals.hx + x*globals.hx;
          double yp = eta*globals.hy + y*globals.hy;
          double sn = sin(-2.*M_PI*xp - 2.*M_PI*yp);
          ic(i,j,0) = material.K0*sn;
          ic(i,j,1) = scaledWavespeed*sn;
          ic(i,j,2) = scaledWavespeed*sn;
        }
      }
      
      krnl.Q = degreesOfFreedom;
      krnl.execute();
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
  
  double area = globals.hx*globals.hy;

  #pragma omp parallel for collapse(2) reduction(+:l2error[:NUMBER_OF_QUANTITIES])
  for (int y = 0; y < globals.Y; ++y) {
    for (int x = 0; x < globals.X; ++x) {
      DegreesOfFreedom& degreesOfFreedom = degreesOfFreedomGrid.get(x, y);
      Material& material = materialGrid.get(x, y);
      
      auto dofs = lina::init::Q::view::create(degreesOfFreedom);

      double scaledWavespeed = sqrt(2.) * material.wavespeed() / 2.;
      double omega = 2.*sqrt(2.) * M_PI * material.wavespeed();
      
      for (int i = 0; i < npoints; ++i) {
        double xi = (points[i]+1.)/2.;
        for (unsigned j = 0; j < npoints; ++j) {
          double eta = (points[j]+1.)/2.;
          double weight = weights[i] * weights[j] / 4.;

          double xp = xi*globals.hx + x*globals.hx;
          double yp = eta*globals.hy + y*globals.hy;
          double sn = sin(omega*time-2.*M_PI*xp - 2.*M_PI*yp);
          double f[] = {material.K0*sn, scaledWavespeed*sn, scaledWavespeed*sn};

          double Q[NUMBER_OF_QUANTITIES] = {};
          for (unsigned l = 0; l < NUMBER_OF_BASIS_FUNCTIONS; ++l) {
            for (unsigned k = 0; k < NUMBER_OF_BASIS_FUNCTIONS; ++k) {
              double bf = (*basisFunctions[k])(xi) * (*basisFunctions[l])(eta);
              for (unsigned p = 0; p < NUMBER_OF_QUANTITIES; ++p) {
                Q[p] += dofs(k,l,p) * bf;
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
  for (unsigned q = 0; q < NUMBER_OF_QUANTITIES; ++q) {
    l2error[q] = sqrt(l2error[q]);
  }
}


void initSourcetermPhi(double xi, double eta, SourceTerm& sourceterm) {
  /*for (unsigned k = 0; k < NUMBER_OF_BASIS_FUNCTIONS; ++k) {
    sourceterm.phi[k] = basisFunctions[k](xi, eta) * GlobalMatrices::Minv[k*NUMBER_OF_BASIS_FUNCTIONS + k];
  }*/
}
