/** This file is generated. Do not edit. */
#ifndef BASISFUNCTIONS_H_
#define BASISFUNCTIONS_H_
#include <cmath>
static double basisFunction0(double xi) {
  double phi = 1;
  return phi;
}
static double basisFunction1(double xi) {
  double phi = 2 * xi - 1;
  return phi;
}
static double basisFunction2(double xi) {
  double phi = 6 * xi * xi - 6 * xi + 1;
  return phi;
}
static double basisFunction3(double xi) {
  double phi = 0.20e2 * pow(xi, 0.3e1) - 0.30e2 * xi * xi + 0.12e2 * xi - 0.1e1;
  return phi;
}
static double basisFunction4(double xi) {
  double phi = 0.70e2 * pow(xi, 0.4e1) - 0.140e3 * pow(xi, 0.3e1) + 0.90e2 * xi * xi - 0.20e2 * xi + 0.1e1;
  return phi;
}
static double basisFunction5(double xi) {
  double phi = 0.252e3 * pow(xi, 0.5e1) - 0.630e3 * pow(xi, 0.4e1) + 0.560e3 * pow(xi, 0.3e1) - 0.210e3 * xi * xi + 0.30e2 * xi - 0.1e1;
  return phi;
}
static double basisFunction6(double xi) {
  double phi = 0.924e3 * pow(xi, 0.6e1) - 0.2772e4 * pow(xi, 0.5e1) + 0.3150e4 * pow(xi, 0.4e1) - 0.1680e4 * pow(xi, 0.3e1) + 0.420e3 * xi * xi - 0.42e2 * xi + 0.1e1;
  return phi;
}
static double basisFunction7(double xi) {
  double phi = 0.3432e4 * pow(xi, 0.7e1) - 0.12012e5 * pow(xi, 0.6e1) + 0.16632e5 * pow(xi, 0.5e1) - 0.11550e5 * pow(xi, 0.4e1) + 0.4200e4 * pow(xi, 0.3e1) - 0.756e3 * xi * xi + 0.56e2 * xi - 0.1e1;
  return phi;
}
static double basisFunction8(double xi) {
  double phi = 0.12870e5 * pow(xi, 0.8e1) - 0.51480e5 * pow(xi, 0.7e1) + 0.84084e5 * pow(xi, 0.6e1) - 0.72072e5 * pow(xi, 0.5e1) + 0.34650e5 * pow(xi, 0.4e1) - 0.9240e4 * pow(xi, 0.3e1) + 0.1260e4 * xi * xi - 0.72e2 * xi + 0.1e1;
  return phi;
}
static double basisFunction9(double xi) {
  double phi = 0.48620e5 * pow(xi, 0.9e1) - 0.218790e6 * pow(xi, 0.8e1) + 0.411840e6 * pow(xi, 0.7e1) - 0.420420e6 * pow(xi, 0.6e1) + 0.252252e6 * pow(xi, 0.5e1) - 0.90090e5 * pow(xi, 0.4e1) + 0.18480e5 * pow(xi, 0.3e1) - 0.1980e4 * xi * xi + 0.90e2 * xi - 0.1e1;
  return phi;
}
static double basisFunction10(double xi) {
  double phi = 0.184756e6 * pow(xi, 0.10e2) - 0.923780e6 * pow(xi, 0.9e1) + 0.1969110e7 * pow(xi, 0.8e1) - 0.2333760e7 * pow(xi, 0.7e1) + 0.1681680e7 * pow(xi, 0.6e1) - 0.756756e6 * pow(xi, 0.5e1) + 0.210210e6 * pow(xi, 0.4e1) - 0.34320e5 * pow(xi, 0.3e1) + 0.2970e4 * xi * xi - 0.110e3 * xi + 0.1e1;
  return phi;
}
static double basisFunction11(double xi) {
  double phi = 0.705432e6 * pow(xi, 0.11e2) - 0.3879876e7 * pow(xi, 0.10e2) + 0.9237800e7 * pow(xi, 0.9e1) - 0.12471030e8 * pow(xi, 0.8e1) + 0.10501920e8 * pow(xi, 0.7e1) - 0.5717712e7 * pow(xi, 0.6e1) + 0.2018016e7 * pow(xi, 0.5e1) - 0.450450e6 * pow(xi, 0.4e1) + 0.60060e5 * pow(xi, 0.3e1) - 0.4290e4 * xi * xi + 0.132e3 * xi - 0.1e1;
  return phi;
}
static double (* const basisFunctions[])(double) = { basisFunction0,basisFunction1,basisFunction2,basisFunction3,basisFunction4,basisFunction5,basisFunction6,basisFunction7,basisFunction8,basisFunction9,basisFunction10,basisFunction11 };
#endif // BASISFUNCTIONS_H_
