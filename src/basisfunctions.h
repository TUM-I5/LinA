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
static double basisFunction12(double xi) {
  double phi = 0.2704156e7 * pow(xi, 0.12e2) - 0.16224936e8 * pow(xi, 0.11e2) + 0.42678636e8 * pow(xi, 0.10e2) - 0.64664600e8 * pow(xi, 0.9e1) + 0.62355150e8 * pow(xi, 0.8e1) - 0.39907296e8 * pow(xi, 0.7e1) + 0.17153136e8 * pow(xi, 0.6e1) - 0.4900896e7 * pow(xi, 0.5e1) + 0.900900e6 * pow(xi, 0.4e1) - 0.100100e6 * pow(xi, 0.3e1) + 0.6006e4 * xi * xi - 0.156e3 * xi + 0.1e1;
  return phi;
}
static double basisFunction13(double xi) {
  double phi = 0.10400600e8 * pow(xi, 0.13e2) - 0.67603900e8 * pow(xi, 0.12e2) + 0.194699232e9 * pow(xi, 0.11e2) - 0.327202876e9 * pow(xi, 0.10e2) + 0.355655300e9 * pow(xi, 0.9e1) - 0.261891630e9 * pow(xi, 0.8e1) + 0.133024320e9 * pow(xi, 0.7e1) - 0.46558512e8 * pow(xi, 0.6e1) + 0.11027016e8 * pow(xi, 0.5e1) - 0.1701700e7 * pow(xi, 0.4e1) + 0.160160e6 * pow(xi, 0.3e1) - 0.8190e4 * xi * xi + 0.182e3 * xi - 0.1e1;
  return phi;
}
static double basisFunction14(double xi) {
  double phi = 0.40116600e8 * pow(xi, 0.14e2) - 0.280816200e9 * pow(xi, 0.13e2) + 0.878850700e9 * pow(xi, 0.12e2) - 0.1622493600e10 * pow(xi, 0.11e2) + 0.1963217256e10 * pow(xi, 0.10e2) - 0.1636014380e10 * pow(xi, 0.9e1) + 0.960269310e9 * pow(xi, 0.8e1) - 0.399072960e9 * pow(xi, 0.7e1) + 0.116396280e9 * pow(xi, 0.6e1) - 0.23279256e8 * pow(xi, 0.5e1) + 0.3063060e7 * pow(xi, 0.4e1) - 0.247520e6 * pow(xi, 0.3e1) + 0.10920e5 * xi * xi - 0.210e3 * xi + 0.1e1;
  return phi;
}
static double basisFunction15(double xi) {
  double phi = 0.155117520e9 * pow(xi, 0.15e2) - 0.1163381400e10 * pow(xi, 0.14e2) + 0.3931426800e10 * pow(xi, 0.13e2) - 0.7909656300e10 * pow(xi, 0.12e2) + 0.10546208400e11 * pow(xi, 0.11e2) - 0.9816086280e10 * pow(xi, 0.10e2) + 0.6544057520e10 * pow(xi, 0.9e1) - 0.3155170590e10 * pow(xi, 0.8e1) + 0.1097450640e10 * pow(xi, 0.7e1) - 0.271591320e9 * pow(xi, 0.6e1) + 0.46558512e8 * pow(xi, 0.5e1) - 0.5290740e7 * pow(xi, 0.4e1) + 0.371280e6 * pow(xi, 0.3e1) - 0.14280e5 * xi * xi + 0.240e3 * xi - 0.1e1;
  return phi;
}
static double basisFunction16(double xi) {
  double phi = 0.601080390e9 * pow(xi, 0.16e2) - 0.4808643120e10 * pow(xi, 0.15e2) + 0.17450721000e11 * pow(xi, 0.14e2) - 0.38003792400e11 * pow(xi, 0.13e2) + 0.55367594100e11 * pow(xi, 0.12e2) - 0.56949525360e11 * pow(xi, 0.11e2) + 0.42536373880e11 * pow(xi, 0.10e2) - 0.23371634000e11 * pow(xi, 0.9e1) + 0.9465511770e10 * pow(xi, 0.8e1) - 0.2804596080e10 * pow(xi, 0.7e1) + 0.597500904e9 * pow(xi, 0.6e1) - 0.88884432e8 * pow(xi, 0.5e1) + 0.8817900e7 * pow(xi, 0.4e1) - 0.542640e6 * pow(xi, 0.3e1) + 0.18360e5 * xi * xi - 0.272e3 * xi + 0.1e1;
  return phi;
}
static double basisFunction17(double xi) {
  double phi = 0.2333606220e10 * pow(xi, 0.17e2) - 0.19835652870e11 * pow(xi, 0.16e2) + 0.76938289920e11 * pow(xi, 0.15e2) - 0.180324117000e12 * pow(xi, 0.14e2) + 0.285028443000e12 * pow(xi, 0.13e2) - 0.321132045780e12 * pow(xi, 0.12e2) + 0.265764451680e12 * pow(xi, 0.11e2) - 0.164068870680e12 * pow(xi, 0.10e2) + 0.75957810500e11 * pow(xi, 0.9e1) - 0.26293088250e11 * pow(xi, 0.8e1) + 0.6731030592e10 * pow(xi, 0.7e1) - 0.1249320072e10 * pow(xi, 0.6e1) + 0.162954792e9 * pow(xi, 0.5e1) - 0.14244300e8 * pow(xi, 0.4e1) + 0.775200e6 * pow(xi, 0.3e1) - 0.23256e5 * xi * xi + 0.306e3 * xi - 0.1e1;
  return phi;
}
static double basisFunction18(double xi) {
  double phi = 0.9075135300e10 * pow(xi, 0.18e2) - 0.81676217700e11 * pow(xi, 0.17e2) + 0.337206098790e12 * pow(xi, 0.16e2) - 0.846321189120e12 * pow(xi, 0.15e2) + 0.1442592936000e13 * pow(xi, 0.14e2) - 0.1767176346600e13 * pow(xi, 0.13e2) + 0.1605660228900e13 * pow(xi, 0.12e2) - 0.1101024156960e13 * pow(xi, 0.11e2) + 0.574241047380e12 * pow(xi, 0.10e2) - 0.227873431500e12 * pow(xi, 0.9e1) + 0.68362029450e11 * pow(xi, 0.8e1) - 0.15297796800e11 * pow(xi, 0.7e1) + 0.2498640144e10 * pow(xi, 0.6e1) - 0.288304632e9 * pow(xi, 0.5e1) + 0.22383900e8 * pow(xi, 0.4e1) - 0.1085280e7 * pow(xi, 0.3e1) + 0.29070e5 * xi * xi - 0.342e3 * xi + 0.1e1;
  return phi;
}
static double basisFunction19(double xi) {
  double phi = 0.35345263800e11 * pow(xi, 0.19e2) - 0.335780006100e12 * pow(xi, 0.18e2) + 0.1470171918600e13 * pow(xi, 0.17e2) - 0.3934071152550e13 * pow(xi, 0.16e2) + 0.7193730107520e13 * pow(xi, 0.15e2) - 0.9521113377600e13 * pow(xi, 0.14e2) + 0.9424940515200e13 * pow(xi, 0.13e2) - 0.7110781013700e13 * pow(xi, 0.12e2) + 0.4128840588600e13 * pow(xi, 0.11e2) - 0.1850332263780e13 * pow(xi, 0.10e2) + 0.638045608200e12 * pow(xi, 0.9e1) - 0.167797708650e12 * pow(xi, 0.8e1) + 0.33145226400e11 * pow(xi, 0.7e1) - 0.4805077200e10 * pow(xi, 0.6e1) + 0.494236512e9 * pow(xi, 0.5e1) - 0.34321980e8 * pow(xi, 0.4e1) + 0.1492260e7 * pow(xi, 0.3e1) - 0.35910e5 * xi * xi + 0.380e3 * xi - 0.1e1;
  return phi;
}
static double basisFunction20(double xi) {
  double phi = 0.137846528820e12 * pow(xi, 0.20e2) - 0.1378465288200e13 * pow(xi, 0.19e2) + 0.6379820115900e13 * pow(xi, 0.18e2) - 0.18132120329400e14 * pow(xi, 0.17e2) + 0.35406640372950e14 * pow(xi, 0.16e2) - 0.50356110752640e14 * pow(xi, 0.15e2) + 0.53952975806400e14 * pow(xi, 0.14e2) - 0.44431862428800e14 * pow(xi, 0.13e2) + 0.28443124054800e14 * pow(xi, 0.12e2) - 0.14221562027400e14 * pow(xi, 0.11e2) + 0.5550996791340e13 * pow(xi, 0.10e2) - 0.1682120239800e13 * pow(xi, 0.9e1) + 0.391527986850e12 * pow(xi, 0.8e1) - 0.68840085600e11 * pow(xi, 0.7e1) + 0.8923714800e10 * pow(xi, 0.6e1) - 0.823727520e9 * pow(xi, 0.5e1) + 0.51482970e8 * pow(xi, 0.4e1) - 0.2018940e7 * pow(xi, 0.3e1) + 0.43890e5 * xi * xi - 0.420e3 * xi + 0.1e1;
  return phi;
}
static double basisFunction21(double xi) {
  double phi = 0.538257874440e12 * pow(xi, 0.21e2) - 0.5651707681620e13 * pow(xi, 0.20e2) + 0.27569305764000e14 * pow(xi, 0.19e2) - 0.82937661506700e14 * pow(xi, 0.18e2) + 0.172255143129300e15 * pow(xi, 0.17e2) - 0.262009138759830e15 * pow(xi, 0.16e2) + 0.302136664515840e15 * pow(xi, 0.15e2) - 0.269764879032000e15 * pow(xi, 0.14e2) + 0.188835415322400e15 * pow(xi, 0.13e2) - 0.104291454867600e15 * pow(xi, 0.12e2) + 0.45508998487680e14 * pow(xi, 0.11e2) - 0.15643718230140e14 * pow(xi, 0.10e2) + 0.4205300599500e13 * pow(xi, 0.9e1) - 0.873408586050e12 * pow(xi, 0.8e1) + 0.137680171200e12 * pow(xi, 0.7e1) - 0.16062686640e11 * pow(xi, 0.6e1) + 0.1338557220e10 * pow(xi, 0.5e1) - 0.75710250e8 * pow(xi, 0.4e1) + 0.2691920e7 * pow(xi, 0.3e1) - 0.53130e5 * xi * xi + 0.462e3 * xi - 0.1e1;
  return phi;
}
static double basisFunction22(double xi) {
  double phi = 0.2104098963720e13 * pow(xi, 0.22e2) - 0.23145088600920e14 * pow(xi, 0.21e2) + 0.118685861314020e15 * pow(xi, 0.20e2) - 0.376780512108000e15 * pow(xi, 0.19e2) + 0.829376615067000e15 * pow(xi, 0.18e2) - 0.1343590116408540e16 * pow(xi, 0.17e2) + 0.1659391212145590e16 * pow(xi, 0.16e2) - 0.1597008083869440e16 * pow(xi, 0.15e2) + 0.1213941955644000e16 * pow(xi, 0.14e2) - 0.734359948476000e15 * pow(xi, 0.13e2) + 0.354590946549840e15 * pow(xi, 0.12e2) - 0.136526995463040e15 * pow(xi, 0.11e2) + 0.41716581947040e14 * pow(xi, 0.10e2) - 0.10028024506500e14 * pow(xi, 0.9e1) + 0.1871589827250e13 * pow(xi, 0.8e1) - 0.266181664320e12 * pow(xi, 0.7e1) + 0.28109701620e11 * pow(xi, 0.6e1) - 0.2125943820e10 * pow(xi, 0.5e1) + 0.109359250e9 * pow(xi, 0.4e1) - 0.3542000e7 * pow(xi, 0.3e1) + 0.63756e5 * xi * xi - 0.506e3 * xi + 0.1e1;
  return phi;
}
static double basisFunction23(double xi) {
  double phi = 0.8233430727600e13 * pow(xi, 0.23e2) - 0.94684453367400e14 * pow(xi, 0.22e2) + 0.509191949220240e15 * pow(xi, 0.21e2) - 0.1701164012167620e16 * pow(xi, 0.20e2) + 0.3956195377134000e16 * pow(xi, 0.19e2) - 0.6800888243549400e16 * pow(xi, 0.18e2) + 0.8957267442723600e16 * pow(xi, 0.17e2) - 0.9245179610525430e16 * pow(xi, 0.16e2) + 0.7585788398379840e16 * pow(xi, 0.15e2) - 0.4990650262092000e16 * pow(xi, 0.14e2) + 0.2643695814513600e16 * pow(xi, 0.13e2) - 0.1128243920840400e16 * pow(xi, 0.12e2) + 0.386826487145280e15 * pow(xi, 0.11e2) - 0.105895938788640e15 * pow(xi, 0.10e2) + 0.22921198872000e14 * pow(xi, 0.9e1) - 0.3867952309650e13 * pow(xi, 0.8e1) + 0.499090620600e12 * pow(xi, 0.7e1) - 0.47951843940e11 * pow(xi, 0.6e1) + 0.3307023720e10 * pow(xi, 0.5e1) - 0.155405250e9 * pow(xi, 0.4e1) + 0.4604600e7 * pow(xi, 0.3e1) - 0.75900e5 * xi * xi + 0.552e3 * xi - 0.1e1;
  return phi;
}
static double basisFunction24(double xi) {
  double phi = 0.32247603683100e14 * pow(xi, 0.24e2) - 0.386971244197200e15 * pow(xi, 0.23e2) + 0.2177742427450200e16 * pow(xi, 0.22e2) - 0.7637879238303600e16 * pow(xi, 0.21e2) + 0.18712804133843820e17 * pow(xi, 0.20e2) - 0.34023280243352400e17 * pow(xi, 0.19e2) + 0.47606217704845800e17 * pow(xi, 0.18e2) - 0.52463995021666800e17 * pow(xi, 0.17e2) + 0.46225898052627150e17 * pow(xi, 0.16e2) - 0.32871749726312640e17 * pow(xi, 0.15e2) + 0.18964470995949600e17 * pow(xi, 0.14e2) - 0.8892431376091200e16 * pow(xi, 0.13e2) + 0.3384731762521200e16 * pow(xi, 0.12e2) - 0.1041455926929600e16 * pow(xi, 0.11e2) + 0.257175851343840e15 * pow(xi, 0.10e2) - 0.50426637518400e14 * pow(xi, 0.9e1) + 0.7735904619300e13 * pow(xi, 0.8e1) - 0.910106425800e12 * pow(xi, 0.7e1) + 0.79919739900e11 * pow(xi, 0.6e1) - 0.5047562520e10 * pow(xi, 0.5e1) + 0.217567350e9 * pow(xi, 0.4e1) - 0.5920200e7 * pow(xi, 0.3e1) + 0.89700e5 * xi * xi - 0.600e3 * xi + 0.1e1;
  return phi;
}
static double basisFunction25(double xi) {
  double phi = 0.126410606437752e15 * pow(xi, 0.25e2) - 0.1580132580471900e16 * pow(xi, 0.24e2) + 0.9287309860732800e16 * pow(xi, 0.23e2) - 0.34117964696719800e17 * pow(xi, 0.22e2) + 0.87835611240491400e17 * pow(xi, 0.21e2) - 0.168415237204594380e18 * pow(xi, 0.20e2) + 0.249504055117917600e18 * pow(xi, 0.19e2) - 0.292438194472624200e18 * pow(xi, 0.18e2) + 0.275435973863750700e18 * pow(xi, 0.17e2) - 0.210584646684190350e18 * pow(xi, 0.16e2) + 0.131486998905250560e18 * pow(xi, 0.15e2) - 0.67237669894730400e17 * pow(xi, 0.14e2) + 0.28159366024288800e17 * pow(xi, 0.13e2) - 0.9633467324098800e16 * pow(xi, 0.12e2) + 0.2678029526390400e16 * pow(xi, 0.11e2) - 0.600076986468960e15 * pow(xi, 0.10e2) + 0.107156604726600e15 * pow(xi, 0.9e1) - 0.15016756025700e14 * pow(xi, 0.8e1) + 0.1617966979200e13 * pow(xi, 0.7e1) - 0.130395365100e12 * pow(xi, 0.6e1) + 0.7571343780e10 * pow(xi, 0.5e1) - 0.300450150e9 * pow(xi, 0.4e1) + 0.7534800e7 * pow(xi, 0.3e1) - 0.105300e6 * xi * xi + 0.650e3 * xi - 0.1e1;
  return phi;
}
static double basisFunction26(double xi) {
  double phi = 0.495918532948104e15 * pow(xi, 0.26e2) - 0.6446940928325352e16 * pow(xi, 0.25e2) + 0.39503314511797500e17 * pow(xi, 0.24e2) - 0.151692727725302400e18 * pow(xi, 0.23e2) + 0.409415576360637600e18 * pow(xi, 0.22e2) - 0.825654745660619160e18 * pow(xi, 0.21e2) + 0.1291183485235223580e19 * pow(xi, 0.20e2) - 0.1603954640043756000e19 * pow(xi, 0.19e2) + 0.1608410069599433100e19 * pow(xi, 0.18e2) - 0.1315971875126808900e19 * pow(xi, 0.17e2) + 0.884455516073599470e18 * pow(xi, 0.16e2) - 0.490087905010479360e18 * pow(xi, 0.15e2) + 0.224125566315768000e18 * pow(xi, 0.14e2) - 0.84478098072866400e17 * pow(xi, 0.13e2) + 0.26147982736839600e17 * pow(xi, 0.12e2) - 0.6605806165096320e16 * pow(xi, 0.11e2) + 0.1350173219555160e16 * pow(xi, 0.10e2) - 0.220616539143000e15 * pow(xi, 0.9e1) + 0.28364983604100e14 * pow(xi, 0.8e1) - 0.2810153174400e13 * pow(xi, 0.7e1) + 0.208632584160e12 * pow(xi, 0.6e1) - 0.11176745580e11 * pow(xi, 0.5e1) + 0.409704750e9 * pow(xi, 0.4e1) - 0.9500400e7 * pow(xi, 0.3e1) + 0.122850e6 * xi * xi - 0.702e3 * xi + 0.1e1;
  return phi;
}
static double basisFunction27(double xi) {
  double phi = 0.1946939425648112e16 * pow(xi, 0.27e2) - 0.26283682246249512e17 * pow(xi, 0.26e2) + 0.167620464136459152e18 * pow(xi, 0.25e2) - 0.671556346700557500e18 * pow(xi, 0.24e2) + 0.1896159096566280000e19 * pow(xi, 0.23e2) - 0.4012272648334248480e19 * pow(xi, 0.22e2) + 0.6605237965284953280e19 * pow(xi, 0.21e2) - 0.8669374829436501180e19 * pow(xi, 0.20e2) + 0.9222739180251597000e19 * pow(xi, 0.19e2) - 0.8042050347997165500e19 * pow(xi, 0.18e2) + 0.5790276250557959160e19 * pow(xi, 0.17e2) - 0.3457417017378616110e19 * pow(xi, 0.16e2) + 0.1715307667536677760e19 * pow(xi, 0.15e2) - 0.706857555303576000e18 * pow(xi, 0.14e2) + 0.241365994493904000e18 * pow(xi, 0.13e2) - 0.67984755115782960e17 * pow(xi, 0.12e2) + 0.15688789642103760e17 * pow(xi, 0.11e2) - 0.2938612301384760e16 * pow(xi, 0.10e2) + 0.441233078286000e15 * pow(xi, 0.9e1) - 0.52251285586500e14 * pow(xi, 0.8e1) + 0.4777260396480e13 * pow(xi, 0.7e1) - 0.327851203680e12 * pow(xi, 0.6e1) + 0.16257084480e11 * pow(xi, 0.5e1) - 0.552210750e9 * pow(xi, 0.4e1) + 0.11875500e8 * pow(xi, 0.3e1) - 0.142506e6 * xi * xi + 0.756e3 * xi - 0.1e1;
  return phi;
}
static double basisFunction28(double xi) {
  double phi = 0.7648690600760440e16 * pow(xi, 0.28e2) - 0.107081668410646160e18 * pow(xi, 0.27e2) + 0.709659420648736824e18 * pow(xi, 0.26e2) - 0.2961294866410778352e19 * pow(xi, 0.25e2) + 0.8730232507107247500e19 * pow(xi, 0.24e2) - 0.19340822784976056000e20 * pow(xi, 0.23e2) + 0.33435605402785404000e20 * pow(xi, 0.22e2) - 0.46236665756994672960e20 * pow(xi, 0.21e2) + 0.52016248976619007080e20 * pow(xi, 0.20e2) - 0.48163193496869451000e20 * pow(xi, 0.19e2) + 0.36993431600786961300e20 * pow(xi, 0.18e2) - 0.23687493752282560200e20 * pow(xi, 0.17e2) + 0.12677195730388259070e20 * pow(xi, 0.16e2) - 0.5673709977236703360e19 * pow(xi, 0.15e2) + 0.2120572665910728000e19 * pow(xi, 0.14e2) - 0.659733718283337600e18 * pow(xi, 0.13e2) + 0.169961887789457400e18 * pow(xi, 0.12e2) - 0.35991929178943920e17 * pow(xi, 0.11e2) + 0.6203737080701160e16 * pow(xi, 0.10e2) - 0.859243362978000e15 * pow(xi, 0.9e1) + 0.94052314055700e14 * pow(xi, 0.8e1) - 0.7962100660800e13 * pow(xi, 0.7e1) + 0.506679132960e12 * pow(xi, 0.6e1) - 0.23325382080e11 * pow(xi, 0.5e1) + 0.736281000e9 * pow(xi, 0.4e1) - 0.14725620e8 * pow(xi, 0.3e1) + 0.164430e6 * xi * xi - 0.812e3 * xi + 0.1e1;
  return phi;
}
static double basisFunction29(double xi) {
  double phi = 0.30067266499541040e17 * pow(xi, 0.29e2) - 0.435975364243345080e18 * pow(xi, 0.28e2) + 0.2998286715498092480e19 * pow(xi, 0.27e2) - 0.13010422711893508440e20 * pow(xi, 0.26e2) + 0.39977480696545507752e20 * pow(xi, 0.25e2) - 0.92540464575336823500e20 * pow(xi, 0.24e2) + 0.167620464136459152000e21 * pow(xi, 0.23e2) - 0.243602267934579372000e21 * pow(xi, 0.22e2) + 0.288979160981216706000e21 * pow(xi, 0.21e2) - 0.283199577761592371880e21 * pow(xi, 0.20e2) + 0.231183328784973364800e21 * pow(xi, 0.19e2) - 0.158062844112453380100e21 * pow(xi, 0.18e2) + 0.90802059383749814100e20 * pow(xi, 0.17e2) - 0.43882600605190127550e20 * pow(xi, 0.16e2) + 0.17831659928458210560e20 * pow(xi, 0.15e2) - 0.6078974975610753600e19 * pow(xi, 0.14e2) + 0.1731801010493761200e19 * pow(xi, 0.13e2) - 0.409908082315750200e18 * pow(xi, 0.12e2) + 0.79982064842097600e17 * pow(xi, 0.11e2) - 0.12733986639333960e17 * pow(xi, 0.10e2) + 0.1632562389658200e16 * pow(xi, 0.9e1) - 0.165711220002900e15 * pow(xi, 0.8e1) + 0.13028891990400e14 * pow(xi, 0.7e1) - 0.771033463200e12 * pow(xi, 0.6e1) + 0.33044291280e11 * pow(xi, 0.5e1) - 0.971890920e9 * pow(xi, 0.4e1) + 0.18123840e8 * pow(xi, 0.3e1) - 0.188790e6 * xi * xi + 0.870e3 * xi - 0.1e1;
  return phi;
}
static double basisFunction30(double xi) {
  double phi = 0.118264581564861424e18 * pow(xi, 0.30e2) - 0.1773968723472921360e19 * pow(xi, 0.29e2) + 0.12643285563057007320e20 * pow(xi, 0.28e2) - 0.56967447594463757120e20 * pow(xi, 0.27e2) + 0.182145917966509118160e21 * pow(xi, 0.26e2) - 0.439752287662000585272e21 * pow(xi, 0.25e2) + 0.832864181178031411500e21 * pow(xi, 0.24e2) - 0.1269126371318905008000e22 * pow(xi, 0.23e2) + 0.1583414741574765918000e22 * pow(xi, 0.22e2) - 0.1637548578893561334000e22 * pow(xi, 0.21e2) + 0.1415997888807961859400e22 * pow(xi, 0.20e2) - 0.1029816646405790443200e22 * pow(xi, 0.19e2) + 0.632251376449813520400e21 * pow(xi, 0.18e2) - 0.328284368541249327900e21 * pow(xi, 0.17e2) + 0.144185687702767561950e21 * pow(xi, 0.16e2) - 0.53494979785374631680e20 * pow(xi, 0.15e2) + 0.16717181182929572400e20 * pow(xi, 0.14e2) - 0.4380437850072454800e19 * pow(xi, 0.13e2) + 0.956452192070083800e18 * pow(xi, 0.12e2) - 0.172592876764526400e18 * pow(xi, 0.11e2) + 0.25467973278667920e17 * pow(xi, 0.10e2) - 0.3031901580793800e16 * pow(xi, 0.9e1) + 0.286228470914100e15 * pow(xi, 0.8e1) - 0.20959521897600e14 * pow(xi, 0.7e1) + 0.1156550194800e13 * pow(xi, 0.6e1) - 0.46262007792e11 * pow(xi, 0.5e1) + 0.1270934280e10 * pow(xi, 0.4e1) - 0.22151360e8 * pow(xi, 0.3e1) + 0.215760e6 * xi * xi - 0.930e3 * xi + 0.1e1;
  return phi;
}
static double basisFunction31(double xi) {
  double phi = 0.465428353255261088e18 * pow(xi, 0.31e2) - 0.7214139475456546864e19 * pow(xi, 0.30e2) + 0.53219061704187640800e20 * pow(xi, 0.29e2) - 0.248651282740121143960e21 * pow(xi, 0.28e2) + 0.826027990119724478240e21 * pow(xi, 0.27e2) - 0.2076463464818203947024e22 * pow(xi, 0.26e2) + 0.4104354684845338795872e22 * pow(xi, 0.25e2) - 0.6543932852113103947500e22 * pow(xi, 0.24e2) + 0.8566603006402608804000e22 * pow(xi, 0.23e2) - 0.9324553478162510406000e22 * pow(xi, 0.22e2) + 0.8515252610246518936800e22 * pow(xi, 0.21e2) - 0.6565081120836914075400e22 * pow(xi, 0.20e2) + 0.4290902693357460180000e22 * pow(xi, 0.19e2) - 0.2383101342003143269200e22 * pow(xi, 0.18e2) + 0.1125546406427140552800e22 * pow(xi, 0.17e2) - 0.451781821468671694110e21 * pow(xi, 0.16e2) + 0.153798066882952066080e21 * pow(xi, 0.15e2) - 0.44251361954813574000e20 * pow(xi, 0.14e2) + 0.10707736966843778400e20 * pow(xi, 0.13e2) - 0.2164602329421768600e19 * pow(xi, 0.12e2) + 0.362445041205505440e18 * pow(xi, 0.11e2) - 0.49723185925018320e17 * pow(xi, 0.10e2) + 0.5512548328716000e16 * pow(xi, 0.9e1) - 0.485343928941300e15 * pow(xi, 0.8e1) + 0.33185909671200e14 * pow(xi, 0.7e1) - 0.1711694288304e13 * pow(xi, 0.6e1) + 0.64055087712e11 * pow(xi, 0.5e1) - 0.1647507400e10 * pow(xi, 0.4e1) + 0.26898080e8 * pow(xi, 0.3e1) - 0.245520e6 * xi * xi + 0.992e3 * xi - 0.1e1;
  return phi;
}
static double (* const basisFunctions[])(double) = { basisFunction0,basisFunction1,basisFunction2,basisFunction3,basisFunction4,basisFunction5,basisFunction6,basisFunction7,basisFunction8,basisFunction9,basisFunction10,basisFunction11,basisFunction12,basisFunction13,basisFunction14,basisFunction15,basisFunction16,basisFunction17,basisFunction18,basisFunction19,basisFunction20,basisFunction21,basisFunction22,basisFunction23,basisFunction24,basisFunction25,basisFunction26,basisFunction27,basisFunction28,basisFunction29,basisFunction30,basisFunction31 };
#endif // BASISFUNCTIONS_H_
