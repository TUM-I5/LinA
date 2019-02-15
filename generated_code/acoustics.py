#!/usr/bin/env python3

import argparse
from yateto import *
from yateto.input import parseJSONMatrixFile, memoryLayoutFromFile
from yateto.ast.transformer import DeduceIndices, EquivalentSparsityPattern

cmdLineParser = argparse.ArgumentParser()
cmdLineParser.add_argument('--matricesDir')
cmdLineParser.add_argument('--outputDir')
cmdLineParser.add_argument('--arch')
cmdLineParser.add_argument('--order')
cmdLineParser.add_argument('--memLayout')
cmdLineParser.add_argument('--libxsmm', action='store_true')
cmdLineParser.add_argument('--pspamm', action='store_true')
cmdLineArgs = cmdLineParser.parse_args()

arch = useArchitectureIdentifiedBy(cmdLineArgs.arch)

order = int(cmdLineArgs.order)
degree = order-1
numberOf1DBasisFunctions = order
numberOfQuadraturePoints = order+1
numberOfQuantities = 4

qShape = (numberOf1DBasisFunctions, numberOf1DBasisFunctions, numberOf1DBasisFunctions, numberOfQuantities)
qShape1 = (numberOf1DBasisFunctions, numberOf1DBasisFunctions, numberOfQuantities)

clones = {
  'kDivM': ['kDivM', 'kDivMT'],
  'kTDivM': ['kTDivM', 'kTDivMT']
}
transpose = {'kDivMT', 'kTDivMT'}
alignStride = {'kDivM', 'kTDivM'}
db = parseJSONMatrixFile('{}/matrices_{}.json'.format(cmdLineArgs.matricesDir, degree), clones, transpose=transpose, alignStride=alignStride)
db.update( parseJSONMatrixFile('{}/star.json'.format(cmdLineArgs.matricesDir)) )
memoryLayoutFromFile(cmdLineArgs.memLayout, db, dict())

Q = Tensor('Q', qShape)
Q1 = Tensor('Q1', qShape1)
Q1Neighbour = Tensor('Q1Neighbour', qShape1)
dQ0 = Tensor('dQ(0)', qShape)
I = Tensor('I', qShape)

initialCond = Tensor('initialCond', (numberOfQuadraturePoints, numberOfQuadraturePoints, numberOfQuadraturePoints, numberOfQuantities))

# Flux solver
fluxSolver = Tensor('fluxSolver', (numberOfQuantities, numberOfQuantities))
fluxSolverNeighbour = Tensor('fluxSolverNeighbour', (numberOfQuantities, numberOfQuantities))
T = Tensor('T', (numberOfQuantities, numberOfQuantities))
TT = Tensor('TT', (numberOfQuantities, numberOfQuantities))
Apm = Tensor('Apm', (numberOfQuantities, numberOfQuantities))

# Kernels
g = Generator(arch)

## Main kernels
volume = (Q['xyzp'] <= Q['xyzp'] + db.kDivM['xl'] * I['lyzq'] * db.star[0]['qp'] + db.kDivMT['my'] * I['xmzq'] * db.star[1]['qp'] + db.kDivMT['nz'] * I['xynq'] * db.star[2]['qp'])
g.add('volume', volume)

def evaluateSide(dim,side):
  if dim == 0:
    return Q1['yzp'] <= db.F[side]['l'] * I['lyzp']
  if dim == 1:
    return Q1['xzp'] <= db.F[side]['m'] * I['xmzp']
  return Q1['xyp'] <= db.F[side]['n'] * I['xynp']

def flux(dim,side):
  if dim == 0:
    return Q['xyzp'] <= Q['xyzp'] + db.FDivM[side]['x'] * (Q1['yzq'] * fluxSolver['qp'] + Q1Neighbour['yzq'] * fluxSolverNeighbour['qp'])
  if dim == 1:
    return Q['xyzp'] <= Q['xyzp'] + db.FDivM[side]['y'] * (Q1['xzq'] * fluxSolver['qp'] + Q1Neighbour['xzq'] * fluxSolverNeighbour['qp'])
  return Q['xyzp'] <= Q['xyzp'] + db.FDivM[side]['z'] * (Q1['xyq'] * fluxSolver['qp'] + Q1Neighbour['xyq'] * fluxSolverNeighbour['qp'])

g.addFamily('evaluateSide', simpleParameterSpace(3,2), evaluateSide)
g.addFamily('flux', simpleParameterSpace(3,2), flux)

power = Scalar('power')
dQcur = Tensor('dQcur', qShape)
dQnext = Tensor('dQnext', qShape)
g.add('derivative', dQnext['xyzp'] <= db.kTDivM['xl'] * dQcur['lyzq'] * db.star[0]['qp'] + db.kTDivMT['my'] * dQcur['xmzq'] * db.star[1]['qp'] + db.kTDivMT['nz'] * dQcur['xynq'] * db.star[2]['qp'])
g.add('derivativeTaylorExpansion(0)', I['xyzp'] <= power * Q['xyzp'])
g.add('derivativeTaylorExpansion(1)', I['xyzp'] <= I['xyzp'] + power * dQnext['xyzp'])


## Initialization kernels
fluxScale = Scalar('fluxScale')
computeFluxSolver = fluxSolver['ij'] <= fluxScale * T['ik'] * Apm['kl'] * TT['lj']
g.add('computeFluxSolver', computeFluxSolver)

quadrature = Q['xyzp'] <= db.quadrature['xl'] * db.quadrature['ym'] * db.quadrature['zn'] * initialCond['lmnp']
g.add('quadrature', quadrature)

class MyPSpaMM(PSpaMM):
  def preference(self, m, n, k, sparseA, sparseB, transA, transB, alpha, beta):
    pref = super().preference(m, n, k, sparseA, sparseB, transA, transB, alpha, beta)
    if pref >= Preference.HIGH and n <= 3:
      return Preference.HIGHEST
    return pref

  def blockSize(self, m, n, k):
    if n <= 4:
      return {'bm': min(m, 32), 'bn': n, 'bk': 1}
    return super().blockSize(m, n, k)

# Generate code
generators = list()
if cmdLineArgs.libxsmm:
  generators.append(LIBXSMM(arch))
if cmdLineArgs.pspamm:
  generators.append(MyPSpaMM(arch))
generators.append(MKL(arch))
gemmTool = GeneratorCollection(generators)
g.generate(cmdLineArgs.outputDir, 'lina', gemmTool)
