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
cmdLineArgs = cmdLineParser.parse_args()

arch = useArchitectureIdentifiedBy(cmdLineArgs.arch)

order = int(cmdLineArgs.order)
degree = order-1
numberOf1DBasisFunctions = order
numberOfQuadraturePoints = order+1
numberOfQuantities = 3

qShape = (numberOf1DBasisFunctions, numberOf1DBasisFunctions, numberOfQuantities)
qShape1 = (numberOf1DBasisFunctions, numberOfQuantities)

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

initialCond = Tensor('initialCond', (numberOfQuadraturePoints, numberOfQuadraturePoints, numberOfQuantities))

# Flux solver
fluxSolver = Tensor('fluxSolver', (numberOfQuantities, numberOfQuantities))
fluxSolverNeighbour = Tensor('fluxSolverNeighbour', (numberOfQuantities, numberOfQuantities))
T = Tensor('T', (numberOfQuantities, numberOfQuantities))
TT = Tensor('TT', (numberOfQuantities, numberOfQuantities))
Apm = Tensor('Apm', (numberOfQuantities, numberOfQuantities))

# Kernels
g = Generator(arch)

## Main kernels
volume = (Q['xyp'] <= Q['xyp'] + db.kDivM['xl'] * I['lyq'] * db.star[0]['qp'] + db.kDivMT['my'] * I['xmq'] * db.star[1]['qp'])
g.add('volume', volume)

def evaluateSide(dim,side):
  if dim == 0:
    return Q1['yp'] <= db.F[side]['l'] * I['lyp']
  return Q1['xp'] <= db.F[side]['m'] * I['xmp']

def flux(dim,side):
  if dim == 0:
    return Q['xyp'] <= Q['xyp'] + db.FDivM[side]['x'] * (Q1['yq'] * fluxSolver['qp'] + Q1Neighbour['yq'] * fluxSolverNeighbour['qp'])
  return Q['xyp'] <= Q['xyp'] + db.FDivM[side]['y'] * (Q1['xq'] * fluxSolver['qp'] + Q1Neighbour['xq'] * fluxSolverNeighbour['qp'])

g.addFamily('evaluateSide', simpleParameterSpace(2,2), evaluateSide)
g.addFamily('flux', simpleParameterSpace(2,2), flux)

power = Scalar('power')
dQcur = Tensor('dQcur', qShape)
dQnext = Tensor('dQnext', qShape)
g.add('derivative', dQnext['xyp'] <= db.kTDivM['xl'] * dQcur['lyq'] * db.star[0]['qp'] + db.kTDivMT['my'] * dQcur['xmq'] * db.star[1]['qp'])
g.add('derivativeTaylorExpansion(0)', I['xyp'] <= power * Q['xyp'])
g.add('derivativeTaylorExpansion(1)', I['xyp'] <= I['xyp'] + power * dQnext['xyp'])



## Initialization kernels
fluxScale = Scalar('fluxScale')
computeFluxSolver = fluxSolver['ij'] <= fluxScale * T['ik'] * Apm['kl'] * TT['lj']
g.add('computeFluxSolver', computeFluxSolver)

quadrature = Q['xyp'] <= db.quadrature['xl'] * db.quadrature['ym'] * initialCond['lmp']
g.add('quadrature', quadrature)

class MyPSpaMM(PSpaMM):
  def preference(self, m, n, k, sparseA, sparseB, transA, transB, alpha, beta):
    pref = super().preference(m, n, k, sparseA, sparseB, transA, transB, alpha, beta)
    if pref >= Preference.HIGH and n <= 3:
      return Preference.HIGHEST
    return pref

  def blockSize(self, m, n, k):
    if n <= 3:
      return {'bm': min(m, 32), 'bn': n, 'bk': 1}
    return super().blockSize(m, n, k)

# Generate code
gemmTool = GeneratorCollection([LIBXSMM(arch), MyPSpaMM(arch), MKL(arch)])
g.generate(cmdLineArgs.outputDir, 'lina', gemmTool)
