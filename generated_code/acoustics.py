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

Fnames = [tuple('{}({},{})'.format(name, i,j) for name in ['FDivM', 'FDivMT']) for i in range(2) for j in range(2)]
clones = {
  'kDivM': ['kDivM', 'kDivMT'],
  'kTDivM': ['kTDivM', 'kTDivMT']
}
clones.update({a: [a, b] for a, b in Fnames})
transpose = {'kDivMT', 'kTDivMT'} | set(b for a, b in Fnames)
alignStride = {'kDivM', 'kTDivM'}
db = parseJSONMatrixFile('{}/matrices_{}.json'.format(cmdLineArgs.matricesDir, degree), clones, transpose=transpose, alignStride=alignStride)
db.update( parseJSONMatrixFile('{}/star.json'.format(cmdLineArgs.matricesDir)) )
memoryLayoutFromFile(cmdLineArgs.memLayout, db, dict())

Q = Tensor('Q', qShape)
dQ0 = Tensor('dQ(0)', qShape)
I = Tensor('I', qShape)

initialCond = Tensor('initialCond', (numberOfQuadraturePoints, numberOfQuadraturePoints, numberOfQuantities))

# Flux solver
fluxSolver = Tensor('fluxSolver', (numberOfQuantities, numberOfQuantities))
T = Tensor('T', (numberOfQuantities, numberOfQuantities))
TT = Tensor('TT', (numberOfQuantities, numberOfQuantities))
Apm = Tensor('Apm', (numberOfQuantities, numberOfQuantities))

# Kernels
g = Generator(arch)

## Main kernels
volume = (Q['xyp'] <= Q['xyp'] + db.kDivM['xl'] * I['lyq'] * db.star[0]['qp'] + db.kDivMT['my'] * I['xmq'] * db.star[1]['qp'])
g.add('volume', volume)

def flux(dim,side1,side2):
  if dim == 0:
    return Q['xyp'] <= Q['xyp'] + db.FDivM[side1,side2]['xl'] * I['lyq'] * fluxSolver['qp']
  return Q['xyp'] <= Q['xyp'] + db.FDivMT[side1,side2]['my'] * I['xmq'] * fluxSolver['qp']
def fluxPrefetch(dim,side1,side2):
  if side1 == side2:
    if dim == 1:
      return Q if side1 == 1 else I
  elif side1 != side2:
    if dim != 1 or side1 != 1:
      return I
    else:
      return Q
  return None
g.addFamily('flux', simpleParameterSpace(2,2,2), flux, fluxPrefetch)

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