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
numberOfQuantities = 4

qShape = (numberOf1DBasisFunctions, numberOf1DBasisFunctions, numberOf1DBasisFunctions, numberOfQuantities)

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

initialCond = Tensor('initialCond', (numberOfQuadraturePoints, numberOfQuadraturePoints, numberOfQuadraturePoints, numberOfQuantities))

# Flux solver
fluxSolver = Tensor('fluxSolver', (numberOfQuantities, numberOfQuantities))
T = Tensor('T', (numberOfQuantities, numberOfQuantities))
TT = Tensor('TT', (numberOfQuantities, numberOfQuantities))
Apm = Tensor('Apm', (numberOfQuantities, numberOfQuantities))

# Kernels
g = Generator(arch)

## Main kernels
volume = (Q['xyzp'] <= Q['xyzp'] + db.kDivM['xl'] * I['lyzq'] * db.star[0]['qp'] + db.kDivMT['my'] * I['xmzq'] * db.star[1]['qp'] + db.kDivMT['nz'] * I['xynq'] * db.star[2]['qp'])
g.add('volume', volume)

def flux(dim,side1,side2):
  if dim == 0:
    return Q['xyzp'] <= Q['xyzp'] + db.FDivM[side1,side2]['xl'] * I['lyzq'] * fluxSolver['qp']
  elif dim == 1:
    return Q['xyzp'] <= Q['xyzp'] + db.FDivMT[side1,side2]['my'] * I['xmzq'] * fluxSolver['qp']
  return Q['xyzp'] <= Q['xyzp'] + db.FDivMT[side1,side2]['nz'] * I['xynq'] * fluxSolver['qp']
def fluxPrefetch(dim,side1,side2):
  if side1 == side2:
    if dim == 1:
      return Q if side1 == 1 else I
  elif side1 != side2:
    if dim != 1 or side1 != 1:
      return I
  return None
g.addFamily('flux', simpleParameterSpace(3,2,2), flux, fluxPrefetch)

power = Scalar('power')
derivatives = [dQ0]
g.add('derivativeTaylorExpansion(0)', I['xyzp'] <= power * dQ0['xyzp'])
for i in range(1,order):
  derivativeSum = db.kTDivM['xl'] * derivatives[-1]['lyzq'] * db.star[0]['qp']    \
                  + db.kTDivMT['my'] * derivatives[-1]['xmzq'] * db.star[1]['qp'] \
                  + db.kTDivMT['nz'] * derivatives[-1]['xynq'] * db.star[2]['qp']
  derivativeSum = DeduceIndices( Q['xyzp'].indices ).visit(derivativeSum)
  derivativeSum = EquivalentSparsityPattern().visit(derivativeSum)
  dQ = Tensor('dQ({})'.format(i), qShape, spp=derivativeSum.eqspp())
  g.add('derivative({})'.format(i), dQ['xyzp'] <= derivativeSum)
  g.add('derivativeTaylorExpansion({})'.format(i), I['xyzp'] <= I['xyzp'] + power * dQ['xyzp'])
  derivatives.append(dQ)


## Initialization kernels
fluxScale = Scalar('fluxScale')
computeFluxSolver = fluxSolver['ij'] <= fluxScale * T['ik'] * Apm['kl'] * TT['lj']
g.add('computeFluxSolver', computeFluxSolver)

quadrature = Q['xyzp'] <= db.quadrature['xl'] * db.quadrature['ym'] * db.quadrature['zn'] * initialCond['lmnp']
g.add('quadrature', quadrature)

# Generate code
gemmTool = DefaultGeneratorCollection(arch)
g.generate(cmdLineArgs.outputDir, 'lina', gemmTool)
