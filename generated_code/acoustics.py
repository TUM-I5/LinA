#!/usr/bin/env python3

import argparse
from yateto import *
from yateto.input import parseJSONMatrixFile
from yateto.ast.transformer import DeduceIndices, EquivalentSparsityPattern

cmdLineParser = argparse.ArgumentParser()
cmdLineParser.add_argument('--matricesDir')
cmdLineParser.add_argument('--outputDir')
cmdLineParser.add_argument('--arch')
cmdLineParser.add_argument('--order')
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
db = parseJSONMatrixFile('{}/matrices_{}.json'.format(cmdLineArgs.matricesDir, degree), clones, transpose=transpose)
db.update( parseJSONMatrixFile('{}/star.json'.format(cmdLineArgs.matricesDir)) )

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
g.addFamily('flux', simpleParameterSpace(2,2,2), flux)

power = Scalar('power')
derivatives = [dQ0]
g.add('derivativeTaylorExpansion(0)', I['xyp'] <= power * dQ0['xyp'])
for i in range(1,order):
  derivativeSum = db.kTDivM['xl'] * derivatives[-1]['lyq'] * db.star[0]['qp'] + db.kTDivMT['my'] * derivatives[-1]['xmq'] * db.star[1]['qp']
  derivativeSum = DeduceIndices( Q['xyp'].indices ).visit(derivativeSum)
  derivativeSum = EquivalentSparsityPattern().visit(derivativeSum)
  dQ = Tensor('dQ({})'.format(i), qShape, spp=derivativeSum.eqspp())
  g.add('derivative({})'.format(i), dQ['xyp'] <= derivativeSum)
  g.add('derivativeTaylorExpansion({})'.format(i), I['xyp'] <= I['xyp'] + power * dQ['xyp'])
  derivatives.append(dQ)


## Initialization kernels
fluxScale = Scalar('fluxScale')
computeFluxSolver = fluxSolver['ij'] <= fluxScale * T['ik'] * Apm['kl'] * TT['lj']
g.add('computeFluxSolver', computeFluxSolver)

quadrature = Q['xyp'] <= db.quadrature['xl'] * db.quadrature['ym'] * initialCond['lmp']
g.add('quadrature', quadrature)

# Generate code
g.generate(cmdLineArgs.outputDir, 'lina')
