#!/usr/bin/env python3

import os
import numpy
import math

with open('nodes.cpp', 'w') as out:
  with open('nodes.c', 'r') as bf:
    out.write('/** This file is generated. Do not edit. */\n')
    out.write('#include "nodes.h"\n')
    out.write('namespace lina {\n')
    out.write('#if !defined(CONVERGENCE_ORDER)\n')
    out.write('#error CONVERGENCE_ORDER must be set.\n')
    for line in bf:
      nodes = line.rstrip().split(' ')
      out.write('#elif CONVERGENCE_ORDER == {}\n'.format(len(nodes)))
      out.write('double const LGLNodes[] = {{{}}};\n'.format(', '.join(nodes)))
    out.write('#endif\n')
    out.write('}\n')
