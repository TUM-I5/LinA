#ifndef BASISFUNCTIONS_H_
#define BASISFUNCTIONS_H_

#include "nodes.h"

namespace lina {
  static double basisFunction(double xi, int i) {
    double num = 1.0;
    double denom = 1.0;
    for (int j = 0; j < CONVERGENCE_ORDER; ++j) {
      if (i != j) {
        num *= xi - LGLNodes[j];
        denom *= LGLNodes[i] - LGLNodes[j];
      }
    }
    return num / denom;
  }
}

#endif // BASISFUNCTIONS_H_
