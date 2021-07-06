//! \file piecewise_polytrope.cpp
//  \brief Implementation of PiecewisePolytrope EOSPolicy

#include <cmath>
#include <stdexcept>

#include <piecewise_polytrope.hpp>
#include <eos_units.hpp>

using namespace EOSUnits;
using namespace Primitive;

/// Constructor
PiecewisePolytrope::PiecewisePolytrope() {
  n_pieces = 0;
  initialized = false;
}

int PiecewisePolytrope::FindPiece(Real n) const {
  int polytrope = -1;
  for (int i = 0; i < n_pieces; i++) {
    if (n > density_pieces[i]) {
      polytrope = i;
    }
    else {
      break;
    }
  }
  if (n == -1) {
    // We need to do something if the density is invalid, 
    // but I'm not sure what it is.
  }

  return polytrope;
}

Real PiecewisePolytrope::TemperatureFromE(Real n, Real e, Real *Y) {
  int p = FindPiece(n);
  
  return 0;
}

Real PiecewisePolytrope::Energy(Real n, Real T, Real *Y) {
  int p = FindPiece(n);

  return a_pieces[p] + n*(mb + T/(gamma_pieces[p] - 1.0));
}
