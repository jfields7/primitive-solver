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
