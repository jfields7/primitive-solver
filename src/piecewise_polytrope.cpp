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

/// Destructor
PiecewisePolytrope::~PiecewisePolytrope() {
  if (initialized) {
    delete[] density_pieces;
    delete[] a_pieces;
    delete[] gamma_pieces;
  }
}

void PiecewisePolytrope::AllocateMemory() {
  if (initialized) {
    delete[] density_pieces;
    delete[] a_pieces;
    delete[] gamma_pieces;
  }
  density_pieces = new Real[n_pieces];
  a_pieces = new Real[n_pieces];
  gamma_pieces = new Real[n_pieces];
}

int PiecewisePolytrope::FindPiece(Real n) const {
  int polytrope = 0;
  for (int i = 1; i < n_pieces; i++) {
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
  
  return (gamma_pieces[p] - 1.0)*(e - mb*n*(1.0 + a_pieces[p]))/n;
}

Real PiecewisePolytrope::TemperatureFromP(Real n, Real p, Real *Y) {
  return p/n;
}


Real PiecewisePolytrope::Energy(Real n, Real T, Real *Y) {
  int p = FindPiece(n);

  return n*(mb*(1.0 + a_pieces[p]) + T/(gamma_pieces[p] - 1.0));
}

Real PiecewisePolytrope::Pressure(Real n, Real T, Real *Y) {
  return n*T;
}

Real PiecewisePolytrope::Entropy(Real n, Real T, Real *Y) {
  throw std::logic_error("PiecewisePolytrope::Entropy not currently implemented.");
}

Real PiecewisePolytrope::Enthalpy(Real n, Real T, Real *Y) {
  int p = FindPiece(n);
  return mb*(1.0 + a_pieces[p]) + gamma_pieces[p]/(gamma_pieces[p] - 1.0)*T;
}

Real PiecewisePolytrope::MinimumEnthalpy() {
  return mb*(1.0 + a_pieces[0]);
}

Real PiecewisePolytrope::SoundSpeed(Real n, Real T, Real *Y) {
  int p = FindPiece(n);
  Real gamma = gamma_pieces[p];
  Real gammam1 = gamma - 1.0;
  return std::sqrt(gammam1*gamma*T/(gammam1*mb*(1.0 + a_pieces[p]) + gamma*T));
}

Real PiecewisePolytrope::SpecificEnergy(Real n, Real T, Real *Y) {
  int p = FindPiece(n);
  return a_pieces[p] + T/(mb*(gamma_pieces[p] - 1.0));
}

bool PiecewisePolytrope::ReadParametersFromFile(std::string fname) {
  return false;
}

bool PiecewisePolytrope::InitializeFromData(Real *densities, 
      Real *kappas, Real *gammas, Real m, int n) {
  // Before we even try to construct anything, we need to make sure that
  // the densities are ordered properly.
  for (int i = 1; i < n; i++) {
    if(densities[i] < densities[i-1]) {
      // The densities must be ordered from smallest to largest.
      return false;
    }
  }
  n_pieces = n;
  mb = m;
  initialized = true;
  AllocateMemory();

  // Now we can construct the different pieces.
  for (int i = 0; i < n; i++) {
    density_pieces[i] = densities[i]/m;
    gamma_pieces[i] = gammas[i];
    // Because we've rewritten the EOS in terms of temperature, we don't need
    // kappa in its current form. However, we can use it to define the a
    // constants that show up in our equations.
  }
  return true;
}
