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
  initialized = true;
}

int PiecewisePolytrope::FindPiece(Real n) const {
  int polytrope = -1;
  for (int i = 0; i < n_pieces; i++) {
    if (n <= density_pieces[i]) {
      polytrope = i;
      break;
    }
  }
  if (polytrope == -1) {
    // We need to do something if the density is invalid, 
    // but I'm not sure what it is. Maybe we'll just throw
    // an exception for now.
    throw std::out_of_range("Input density does not correspond to a valid polytrope.");
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
      Real *gammas, Real rho_min, Real kappa0, Real m, int n) {
  // Make sure that we actually *have* polytropes.
  if (n <= 0) {
    return false;
  }
  // Before we even try to construct anything, we need to make sure that
  // the densities are ordered properly.
  for (int i = 1; i < n; i++) {
    if(densities[i] <= densities[i-1]) {
      // The densities must be ordered from smallest to largest and strictly
      // increasing.
      return false;
    }
  }
  if (rho_min >= densities[0]) {
    return false;
  }

  // Initialize (most of) the member variables.
  n_pieces = n;
  mb = m;
  min_n = rho_min/mb;
  max_n = densities[n-1]/mb;
  AllocateMemory();

  // Now we can construct the different pieces.
  density_pieces[0] = densities[0]/m;
  gamma_pieces[0] = gammas[0];
  a_pieces[0] = 0.0; // TODO: Figure out what this actually is.
  Real kappa = kappa0;
  for (int i = 1; i < n; i++) {
    density_pieces[i] = densities[i]/m;
    gamma_pieces[i] = gammas[i];
    // Because we've rewritten the EOS in terms of temperature, we don't need
    // kappa in its current form. However, we can use it to define the a
    // constants that show up in our equations.
    a_pieces[i] = a_pieces[i-1] + kappa*std::pow(densities[i-1],gammas[i-1]-1.0)*
                  (1.0/(gammas[i-1] - 1.0) - 1.0/(gammas[i] - 1.0));
    kappa = kappa*std::pow(densities[i-1],gammas[i-1] - gammas[i]);
  }

  // Find the minimum and maximum energies allowed by the EOS.
  Real T_max = kappa/mb*std::pow(densities[n-1],gammas[n-1]-1.0);
  min_e = Energy(min_n, 0.0, nullptr);
  max_e = Energy(max_n, T_max, nullptr);
  return true;
}
