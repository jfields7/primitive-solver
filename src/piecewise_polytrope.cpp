//! \file piecewise_polytrope.cpp
//  \brief Implementation of PiecewisePolytrope EOSPolicy

#include <cmath>
#include <stdexcept>
#include <iostream>

#include <piecewise_polytrope.hpp>
#include <eos_units.hpp>

using namespace EOSUnits;
using namespace Primitive;

/// Constructor
PiecewisePolytrope::PiecewisePolytrope() {
  n_pieces = 0;
  initialized = false;
  n_species = 0;
}

/// Destructor
PiecewisePolytrope::~PiecewisePolytrope() {
  if (initialized) {
    delete[] density_pieces;
    delete[] a_pieces;
    delete[] gamma_pieces;
    delete[] kappa_pieces;
  }
}

void PiecewisePolytrope::AllocateMemory() {
  if (initialized) {
    delete[] density_pieces;
    delete[] a_pieces;
    delete[] gamma_pieces;
    delete[] kappa_pieces;
  }
  // We store one extra polytrope just in case a
  // minimum isn't specified.
  density_pieces = new Real[n_pieces+1];
  a_pieces = new Real[n_pieces+1];
  gamma_pieces = new Real[n_pieces+1];
  kappa_pieces = new Real[n_pieces+1];
}

int PiecewisePolytrope::FindPiece(Real n) const {
  // In case the density is below the minimum, we
  // implement a default case that is stored just
  // past the current polytrope.
  int polytrope = n_pieces-1;
  if (n < density_pieces[n_pieces]) {
    return n_pieces;
  }
  for (int i = 0; i < n_pieces; i++) {
    if (n <= density_pieces[i]) {
      polytrope = i;
      break;
    }
  }

  return polytrope;
}

Real PiecewisePolytrope::TemperatureFromE(Real n, Real e, Real *Y) {
  int p = FindPiece(n);
  
  //return (gamma_pieces[p] - 1.0)*(e - mb*n*(1.0 + a_pieces[p]))/n;
  return mb*kappa_pieces[p]*std::pow(mb*n, gamma_pieces[p] - 1.0);
}

Real PiecewisePolytrope::TemperatureFromP(Real n, Real p, Real *Y) {
  int i = FindPiece(n);
  return mb*kappa_pieces[i]*std::pow(mb*n, gamma_pieces[i] - 1.0);
  //return p/n;
}


Real PiecewisePolytrope::Energy(Real n, Real T, Real *Y) {
  int p = FindPiece(n);

  return n*(mb*(1.0 + a_pieces[p]) + T/(gamma_pieces[p] - 1.0));
  //Real rho = n*mb;
  //return rho*(1.0 + a_pieces[p]) + kappa_pieces[p]*std::pow(rho,gamma_pieces[p])/(gamma_pieces[p] - 1.0);
}

Real PiecewisePolytrope::Pressure(Real n, Real T, Real *Y) {
  //int p = FindPiece(n);
  return n*T;
  //return kappa_pieces[p]*std::pow(mb*n,gamma_pieces[p]);
}

Real PiecewisePolytrope::Entropy(Real n, Real T, Real *Y) {
  throw std::logic_error("PiecewisePolytrope::Entropy not currently implemented.");
}

Real PiecewisePolytrope::Enthalpy(Real n, Real T, Real *Y) {
  int p = FindPiece(n);
  return mb*(1.0 + a_pieces[p]) + gamma_pieces[p]/(gamma_pieces[p] - 1.0)*T;
  //return mb*(1.0 + a_pieces[p]) + mb*gamma_pieces[p]/(gamma_pieces[p] - 1.0)*kappa_pieces[p]*std::pow(mb*n, gamma_pieces[p]-1.0);
}

Real PiecewisePolytrope::MinimumEnthalpy() {
  return mb;
}

Real PiecewisePolytrope::SoundSpeed(Real n, Real T, Real *Y) {
  int p = FindPiece(n);
  Real gamma = gamma_pieces[p];
  Real gammam1 = gamma - 1.0;
  return std::sqrt(gammam1*gamma*T/(gammam1*mb*(1.0 + a_pieces[p]) + gamma*T));
  //Real P = Pressure(n, T, Y);
  //Real e = Energy(n, T, Y);
  //return std::sqrt(gamma*P/(e + P));
}

Real PiecewisePolytrope::SpecificEnergy(Real n, Real T, Real *Y) {
  int p = FindPiece(n);
  return a_pieces[p] + T/(mb*(gamma_pieces[p] - 1.0));
  //return a_pieces[p] + kappa_pieces[p]/(gamma_pieces[p] - 1.0)*std::pow(mb*n, gamma_pieces[p] - 1.0);
}

bool PiecewisePolytrope::ReadParametersFromFile(std::string fname) {
  return false;
}

bool PiecewisePolytrope::InitializeFromData(Real *densities, 
      Real *gammas, Real rho_min, Real P0, Real m, int n) {
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
  min_n = 0.0;
  max_n = densities[n-1]/mb;
  AllocateMemory();

  // Now we can construct the different pieces.
  density_pieces[0] = densities[0]/m;
  gamma_pieces[0] = gammas[0];
  //Real T = P0/(densities[0]/m);
  kappa_pieces[0] = P0/std::pow(densities[0],gammas[0]);
  if (n > 1){
    //a_pieces[0] = (T/mb)*(1.0/(gammas[0]-1.0) - 1.0/(gammas[1] - 1.0));
    a_pieces[0] = P0/densities[0]*(1.0/(gammas[0] - 1.0) - 1.0/(gammas[1] - 1.0));
  }
  else {
    a_pieces[0] = 0.0;
  }
  for (int i = 1; i < n; i++) {
    density_pieces[i] = densities[i]/m;
    gamma_pieces[i] = gammas[i];
    kappa_pieces[i] = kappa_pieces[i-1]*std::pow(densities[i-1],gammas[i-1]-gammas[i]);
    // Because we've rewritten the EOS in terms of temperature, we don't need
    // kappa in its current form. However, we can use it to define the a
    // constants that show up in our equations.
    //a_pieces[i] = 1.0 + a_pieces[i-1] + (T/mb)*(1.0/(gammas[i-1] - 1.0) - 1.0/(gammas[i] - 1.0));
    a_pieces[i] = a_pieces[i-1] + kappa_pieces[i-1]/(gammas[i-1] - 1.0)*std::pow(densities[i-1],gammas[i-1]-1.0)
                  - kappa_pieces[i]/(gammas[i] - 1.0)*std::pow(densities[i-1],gammas[i]-1.0);
    // Let's double-check that the density is physical.
    /*if (gamma_pieces[i] > 2.0) {
      //Real rho_max = std::pow((a_pieces[i] + 1.0)*std::pow(densities[i-1],gammas[i])
      //                *(gammas[i]-1.0)/(P*gammas[i]*(gammas[i]-2.0)), 1.0/(gammas[i]-1.0));
      Real rho_max = std::pow(mb*(a_pieces[i] + 1.0)*std::pow(densities[i-1],gammas[i]-1.0)
                        *(gammas[i]-1.0)/(T*gammas[i]*(gammas[i]-2.0)), 1.0/(gammas[i]-1.0));
      if (densities[i] > rho_max) {
        std::cout << "The i = " << i 
                  << " piece of the polytrope permits superluminal densities: \n";
        std::cout << "  rho[i]     = " << densities[i] << "\n";
        std::cout << "  rho_max[i] = " << rho_max << "\n";
        return false;
      }
    }*/
    //T = T*std::pow(densities[i]/densities[i-1],gammas[i]-1.0);
  }

  // Set up the default case. Because of thermodynamic constraints,
  // we must force a to be zero. Gamma is fixed by continuity.
  a_pieces[n] = 0.0;
  Real factor = (gammas[0] - 1.0)*rho_min*a_pieces[0];
  Real P_min = P0*std::pow(rho_min/densities[0],gammas[0]);
  gamma_pieces[n] = (gammas[0]*P_min + factor)/(P_min + factor);
  density_pieces[n] = rho_min/mb;
  kappa_pieces[n] = P_min/std::pow(rho_min, gamma_pieces[n]);

  // Find the maximum energies allowed by the EOS. Because of the
  // extension down to zero density, energy is just capped at 0.
  //Real T_max = T;
  Real T_max = kappa_pieces[n-1]*std::pow(densities[n-1],gammas[n-1]-1.0)/mb;
  min_e = 0.0;
  max_e = Energy(max_n, T_max, nullptr);

  // DEBUG ONLY:
  for (int i = 0; i <= n; i++) {
    std::cout << "Polytrope: i = " << i << "\n";
    std::cout << "  n = " << density_pieces[i] << "\n";
    std::cout << "  gamma = " << gamma_pieces[i] << "\n";
    std::cout << "  a = " << a_pieces[i] << "\n";
    std::cout << "  kappa = " << kappa_pieces[i] << "\n";
  }
  std::cout << "  P0 = " << P0 << "\n";

  initialized = true;
  return true;
}
