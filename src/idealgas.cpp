#include <idealgas.hpp>
#include <eos_units.hpp>
#include <cmath>
#include <stdexcept>

//! \file idealgas.cpp
//  \brief Implementation of IdealGas.

using namespace EOSUnits;

/// Constructor
IdealGas::IdealGas() {
  gamma = 5.0/3.0;
  mb = 1.0;
}

Real IdealGas::Temperature(Real n, Real e, Real *Y) {
  return gammam1*e/n;
}

Real IdealGas::Energy(Real n, Real T, Real *Y) {
  return n*T/gammam1;
}

Real IdealGas::Pressure(Real n, Real T, Real *Y) {
  return n*T;
}

Real IdealGas::Entropy(Real n, Real T, Real *Y) {
  throw std::logic_error("IdealGas: Entropy not currently implemented.");
}

Real IdealGas::Enthalpy(Real n, Real T, Real *Y) {
  Real e = Energy(n, T, Y);
  Real p = Pressure(n, T, Y);
  return (e + p)/n;
}

Real IdealGas::SoundSpeed(Real n, Real T, Real *Y) {
  return sqrt(gamma*gammam1*T/(gammam1*mb + gamma*T));
}

Real IdealGas::SpecificEnergy(Real n, Real T, Real *Y) {
  Real e = Energy(n, T, Y);
  return e/n - 1.0;
}
