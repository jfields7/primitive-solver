#ifndef TEST_FUNCTIONS_HPP
#define TEST_FUNCTIONS_HPP

#include <eos.hpp>
#include <ps_types.hpp>
#include <testing.hpp>

/// Check that the temperature and energy density equations are consistent.
template<typename EOSPolicy, typename ErrorPolicy>
bool TestTemperatureFromEnergy(Primitive::EOS<EOSPolicy, ErrorPolicy>* eos, 
			Real n, Real T, Real *Y, const Real tol) {
  Real e = eos->GetEnergy(n, T, Y);
  Real T2 = eos->GetTemperature(n, e, Y);

  Real err = GetError(T, T2);
  if (err > tol) {
    PrintError(T, T2);
    return false;
  }

  return true;
}

/// Check that the temperature and pressure equations are consistent.
template<typename EOSPolicy, typename ErrorPolicy>
bool TestTemperatureFromPressure(Primitive::EOS<EOSPolicy, ErrorPolicy>* eos, 
			Real n, Real T, Real *Y, const Real tol) {
  Real p = eos->GetPressure(n, T, Y);
  Real T2 = eos->GetTemperatureFromP(n, p, Y);

  Real err = GetError(T, T2);
  if (err > tol) {
    PrintError(T, T2);
    return false;
  }

  return true;
}

/// Check that the enthalpy is consistent with the enthalpy calculated
/// directly from pressure and energy density.
template<typename EOSPolicy, typename ErrorPolicy>
bool TestEnthalpy(Primitive::EOS<EOSPolicy, ErrorPolicy>* eos, Real n, Real T, Real *Y, const Real tol) {
  Real h = eos->GetEnthalpy(n, T, Y);
  Real p = eos->GetPressure(n, T, Y);
  Real e = eos->GetEnergy(n, T, Y);

  Real expected = (e + p)/n;
  Real err = GetError(expected, h);
  if (err > tol) {
    PrintError(expected, h);
    return false;
  }

  return true;
}

/// Check that the specific energy is consistent with the energy density.
template<typename EOSPolicy, typename ErrorPolicy>
bool TestSpecificEnergy(Primitive::EOS<EOSPolicy, ErrorPolicy>* eos, Real n, Real T, Real *Y, const Real tol) {
  Real eps = eos->GetSpecificEnergy(n, T, Y);
  Real e = eos->GetEnergy(n, T, Y);
  Real mb = eos->GetBaryonMass();

  Real expected = (e/(n*mb) - 1.0);

  Real err = GetError(expected, eps);
  if (err > tol) {
    PrintError(expected, eps);
    return false;
  }

  return true;
}

#endif
