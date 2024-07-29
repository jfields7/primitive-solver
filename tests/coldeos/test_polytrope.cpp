//! \file test_polytrope.cpp
//  \brief Unit tests for the Polytrope EOS

#include <iostream>
#include <sstream>
#include <cmath>

#include "../../include/coldeos.hpp"
#include "../../include/ps_types.hpp"
#include "../../include/cold_idealgas.hpp"

#include "../include/testing.hpp"

using namespace Primitive;

bool TestConstruction(ColdEOS<Polytrope> * eos,
                      Real k,
                      Real gamma,
                      Real tol) {
  try {
    eos->SetGamma(gamma);
    eos->SetK(k);
  } catch (std::exception & e) {
    std::cout << "  Failed to initialize the EOS" << std::endl;
    std::cout << e.what() << std::endl;
    return false;
  }

  Real K = eos->GetK();
  Real err = GetError(k, K);
  if (err > tol) {
    std::cout << "  Could not read the polytropic constant" << std::endl;
    PrintError(k, K);
    return false;
  }

  return true;
}

bool TestPressure(ColdEOS<Polytrope>* peos,
                  Real p_want, Real rho, Real tol) {
  Real p = peos->GetPressure(rho);
  Real err = GetError(p_want, p);
  if (err > tol) {
    std::cout << "  Pressure not interpolated correctly" << std::endl;
    PrintError(p_want, p);
    return false;
  }
  return true;
}

bool TestEnergy(ColdEOS<Polytrope>* peos,
                  Real e_want, Real rho, Real tol) {
  Real e = peos->GetEnergy(rho);
  Real err = GetError(e_want, e);
  if (err > tol) {
    std::cout << "  Energy not interpolated correctly" << std::endl;
    PrintError(e_want, e);
    return false;
  }
  return true;
}

bool TestSpecificInternalEnergy(ColdEOS<Polytrope>* peos,
                                Real eps_want, Real rho, Real tol) {
  Real eps = peos->GetSpecificInternalEnergy(rho);
  Real err = GetError(eps_want, eps);
  if (err > tol) {
    std::cout << "  Specific internal energy not interpolated correctly" << std::endl;
    PrintError(eps_want, eps);
    return false;
  }
  return true;
}

bool TestEnthalpy(ColdEOS<Polytrope>* peos,
                  Real h_want, Real rho, Real tol) {
  Real h = peos->GetEnthalpy(rho);
  Real err = GetError(h_want, h);
  if (err > tol) {
    std::cout << "  Enthalpy not interpolated correctly" << std::endl;
    PrintError(h_want, h);
    return false;
  }
  return true;
}

int main(int argc, char ** argv) {
  UnitTests tester("ColdEOSCompOSE Tabulated Cold slice");
  ColdEOS<Polytrope> eos;

  Real const tol = 1e-12;
  Real gamma = 5.0/3.0;
  Real K = 100;

  tester.RunTest(&TestConstruction, "Initialize EOS", &eos, K, gamma, tol);

  Real rhop = 1e-5;
  Real p_want = K * std::pow(rhop, gamma);
  Real e_want = rhop + K / (gamma - 1.0) * std::pow(rhop, gamma);
  Real eps_want = e_want / rhop - 1;
  Real h_want = 1.0 + eps_want + p_want / rhop;

  tester.RunTest(&TestPressure, "Evaluate Pressure",
                 &eos, p_want, rhop, tol);

  tester.RunTest(&TestEnergy, "Evaluate Energy",
                 &eos, e_want, rhop, tol);

  tester.RunTest(&TestSpecificInternalEnergy, "Evaluate Specific Internal Energy",
                  &eos, eps_want, rhop, tol);

  tester.RunTest(&TestEnthalpy, "Evaluate Enthalpy",
                  &eos, h_want, rhop, tol);

  tester.PrintSummary();

  return 0;
}
