//! \file test_pwp.cpp
//  \brief Unit tests for the ColdPiecewisePolytrope EOS

#include <iostream>
#include <sstream>
#include <cmath>

#include "../../include/coldeos.hpp"
#include "../../include/ps_types.hpp"
#include "../../include/cold_piecewise_polytrope.hpp"

#include "../include/testing.hpp"

using namespace Primitive;

bool TestConstruction(ColdEOS<ColdPiecewisePolytrope> * eos, Real tol) {
  try {
    const int N = 3;
    Real gamma_pieces[N] = {2.0, 1.7, 1.4};
    Real mb_nuc = 1.0;
    //Real rho_nuc = 1.0;
    Real density_pieces[N] = {0.0, 1.0, 3.0};
    Real P0 = 10.0;

    eos->InitializeFromData(density_pieces, gamma_pieces, P0, mb_nuc, N);
    eos->SetCodeUnitSystem(&Nuclear);
    eos->SetTemperature(0.01);
    eos->SetThermalGamma(5.0/3.0);
  } catch (std::exception & e) {
    std::cout << "  Failed to initialize the EOS" << std::endl;
    std::cout << e.what() << std::endl;
    return false;
  }

  return true;
}

bool TestPressure(ColdEOS<ColdPiecewisePolytrope>* peos,
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

bool TestEnergy(ColdEOS<ColdPiecewisePolytrope>* peos,
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

bool TestSpecificInternalEnergy(ColdEOS<ColdPiecewisePolytrope>* peos,
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

bool TestEnthalpy(ColdEOS<ColdPiecewisePolytrope>* peos,
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
  ColdEOS<ColdPiecewisePolytrope> eos;

  Real const tol = 1e-12;

  tester.RunTest(&TestConstruction, "Initialize EOS", &eos, tol);

  Real rhop = 0.5;
  Real gamma = 2.0;
  Real K = 10;
  Real T = 0.01;

  Real p_want = K * std::pow(rhop, gamma) + rhop*T;
  Real e_want = rhop + K / (gamma - 1.0) * std::pow(rhop, gamma) + rhop*T/(5.0/3.0 - 1.0);
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
