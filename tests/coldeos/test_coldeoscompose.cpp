//! \file test_coldeoscompose.cpp
//  \brief Unit tests for the ColdEOSCompOSE EOS

#include <iostream>
#include <sstream>
#include <cmath>

#include "../../include/coldeos.hpp"
#include "../../include/ps_types.hpp"
#include "../../include/cold_eos_compose.hpp"

#include "../../include/eos.hpp"
#include "../../include/eos_compose.hpp"
#include "../../include/do_nothing.hpp"

#include "../include/testing.hpp"

using namespace Primitive;

bool TestConstruction(ColdEOS<ColdEOSCompOSE> * eos, Real tol) {
  try {
    std::string species_names[NSCALARS] = {"e"};
    eos->ReadColdSliceFromFile("../data/SFHo_full.h5", species_names);
    // eos->SetCodeUnitSystem(&GeometricSolar);
  } catch (std::exception & e) {
    std::cout << "  Failed to initialize the EOS" << std::endl;
    std::cout << e.what() << std::endl;
    return false;
  }

  Real mb = eos->GetBaryonMass();
  Real mb_want = 939.56535;
  Real err = GetError(mb_want, mb);
  if (err > tol) {
    std::cout << "  Could not read the baryon mass" << std::endl;
    PrintError(mb_want, mb);
    return false;
  }

  Real min_n = eos->GetMinimumDensity();
  Real min_n_want = 1e-12;
  err = GetError(min_n_want, min_n);
  if (err > tol) {
    std::cout << "  Wrong minimum density" << std::endl;
    PrintError(min_n_want, min_n);
    return false;
  }

  Real max_n = eos->GetMaximumDensity();
  Real max_n_want = 1.4454398;
  err = GetError(max_n_want, max_n);
  if (err > tol) {
    std::cout << "  Wrong maximum density" << std::endl;
    PrintError(max_n_want, max_n);
    return false;
  }

  Real const * table = eos->GetRawTable();
  Real log_p = table[eos->index(ColdEOSCompOSE::ECLOGP, 3)];
  Real log_p_want = -2.647301149233892e+01;
  err = GetError(log_p_want, log_p);
  if (err > tol) {
    std::cout << "  Data not read correctly" << std::endl;
    PrintError(log_p_want, log_p);
    return false;
  }

  return true;
}

bool TestPressure(ColdEOS<ColdEOSCompOSE>* peos,
                  EOS<EOSCompOSE,DoNothing>* pfull_eos,
                  Real rho, Real tol) {
  Real n = rho/pfull_eos->GetBaryonMass();
  Real ye = peos->GetY(rho, 0);
  Real T = peos->GetTemperature(); // * GeometricSolar.TemperatureConversion(Nuclear);
  printf("n = %e\n", n);
  printf("Ye = %e\n", ye);
  printf("T = %e\n", T);
  Real p_want = pfull_eos->GetPressure(n, T, &ye);
  Real p = peos->GetPressure(rho);
  Real err = GetError(p_want, p);
  if (err > tol) {
    std::cout << "  Pressure not interpolated correctly" << std::endl;
    PrintError(p_want, p);
    return false;
  }
  return true;
}

bool TestSpecificInternalEnergy(ColdEOS<ColdEOSCompOSE>* peos,
                                EOS<EOSCompOSE,DoNothing>* pfull_eos,
                                Real rho, Real tol) {
  Real n = rho/pfull_eos->GetBaryonMass();
  Real ye = peos->GetY(rho, 0);
  Real T = peos->GetTemperature(); // * GeometricSolar.TemperatureConversion(Nuclear);
  printf("n = %e\n", n);
  printf("Ye = %e\n", ye);
  printf("T = %e\n", T);
  Real eps_want = pfull_eos->GetSpecificInternalEnergy(n, T, &ye);
  Real eps = peos->GetSpecificInternalEnergy(rho);
  Real err = GetError(eps_want, eps);
  if (err > tol) {
    std::cout << "  Specific internal energy not interpolated correctly" << std::endl;
    PrintError(eps_want, eps);
    return false;
  }
  return true;
}

bool TestEnthalpy(ColdEOS<ColdEOSCompOSE>* peos,
                  EOS<EOSCompOSE,DoNothing>* pfull_eos,
                  Real rho, Real tol) {
  Real n = rho/peos->GetBaryonMass();
  Real ye = peos->GetY(rho, 0);
  Real T = peos->GetTemperature(); // * GeometricSolar.TemperatureConversion(Nuclear);
  printf("n = %e\n", n);
  printf("Ye = %e\n", ye);
  printf("T = %e\n", T);
  Real h_want = pfull_eos->GetEnthalpy(n, T, &ye);
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
  ColdEOS<ColdEOSCompOSE> eos;
  EOS<EOSCompOSE,DoNothing> full_eos;

  Real const tol = 1e-12;

  full_eos.ReadTableFromFile("../data/SFHo_full.h5");

  tester.RunTest(&TestConstruction, "Read HDF5 Table", &eos, tol);

  Real rhop = 0.16 * eos.GetBaryonMass();
  printf("rhop = %e\n", rhop);

  tester.RunTest(&TestPressure, "Evaluate Pressure", &eos, &full_eos, rhop, tol);

  tester.RunTest(&TestSpecificInternalEnergy, "Evaluate Specific Internal Energy",
                  &eos, &full_eos, rhop, tol);

  tester.RunTest(&TestEnthalpy, "Evaluate Enthalpy", &eos, &full_eos, rhop, tol);

  tester.PrintSummary();

  return 0;
}
