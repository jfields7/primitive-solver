//! \file test_piecewise.cpp
//  \brief Unit tests for the PiecewisePolytrope EOS

#include <iostream>
#include <sstream>
#include <cmath>

#include <eos.hpp>
#include <ps_types.hpp>
#include <piecewise_polytrope.hpp>
#include <do_nothing.hpp>
#include <eos_units.hpp>

#include <testing.hpp>
#include <testfunctions.hpp>

using namespace Primitive;

// Functions specific to PiecewisePolytrope
void Initialize(EOS<PiecewisePolytrope, DoNothing>& eos) {
  // This EOS is just a collection of goofy numbers
  // and in no way pertains to anything physical.
  int N = 3;
  Real gamma_pieces[N] = {2.0, 1.7, 1.4};
  Real mb_nuc = 1.0;
  Real rho_nuc = 1.0;
  Real density_pieces[N] = {1.0, 3.0, 10.0};
  Real rho_min = 0.1;
  Real kappa0 = 10.0;

  eos.InitializeFromData(density_pieces, gamma_pieces, rho_min, kappa0, mb_nuc, N);
}

bool TestConstruction() {
  EOS<PiecewisePolytrope, DoNothing> eos;
  if (eos.IsInitialized() || eos.GetNPieces() > 0) {
    return false;
  }
  Initialize(eos);
  if (!eos.IsInitialized() || eos.GetNPieces() != 3) {
    return false;
  }
  // Check that all the pieces were constructed correctly.
  if (eos.GetGamma(0.5) != 2.0) {
    return false;
  }
  else if (eos.GetGamma(1.5) != 1.7) {
    return false;
  }
  else if (eos.GetGamma(6.0) != 1.4) {
    return false;
  }

  // Check that the min and max densities are correct.
  if (eos.GetMinimumDensity() != 0.1) {
    std::cout << "  Incorrect minimum density.\n";
    return false;
  }
  if (eos.GetMaximumDensity() != 10.0) {
    std::cout << "  Incorrect maximum density.\n";
    return false;
  }
  return true;
}

bool TestSoundSpeed(EOS<PiecewisePolytrope, DoNothing>* peos, Real n,
    Real *Y, Real tol) {
  bool success = true;
  Real gamma = peos->GetGamma(n);

  // Test a range of temperatures
  for (Real T = 0; T < 1000.0; T += 50.0) {
    Real p = peos->GetPressure(n, T, Y);
    Real e = peos->GetEnergy(n, T, Y);
    Real expected = sqrt(gamma*p/(e + p));
    Real cs = peos->GetSoundSpeed(n, T, Y);

    Real err = GetError(expected, cs);

    if (err > tol) {
      success = false;
      std::cout << "  Temperature: " << T << "\n";
      PrintError(expected, cs);
    }
  }

  return success;
}

void RunTestSuite(UnitTests& tester, EOS<PiecewisePolytrope, DoNothing>* peos,
      Real n, Real T, Real* Y, int polytrope, Real tol) {
  std::stringstream ss;

  ss << "Temperature From Energy Test - Polytrope " << polytrope;
  tester.RunTest(&TestTemperatureFromEnergy<PiecewisePolytrope, DoNothing>,
                 ss.str(),
                 peos, n, T, Y, tol);
  ss.str(std::string());

  ss << "Temperature From Pressure Test - Polytrope " << polytrope;
  tester.RunTest(&TestTemperatureFromPressure<PiecewisePolytrope, DoNothing>,
                 ss.str(),
                 peos, n, T, Y, tol);
  ss.str(std::string());

  ss << "Enthalpy Test - Polytrope " << polytrope;
  tester.RunTest(&TestEnthalpy<PiecewisePolytrope, DoNothing>,
                 ss.str(),
                 peos, n, T, Y, tol);
  ss.str(std::string());

  ss << "Sound Speed Test - Polytrope " << polytrope;
  tester.RunTest(&TestSoundSpeed, ss.str(), peos, n, Y, tol);
  ss.str(std::string());

  ss << "Specific Energy Test - Polytrope " << polytrope;
  tester.RunTest(&TestSpecificEnergy<PiecewisePolytrope, DoNothing>,
                 ss.str(),
                 peos, n, T, Y, tol);
  ss.str(std::string());
}

int main(int argc, char* argv[]) {
  UnitTests tester("Piecewise Polytropic EOS");
  // Validate that the EOS gets constructed as expected.
  tester.RunTest(&TestConstruction, "Construction Test");

  EOS<PiecewisePolytrope, DoNothing> eos;
  Initialize(eos);
  const Real tol = 1e-12;

  // 1st polytrope tests
  Real n = 0.5; // number density
  Real T = 100.0; // temperature
  Real *Y = nullptr; // crap
  RunTestSuite(tester, &eos, n, T, Y, 0, tol);

  // 2nd polytrope tests
  n = 1.5; // number density
  RunTestSuite(tester, &eos, n, T, Y, 1, tol);

  // 3rd polytrope tests
  n = 6.0;
  RunTestSuite(tester, &eos, n, T, Y, 2, tol);

  tester.PrintSummary();
  return 0;
}
