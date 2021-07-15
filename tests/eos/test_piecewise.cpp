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

bool TestReinitialization() {
  EOS<PiecewisePolytrope, DoNothing> eos;
  Initialize(eos);

  // We should have data in here, but the wrong data.
  // Let's reinitialize with a different set of data.
  int N = 4;
  Real gamma_pieces[N] = {1.3333, 1.7, 1.9, 2.1};
  Real mb_nuc = 2.0;
  Real rho_nuc = 1.0;
  Real density_pieces[N] = {0.78, 1.4, 3.6, 6.4};
  Real rho_min = 1e-10;
  Real kappa0 = 3.0;

  eos.InitializeFromData(density_pieces, gamma_pieces, rho_min, kappa0, mb_nuc, N);
  
  if (!eos.IsInitialized()) {
    std::cout << "  Failed to initialize EOS properly.\n";
    return false;
  }
  if (eos.GetNPieces() != 4) {
    std::cout << "  Wrong number of polytropes.\n";
    if (eos.GetNPieces() == 3){
      std::cout << "  n_pieces = 3 suggests a failed reinitialization.\n";
    }
    return false;
  }

  // Make sure that all the pieces were constructed correctly.
  if (eos.GetGamma(0.5/mb_nuc) != gamma_pieces[0]) {
    std::cout << "  First polytrope is incorrect.\n";
    return false;
  }
  if (eos.GetGamma(1.25/mb_nuc) != gamma_pieces[1]) {
    std::cout << "  Second polytrope is incorrect.\n";
    return false;
  }
  if (eos.GetGamma(2.5/mb_nuc) != gamma_pieces[2]) {
    std::cout << "  Third polytrope is incorrect.\n";
    return false;
  }
  if (eos.GetGamma(5.0/mb_nuc) != gamma_pieces[3]) {
    std::cout << "  Fourth polytrope is incorrect.\n";
    return false;
  }
  return true;
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
  if (eos.GetMinimumDensity() != 0.0) {
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

bool ContinuityTest(EOS<PiecewisePolytrope, DoNothing>* peos, Real n, Real T, Real* Y, Real tol) {
  Real np = n*(1.0 + tol);
  Real nm = n*(1.0 - tol);

  Real ep = peos->GetEnergy(np, T, Y);
  Real em = peos->GetEnergy(nm, T, Y);

  return true;
}

bool TestLowDensity() {
  EOS<PiecewisePolytrope, DoNothing> eos;
  int N = 3;
  Real gamma_pieces[N] = {1.8, 2.3, 1.9};
  Real density_pieces[N] = {1.0, 1.5, 3.0};
  Real mb_nuc = 1.0;
  Real rho_nuc = 1.0;
  Real rho_min = 0.5;
  Real P0 = 10.0;

  eos.InitializeFromData(density_pieces, gamma_pieces, rho_min, P0, mb_nuc, N);

  if (!eos.IsInitialized()) {
    std::cout << "  There was an error initializing the EOS.\n";
    return false;
  }

  Real min_gamma = eos.GetGamma(rho_min/2.0);
  for (int i = 0; i < N; i++) {
    if (min_gamma == gamma_pieces[i]) {
      std::cout << "  Wrong piece retrieved for low densities.\n";
      std::cout << "  Returning the i = " << i << " piece.\n";
      return false;
    }
  }

  return true;
}

int main(int argc, char* argv[]) {
  UnitTests tester("Piecewise Polytropic EOS");
  // Validate that the EOS gets constructed as expected.
  tester.RunTest(&TestConstruction, "Construction Test");

  // Validate that reinitialization works as expected.
  tester.RunTest(&TestReinitialization, "Reinitialization Test");

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

  // Low density test
  tester.RunTest(&TestLowDensity, "Low Density Test");

  tester.PrintSummary();
  return 0;
}
