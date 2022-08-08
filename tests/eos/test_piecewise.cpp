//! \file test_piecewise.cpp
//  \brief Unit tests for the PiecewisePolytrope EOS

#include <iostream>
#include <sstream>
#include <cmath>

#include <eos.hpp>
#include <ps_types.hpp>
#include <piecewise_polytrope.hpp>
#include <do_nothing.hpp>

#include <testing.hpp>
#include <testfunctions.hpp>

using namespace Primitive;

// Functions specific to PiecewisePolytrope
void Initialize(EOS<PiecewisePolytrope, DoNothing>& eos) {
  // This EOS is just a collection of goofy numbers
  // and in no way pertains to anything physical.
  const int N = 3;
  Real gamma_pieces[N] = {2.0, 1.7, 1.4};
  Real mb_nuc = 1.0;
  //Real rho_nuc = 1.0;
  Real density_pieces[N] = {1.0, 3.0, 10.0};
  Real rho_min = 0.1;
  Real kappa0 = 10.0;

  eos.InitializeFromData(density_pieces, gamma_pieces, rho_min, kappa0, mb_nuc, N);
}

void Reinitialize(EOS<PiecewisePolytrope, DoNothing>& eos) {
  const int N = 4;
  Real gamma_pieces[N] = {1.3333, 1.7, 1.9, 2.1};
  Real mb_nuc = 2.0;
  Real density_pieces[N] = {0.78, 1.4, 3.6, 6.4};
  Real rho_min = 1e-10;
  Real kappa0 = 3.0;

  eos.InitializeFromData(density_pieces, gamma_pieces, rho_min, kappa0, mb_nuc, N);
}

bool TestReinitialization() {
  EOS<PiecewisePolytrope, DoNothing> eos;
  Initialize(eos);

  // We should have data in here, but the wrong data.
  // Let's reinitialize with a different set of data.
  Reinitialize(eos);
  const int N = 4;
  Real gamma_pieces[N] = {1.3333, 1.7, 1.9, 2.1};
  Real mb_nuc = 2.0;
  
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

  // Test a range of temperatures
  int p = peos->FindPiece(n);
  Real rho = n*peos->GetBaryonMass();
  Real P_cold = peos->GetColdPressure(n, p) *
                peos->GetEOSUnitSystem()->PressureConversion(*peos->GetCodeUnitSystem());
  Real e_cold = peos->GetColdEnergy(n, p) *
                peos->GetEOSUnitSystem()->PressureConversion(*peos->GetCodeUnitSystem());
  Real h_cold = (e_cold + P_cold)/rho;
  Real csq_cold = peos->GetGamma(n)*P_cold/(e_cold + P_cold);
  for (Real T = 0; T < 1000.0; T += 50.0) {
    Real h = peos->GetEnthalpy(n, T, Y);
    Real P_th = peos->GetPressure(n, T, Y) - P_cold;
    Real e_th = peos->GetEnergy(n, T, Y) - e_cold;
    Real h_th = h - h_cold;
    Real csq_th = peos->GetThermalGamma()*P_th/(e_th + P_th);

    Real expected = sqrt((h_th*csq_th + h_cold*csq_cold)/h);
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
  tester.RunTest(&TestSpecificInternalEnergy<PiecewisePolytrope, DoNothing>,
                 ss.str(),
                 peos, n, T, Y, tol);
  ss.str(std::string());
}

bool TestContinuity(EOS<PiecewisePolytrope, DoNothing>* peos, Real n, Real T, Real* Y, Real tol) {
  // The energy density should be continuous, and to first order its dependence on the change
  // in the number density is linear. Therefore, we can estimate the slope at some points 
  // both above and below the transition density and predict exactly what it should be. If the 
  // two estimates for the energy density are not in reasonable agreement, there's probably a 
  // discontinuity.

  // Get the densities just above and just below the transition density.
  Real np1 = n*(1.0 + tol);
  Real np2 = n*(1.0 + 2.0*tol);
  Real nm1 = n*(1.0 - tol);
  Real nm2 = n*(1.0 - 2.0*tol);

  // Get two energy densities just beyond the transition density and
  // just below the transition density.
  Real ep1 = peos->GetEnergy(np1, T, Y);
  Real ep2 = peos->GetEnergy(np2, T, Y);
  Real em1 = peos->GetEnergy(nm1, T, Y);
  Real em2 = peos->GetEnergy(nm2, T, Y);
  
  // Get estimates for the transition density's energy density.
  Real ep = 2.0*ep1 - ep2;
  Real em = 2.0*em1 - em2;

  // The actual energy density at the transition density.
  Real ea = peos->GetEnergy(n, T, Y);

  // Now estimate the errors
  Real error = GetError(ep, em);
  
  if (error > tol) {
    std::cout << "  The energy density does not seem to be continuous:\n";
    std::cout << "  e+ = " << ep << "\n";
    std::cout << "  e- = " << em << "\n";
    std::cout << "  Error: " << error << "\n";
    std::cout << "  Actual energy density: " << ea << "\n";
    return false;
  }

  return true;
}

bool TestLowDensity() {
  EOS<PiecewisePolytrope, DoNothing> eos;
  const int N = 3;
  Real gamma_pieces[N] = {1.8, 2.3, 1.9};
  Real density_pieces[N] = {1.0, 1.5, 3.0};
  Real mb_nuc = 1.0;
  //Real rho_nuc = 1.0;
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
  Real T = 5.0;
  //Real T = eos.GetTemperatureFromE(n, 0, nullptr); // temperature
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

  // Test the transition densities.
  n = 0.1;
  tester.RunTest(&TestContinuity, "Default Transition Continuity Test", &eos, n, T, Y, tol);

  n = 1.0;
  tester.RunTest(&TestContinuity, "1st and 2nd Polytrope Continuity Test", &eos, n, T, Y, tol);

  n = 3.0;
  tester.RunTest(&TestContinuity, "2nd and 3rd Polytrope Continuity Test", &eos, n, T, Y, tol);

  // Consistency test for different baryon mass.
  Reinitialize(eos);
  // 1st polytrope tests
  n = 0.5;
  std::cout << "Baryon Mass Consistency Tests\n";
  RunTestSuite(tester, &eos, n, T, Y, 0, tol);
  n = 1.0;
  RunTestSuite(tester, &eos, n, T, Y, 1, tol);
  n = 2.0;
  RunTestSuite(tester, &eos, n, T, Y, 2, tol);
  n = 5.0;
  RunTestSuite(tester, &eos, n, T, Y, 3, tol);

  tester.PrintSummary();
  return 0;
}
