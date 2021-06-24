//! \file test_idealgas.cpp
//  \brief Unit test for IdealGas EOS.

#include <iostream>
#include <cmath>

#include <eos.hpp>
#include <ps_types.hpp>
#include <idealgas.hpp>

#include <testing.hpp>
#include <testfunctions.hpp>

// Functions specific to IdealGas
bool TestConstruction() {
  EOS<IdealGas> eos;
  Real gamma = 5.0/3.0;
  Real mb = 1.0;
  return (gamma == eos.GetGamma() && mb == eos.GetBaryonMass());
}

bool TestSoundSpeed(EOS<IdealGas>* eos, Real n, Real *Y, Real tol) {
  bool success = true;
  Real gamma = eos->GetGamma();
  Real gammam1 = gamma - 1.0;
  
  for (Real T = 0; T < 1000.0; T += 50.0) {
    Real eps = eos->GetSpecificEnergy(n, T, Y);
    Real expected = sqrt(gammam1*eps/(eps + 1.0/gamma));
    Real cs = eos->GetSoundSpeed(n, T, Y);
    Real err = GetError(expected, cs);

    if (err > tol) {
      success = false;
      std::cout << "  Temperature: " << T << "\n";
      PrintError(expected, cs);
    }
  }

  return success;
}

int main(int argc, char *argv[]) {
  UnitTests tester{"Ideal Gas EOS"};
  // Validate that the gas was constructed as expected.
  tester.RunTest(&TestConstruction, "Construction Test");

  EOS<IdealGas> eos;
  const Real tol = 1e-12; // error tolerance for floating-point arithmetic.
  Real n = 1.345e2; // number density
  Real T = 4.985e3; // temperature
  Real *Y = nullptr; // crap

  // Check that we can get the temperature back from the energy density
  // in a consistent way.
  tester.RunTest(&TestTemperatureFromEnergy<IdealGas>,
                 "Temperature from Energy Test", 
                 &eos, n, T, Y, tol);

  // Now check the same thing, but using pressure instead.
  tester.RunTest(&TestTemperatureFromPressure<IdealGas>,
                 "Temperature from Pressure Test",
                 &eos, n, T, Y, tol);

  // Check that the enthalpy is correct and consistent with its calculation
  // via pressure and density.
  tester.RunTest(&TestEnthalpy<IdealGas>, "Enthalpy Test",
                 &eos, n, T, Y, tol);

  // Check that the sound speed is consistent at a variety of temperatures with
  // the speed calculated via the energy per baryon.
  tester.RunTest(&TestSoundSpeed, "Sound Speed Test",
                 &eos, n, Y, tol);

  tester.PrintSummary();

  return 0;
}
