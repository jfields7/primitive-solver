//! \file test_idealgas.cpp
//  \brief Unit test for IdealGas EOS.

#include <iostream>
#include <cmath>

#include <eos.hpp>
#include <ps_types.hpp>
#include <idealgas.hpp>
#include <do_nothing.hpp>

#include <testing.hpp>
#include <testfunctions.hpp>

using namespace Primitive;

// Functions specific to IdealGas
bool TestConstruction() {
  EOS<IdealGas, DoNothing> eos;
  Real gamma = 5.0/3.0;
  Real mb = 1.0;
  return (gamma == eos.GetGamma() && mb == eos.GetBaryonMass());
}

bool TestSoundSpeed(EOS<IdealGas, DoNothing>* eos, Real n, Real *Y, Real tol) {
  bool success = true;
  Real gamma = eos->GetGamma();
  Real gammam1 = gamma - 1.0;
  
  for (Real T = 0; T < 1000.0; T += 50.0) {
    Real eps = eos->GetSpecificInternalEnergy(n, T, Y);
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

  EOS<IdealGas, DoNothing> eos;
  const Real tol = 1e-12; // error tolerance for floating-point arithmetic.
  Real n = 1.345e2; // number density
  Real T = 4.985e3; // temperature
  Real *Y = nullptr; // crap

  // Check that we can get the temperature back from the energy density
  // in a consistent way.
  tester.RunTest(&TestTemperatureFromEnergy<IdealGas, DoNothing>,
                 "Temperature from Energy Test", 
                 &eos, n, T, Y, tol);

  // Now check the same thing, but using pressure instead.
  tester.RunTest(&TestTemperatureFromPressure<IdealGas, DoNothing>,
                 "Temperature from Pressure Test",
                 &eos, n, T, Y, tol);

  // Check that the enthalpy is correct and consistent with its calculation
  // via pressure and density.
  tester.RunTest(&TestEnthalpy<IdealGas, DoNothing>, "Enthalpy Test",
                 &eos, n, T, Y, tol);

  // Check that the sound speed is consistent at a variety of temperatures with
  // the speed calculated via the energy per baryon.
  tester.RunTest(&TestSoundSpeed, "Sound Speed Test",
                 &eos, n, Y, tol);

  // Check that the specific energy is consistent with a direct calculation
  // based on the energy density.
  tester.RunTest(&TestSpecificInternalEnergy<IdealGas, DoNothing>, "Specific Energy Test",
                 &eos, n, T, Y, tol);

  // A few quantities in the ideal gas depend on the baryon mass. We need to make sure
  // that these are treated consistently if the baryon mass is not one.
  eos.SetBaryonMass(1.5);

  // Energy density
  tester.RunTest(&TestTemperatureFromEnergy<IdealGas, DoNothing>,
                 "Baryon Mass Consistency Test -- Energy",
                 &eos, n, T, Y, tol);

  // Enthalpy
  tester.RunTest(&TestEnthalpy<IdealGas, DoNothing>,
                 "Baryon Mass Consistency Test -- Enthalpy",
                 &eos, n, T, Y, tol);

  // Specific Energy
  tester.RunTest(&TestSpecificInternalEnergy<IdealGas, DoNothing>,
                 "Baryon Mass Consistency Test -- Specific Energy",
                 &eos, n, T, Y, tol);

  // Sound Speed
  tester.RunTest(&TestSoundSpeed,
                 "Baryon Mass Consistency Test -- Sound Speed",
                 &eos, n, Y, tol);

  tester.PrintSummary();

  return 0;
}
