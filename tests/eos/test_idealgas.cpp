//! \file test_idealgas.cpp
//  \brief Unit test for IdealGas EOS.

#include <iostream>
#include <cmath>

#include <eos.hpp>
#include <ps_types.hpp>
#include <idealgas.hpp>

#include <testing.hpp>

int main(int argc, char *argv[]) {
  UnitTests tester{"Ideal Gas EOS"};
  // Validate that the gas was constructed as expected.
  tester.RunTest([]() -> bool {
    EOS<IdealGas> eos;
    Real gamma = 5.0/3.0;
    Real mb = 1.0;
    return (gamma == eos.GetGamma() && mb == eos.GetBaryonMass());
  }, "Construction Test");

  // Check that we can get the temperature back from the energy density
  // in a consistent way.
  tester.RunTest([]() -> bool {
    EOS<IdealGas> eos;
    const Real tol = 1e-12; // error tolerance for floating-point arithmetic.

    Real n = 1.345e2; // number density
    Real T = 4.985e3; // temperature
    Real *Y = nullptr; // crap

    Real e = eos.GetEnergy(n, T, Y);
    Real T2 = eos.GetTemperature(n, e, Y);

    Real err = (T - T2)/T;
    std::cout << "  Expected: " << T << "\n";
    std::cout << "  Actual: " << T2 << "\n";
    std::cout << "  Error: " << err << "\n";

    return fabs(err) <= tol;
  }, "Temperature from Energy Test");

  // Now check the same thing, but using pressure instead.
  tester.RunTest([]() -> bool {
    EOS<IdealGas> eos;
    const Real tol = 1e-12; // error tolerance for floating-point arithmetic.

    Real n = 1.345e2; // number density
    Real T = 4.985e3; // temperature
    Real *Y = nullptr; // crap

    Real p = eos.GetPressure(n, T, Y);
    Real T2 = eos.GetTemperatureFromP(n, p, Y);

    Real err = (T - T2)/T;
    std::cout << "  Expected: " << T << "\n";
    std::cout << "  Actual: " << T2 << "\n";
    std::cout << "  Error: " << err << "\n";

    return fabs(err) <= tol;
  }, "Temperature from Pressure Test");

  // Check that the enthalpy is correct and consistent with its calculation
  // via pressure and density.
  tester.RunTest([]() -> bool {
    EOS<IdealGas> eos;
    const Real tol = 1e-12; // error tolerance for floating-point arithmetic.

    Real n = 1.345e2; // number density
    Real T = 4.895e3; // temperature
    Real *Y = nullptr; // crap

    Real h = eos.GetEnthalpy(n, T, Y);
    Real p = eos.GetPressure(n, T, Y);
    Real e = eos.GetEnergy(n, T, Y);

    Real expected = (e + p)/n;
    Real err = (expected - h)/expected;
    std::cout << "  Expected: " << expected << "\n";
    std::cout << "  Actual: " << h << "\n";
    std::cout << "  Error: " << err << "\n";

    return fabs(err) <= tol;
  }, "Enthalpy Test");

  // Check that the sound speed is consistent at a variety of temperatures with
  // the speed calculated via the energy per baryon.
  tester.RunTest([]() -> bool {
    EOS<IdealGas> eos;
    const Real tol = 1e-12; // error tolerance for floating-point arithmetic.
    
    bool success = true;
    Real n = 1.345e2; // number density
    Real *Y = nullptr; // crap
    Real gamma = eos.GetGamma();
    Real gammam1 = gamma - 1.0;

    for(Real T = 0; T < 1000.0; T += 50.0) {
      Real eps = eos.GetSpecificEnergy(n, T, Y);
      Real expected = sqrt(gammam1*eps/(eps + 1.0/gamma));
      Real cs = eos.GetSoundSpeed(n, T, Y);
      Real err = (expected - cs)/expected;

      if (fabs(err) > tol) {
        success = false;
        std::cout << "  Temperature: " << T << "\n";
        std::cout << "  Expected: " << expected << "\n";
        std::cout << "  Actual: " << cs << "\n";
        std::cout << "  Error: " << err << "\n";
      }
    }

    return success;
  }, "Sound Speed Test");

  tester.PrintSummary();

  return 0;
}
