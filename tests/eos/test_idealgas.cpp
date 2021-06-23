//! \file test_idealgas.cpp
//  \brief Unit test for IdealGas EOS.

#include <iostream>

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
  }, "Ideal Gas Construction");

  // Number density
  Real n = 100.0;
  // Temperature
  Real T = 100.0;
  // Array of crap for particle fractions
  Real *Y = nullptr;

  tester.PrintSummary();

  return 0;
}
