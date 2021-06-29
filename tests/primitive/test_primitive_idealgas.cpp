//! \file test_primitive_solver.cpp
//  \brief Unit test for PrimitiveSolver object

#include <iostream>
#include <cmath>

#include <eos.hpp>
#include <ps_types.hpp>
#include <idealgas.hpp>
#include <do_nothing.hpp>
#include <primitive_solver.hpp>

#include <testing.hpp>
#include <testfunctions.hpp>

// Test functions specific to this particular test suite.
bool TestConstruction() {
  EOS<IdealGas, DoNothing> eos;
  PrimitiveSolver<IdealGas, DoNothing> ps{&eos};
  const int n_species = eos.GetNSpecies();
  return (&eos == ps.GetEOS() && ps.GetNSpecies() == n_species);
}

int main(int argc, char *argv[]) {
  UnitTests tester{"Primitive Solver"};

  // Validate that the primitive solver gets constructed as expected.
  tester.RunTest(&TestConstruction, "Construction Test");

  tester.PrintSummary();
}
