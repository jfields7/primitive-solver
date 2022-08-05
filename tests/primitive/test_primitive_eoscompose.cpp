//! \file test_primitive_eos_compose.cpp
//  \brief Unit tests for a CompOSE EOS primitive solver.
#include <iostream>
#include <cmath>

#include <eos.hpp>
#include <ps_types.hpp>
#include <eos_compose.hpp>
#include <reset_floor.hpp>
#include <primitive_solver.hpp>

#include <testing.hpp>
#include <testfunctions.hpp>
#include <primitive_utility.hpp>

using namespace Primitive;

// Test functions specific to this particular test suite.
// TestConstruction {{{
bool TestConstruction(EOS<EOSCompOSE, ResetFloor>* eos) {
  try {
    eos->ReadTableFromFile("../data/SFHo.h5");
  } catch (std::exception& e) {
    std::cout << "  Failed to initialize the EOS" << std::endl;
    std::cout << e.what() << std::endl;
    return false;
  }

  PrimitiveSolver<EOSCompOSE, ResetFloor> ps{eos};
  return (eos == ps.GetEOS());
}
// }}}

// InitializeVariables {{{
void InitializeVariables(EOS<EOSCompOSE, ResetFloor>& eos, Real prim[NPRIM]) {
  prim[IDN] = 0.5;
  prim[ITM] = 10.0;
  prim[IYF] = 0.3;
  prim[IPR] = eos.GetPressure(prim[IDN], prim[ITM], &prim[IYF]);
}
// }}}

int main(int argc, char* argv[]) {
  (void)argv[argc-1];

  UnitTests tester("Primitive Solver, CompOSE Tabulated EOS");
  EOS<EOSCompOSE, ResetFloor> eos;

  // Validate that the table is valid and the primitive solver is
  // constructed as expected.
  tester.RunTest(&TestConstruction, "Construction Test", &eos);

  // Construct a simple test case of primitive and conserved variables
  // with Y_e.
  PrimitiveSolver<EOSCompOSE, ResetFloor> ps{&eos};
  Real tol = 1e-10;
  Real prim[NPRIM] = {0.0};
  Real cons[NCONS] = {0.0};
  Real bu[NMAG] = {0.0};
  Real gd[NSPMETRIC] = {0.0};
  Real gu[NSPMETRIC] = {0.0};
  // Initialize the variables on a Minkowski metric
  InitializeVariables(eos, prim);
  MinkowskiMetric(gd, gu);

  tester.RunTest(&TestConToPrim<EOSCompOSE, ResetFloor>,
                 "Conserved to Primitive Test - Flat, Fieldless, Static",
                 &ps, prim, cons, bu, gd, gu, tol);

  // Test a strong field at zero velocity and curvature.
  StrongField(bu);
  tester.RunTest(&TestConToPrim<EOSCompOSE, ResetFloor>,
                 "Conserved to Primitive Test - Flat, Static, Strong Field",
                 &ps, prim, cons, bu, gd, gu, tol);

  // Test a strongly relativistic flow with no field or curvature.
  ZeroField(bu);
  StrongVelocity(prim);
  tester.RunTest(&TestConToPrim<EOSCompOSE, ResetFloor>,
                 "Conserved to Primitive Test - Flat, Fieldless, Relativistic Flow",
                 &ps, prim, cons, bu, gd, gu, tol);

  // Test static flow with no field but strong gravity.
  ZeroVelocity(prim);
  SchwarzschildMetric(gd, gu);
  tester.RunTest(&TestConToPrim<EOSCompOSE, ResetFloor>,
                 "Conserved to Primitive Test - Fieldless, Static, Strong Gravity",
                 &ps, prim, cons, bu, gd, gu, tol);
  
  // Test a strongly relativistic flow with a strong field but no curvature.
  StrongField(bu);
  StrongVelocity(prim);
  MinkowskiMetric(gd, gu);
  tester.RunTest(&TestConToPrim<EOSCompOSE, ResetFloor>,
                 "Conserved to Primitive Test - Flat, Strong Field and Flow",
                 &ps, prim, cons, bu, gd, gu, tol);

  // Test a strongly relativistic flow with strong gravity but no magnetic field.
  ZeroField(bu);
  SchwarzschildMetric(gd, gu);
  tester.RunTest(&TestConToPrim<EOSCompOSE, ResetFloor>,
                 "Conserved to Primitive Test - Fieldless, Strong Flow and Gravity",
                 &ps, prim, cons, bu, gd, gu, tol);

  // Test a static fluid with strong gravity and a strong magnetic field.
  ZeroVelocity(prim);
  StrongField(bu);
  tester.RunTest(&TestConToPrim<EOSCompOSE, ResetFloor>,
                 "Conserved to Primitive Test - Static, Strong Field and Gravity",
                 &ps, prim, cons, bu, gd, gu, tol);

  // Test a strongly relativistic flow with strong gravity and a strong magnetic field.
  StrongVelocity(prim);
  tester.RunTest(&TestConToPrim<EOSCompOSE, ResetFloor>,
                 "Conserved to Primitive Test - Strong Flow, Field, and Gravity",
                 &ps, prim, cons, bu, gd, gu, tol);


  tester.PrintSummary();

  return 0;
}
