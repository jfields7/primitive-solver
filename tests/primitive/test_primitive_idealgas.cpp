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
#include <primitive_utility.hpp>

using namespace Primitive;

// Test functions specific to this particular test suite.
// TestConstruction {{{
bool TestConstruction() {
  EOS<IdealGas, DoNothing> eos;
  PrimitiveSolver<IdealGas, DoNothing> ps{&eos};
  return (&eos == ps.GetEOS());
}
// }}}

// InitializeVariables {{{
void InitializeVariables(EOS<IdealGas, DoNothing>& eos, Real prim[NPRIM]) {
  Real mb = eos.GetBaryonMass();
  prim[IDN] = 10.0;
  prim[ITM] = 5.0;
  prim[IPR] = eos.GetPressure(prim[IDN]/mb, prim[ITM], nullptr);
}
// }}}

int main(int argc, char *argv[]) {
  UnitTests tester{"Primitive Solver, Ideal Gas"};

  // Validate that the primitive solver gets constructed as expected.
  tester.RunTest(&TestConstruction, "Construction Test");

  // Construct a simple test case of primitive and conserved variables
  // with no additional species, no magnetic field, zero velocity,
  // and no curvature.
  EOS<IdealGas, DoNothing> eos;
  PrimitiveSolver<IdealGas, DoNothing> ps{&eos};
  Real tol = 1e-10;
  Real prim[NPRIM] = {0.0};
  Real cons[NCONS] = {0.0};
  Real bu[NMAG] = {0.0};
  Real gd[NSPMETRIC] = {0.0};
  Real gu[NSPMETRIC] = {0.0};
  //Real mb = eos.GetBaryonMass();
  // Initialize the new variables on a Minkowski metric.
  InitializeVariables(eos, prim);
  ParticleFractions(prim,MAX_SPECIES);
  MinkowskiMetric(gd, gu);

  tester.RunTest(&TestConToPrim<IdealGas, DoNothing>, 
                 "Conserved to Primitive Test - Flat, Fieldless, Static", 
                 &ps, prim, cons, bu, gd, gu, tol);

  // Test a strong field at zero velocity and curvature.
  StrongField(bu);
  tester.RunTest(&TestConToPrim<IdealGas, DoNothing>,
                 "Conserved to Primitive Test - Flat, Static, Strong Field",
                 &ps, prim, cons, bu, gd, gu, tol);

  // Test a strongly relativistic flow with no field or curvature.
  ZeroField(bu);
  StrongVelocity(prim);
  tester.RunTest(&TestConToPrim<IdealGas, DoNothing>,
                 "Conserved to Primitive Test - Flat, Fieldless, Relativistic Flow",
                 &ps, prim, cons, bu, gd, gu, tol);

  // Test static flow with no field but strong gravity.
  ZeroVelocity(prim);
  SchwarzschildMetric(gd, gu);
  tester.RunTest(&TestConToPrim<IdealGas, DoNothing>,
                 "Conserved to Primitive Test - Fieldless, Static, Strong Gravity",
                 &ps, prim, cons, bu, gd, gu, tol);

  // Test a strongly relativistic flow with a strong field but no curvature.
  StrongField(bu);
  StrongVelocity(prim);
  MinkowskiMetric(gd, gu);
  tester.RunTest(&TestConToPrim<IdealGas, DoNothing>,
                 "Conserved to Primitive Test - Flat, Strong Field and Flow",
                 &ps, prim, cons, bu, gd, gu, tol);

  // Test a strongly relativistic flow with strong gravity but no magnetic field.
  ZeroField(bu);
  SchwarzschildMetric(gd, gu);
  tester.RunTest(&TestConToPrim<IdealGas, DoNothing>,
                 "Conserved to Primitive Test - Fieldless, Strong Flow and Gravity",
                 &ps, prim, cons, bu, gd, gu, tol);

  // Test a static fluid with strong gravity and a strong magnetic field.
  ZeroVelocity(prim);
  StrongField(bu);
  tester.RunTest(&TestConToPrim<IdealGas, DoNothing>,
                 "Conserved to Primitive Test - Static, Strong Field and Gravity",
                 &ps, prim, cons, bu, gd, gu, tol);

  // Test a strongly relativistic flow with strong gravity and a strong magnetic field.
  StrongVelocity(prim);
  tester.RunTest(&TestConToPrim<IdealGas, DoNothing>,
                 "Conserved to Primitive Test - Strong Flow, Field, and Gravity",
                 &ps, prim, cons, bu, gd, gu, tol);

  // Now we need to test how the code behaves in a screwball coordinate system.
  // Static flow, no magnetic field.
  ZeroField(bu);
  ZeroVelocity(prim);
  ScrewballMinkowskiMetric(gd, gu);
  tester.RunTest(&TestConToPrim<IdealGas, DoNothing>,
                 "Conserved to Primitive Test - Screwball Minkowski, Static and Fieldless",
                 &ps, prim, cons, bu, gd, gu, tol);

  // Strong flow, no magnetic field.
  StrongVelocity(prim);
  tester.RunTest(&TestConToPrim<IdealGas, DoNothing>,
                 "Conserved to Primitive Test - Screwball Minkowski, Strong Flow and Fieldless",
                 &ps, prim, cons, bu, gd, gu, tol);

  // Strong field, no flow.
  ZeroVelocity(prim);
  StrongField(bu);
  tester.RunTest(&TestConToPrim<IdealGas, DoNothing>,
                 "Conserved to Primitive Test - Screwball Minkowski, Static and Strong Flow",
                 &ps, prim, cons, bu, gd, gu, tol);

  // Strong field and flow
  StrongVelocity(prim);
  tester.RunTest(&TestConToPrim<IdealGas, DoNothing>,
                 "Conserved to Primitive Test - Screwball Minkowski, Strong Flow and Field",
                 &ps, prim, cons, bu, gd, gu, tol);

  
  // Now we can try the same thing for the Schwarzschild metric.
  ZeroField(bu);
  ZeroVelocity(prim);
  ScrewballSchwarzschildMetric(gd, gu);
  tester.RunTest(&TestConToPrim<IdealGas, DoNothing>,
                 "Conserved to Primitive Test - Screwball Schwarzschild, Static and Fieldless",
                 &ps, prim, cons, bu, gd, gu, tol);

  // Strong flow, no magnetic field.
  StrongVelocity(prim);
  tester.RunTest(&TestConToPrim<IdealGas, DoNothing>,
                 "Conserved to Primitive Test - Screwball Schwarzschild, Strong Flow and Fieldless",
                 &ps, prim, cons, bu, gd, gu, tol);

  // Strong field, no flow.
  ZeroVelocity(prim);
  StrongField(bu);
  tester.RunTest(&TestConToPrim<IdealGas, DoNothing>,
                 "Conserved to Primitive Test - Screwball Schwarzschild, Static and Strong Flow",
                 &ps, prim, cons, bu, gd, gu, tol);

  // Strong field and flow
  StrongVelocity(prim);
  tester.RunTest(&TestConToPrim<IdealGas, DoNothing>,
                 "Conserved to Primitive Test - Screwball Schwarzschild, Strong Flow and Field",
                 &ps, prim, cons, bu, gd, gu, tol);


  // Baryon mass consistency tests.
  // Zero field or flow
  MinkowskiMetric(gd, gu);
  ZeroField(bu);
  ZeroVelocity(prim);
  eos.SetBaryonMass(10.3);
  InitializeVariables(eos, prim);
  tester.RunTest(&TestConToPrim<IdealGas, DoNothing>,
                 "Baryon Mass Consistency Test - Static and Fieldless",
                 &ps, prim, cons, bu, gd, gu, tol);

  // Relativistic flow, no field
  StrongVelocity(prim);
  tester.RunTest(&TestConToPrim<IdealGas, DoNothing>,
                 "Baryon Mass Consistency Test - Strong Flow and Fieldless",
                 &ps, prim, cons, bu, gd, gu, tol);
  
  // Strong field, no flow
  ZeroVelocity(prim);
  StrongField(bu);
  tester.RunTest(&TestConToPrim<IdealGas, DoNothing>,
                 "Baryon Mass Consistency Test - Static and Strong Field",
                 &ps, prim, cons, bu, gd, gu, tol);

  // Relativistic flow and strong field
  StrongVelocity(prim);
  tester.RunTest(&TestConToPrim<IdealGas, DoNothing>,
                 "Baryon Mass Consistency Test - Strong Flow and Field",
                 &ps, prim, cons, bu, gd, gu, tol);

  // One species test
  // Flat, fieldless, static
  eos.SetBaryonMass(1.0);
  InitializeVariables(eos, prim);
  ZeroField(bu);
  ZeroVelocity(prim);
  eos.SetNSpecies(1);
  tester.RunTest(&TestConToPrim<IdealGas, DoNothing>,
                 "One Species Test - Flat, Fieldless, Static",
                 &ps, prim, cons, bu, gd, gu, tol);
  
  // Flat, static, strong field
  StrongField(bu);
  tester.RunTest(&TestConToPrim<IdealGas, DoNothing>,
                 "One Species Test - Flat, Static, Strong Field",
                 &ps, prim, cons, bu, gd, gu, tol);

  // Flat, no field, relativistic flow
  ZeroField(bu);
  StrongVelocity(prim);
  tester.RunTest(&TestConToPrim<IdealGas, DoNothing>,
                 "One Species Test - Flat, Fieldless, Relativistic Flow",
                 &ps, prim, cons, bu, gd, gu, tol);

  // Flat, strong field, strong flow
  StrongField(bu);
  tester.RunTest(&TestConToPrim<IdealGas, DoNothing>,
                 "One Species Test - Flat, Strong Flow and Field",
                 &ps, prim, cons, bu, gd, gu, tol);

  // Static and fieldless, strong gravity
  ZeroField(bu);
  ZeroVelocity(prim);
  SchwarzschildMetric(gd, gu);
  tester.RunTest(&TestConToPrim<IdealGas, DoNothing>,
                 "One Species Test - Static and Fieldless, Strong Gravity",
                 &ps, prim, cons, bu, gd, gu, tol);

  // Static, strong gravity and field
  StrongField(bu);
  tester.RunTest(&TestConToPrim<IdealGas, DoNothing>,
                 "One Species Test - Static, Strong Field and Gravity",
                 &ps, prim, cons, bu, gd, gu, tol);

  // Fieldless, strong flow and gravity
  ZeroField(bu);
  StrongVelocity(prim);
  tester.RunTest(&TestConToPrim<IdealGas, DoNothing>,
                 "One Species Test - Fieldless, Strong Flow and Gravity",
                 &ps, prim, cons, bu, gd, gu, tol);

  // Strong flow, field, and gravity
  StrongField(bu);
  tester.RunTest(&TestConToPrim<IdealGas, DoNothing>,
                 "One Species Test - Strong Flow, Field, and Gravity",
                 &ps, prim, cons, bu, gd, gu, tol);

  tester.PrintSummary();
  return 0;
}
