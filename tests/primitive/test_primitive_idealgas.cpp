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

using namespace Primitive;

// Test functions specific to this particular test suite.
// TestConstruction {{{
bool TestConstruction() {
  EOS<IdealGas, DoNothing> eos;
  PrimitiveSolver<IdealGas, DoNothing> ps{&eos};
  return (&eos == ps.GetEOS());
}
// }}}

// MinkowskiMetric {{{
void MinkowskiMetric(Real gd[NMETRIC], Real gu[NMETRIC]) {
  for (int i = 0; i < NMETRIC; i++) {
    gu[i] = 0.0;
    gd[i] = 0.0;
  }
  
  gd[I00] = -1.0;
  gd[I11] = 1.0;
  gd[I22] = 1.0;
  gd[I33] = 1.0;

  gu[I00] = -1.0;
  gu[I11] = 1.0;
  gu[I22] = 1.0;
  gu[I33] = 1.0;
}
// }}}

// SchwarzschildMetric {{{
void SchwarzschildMetric(Real gd[NMETRIC], Real gu[NMETRIC]) {
  for (int i = 0; i < NMETRIC; i++) {
    gu[i] = 0.0;
    gd[i] = 0.0;
  }
  Real R = 1.0;
  Real rs = 1.0;
  Real hp = 1.0 + rs/(4.0*R);
  Real hm = 1.0 - rs/(4.0*R);
  Real gt = hm*hm/(hp*hp);
  Real gx = hp*hp*hp*hp;

  gd[I00] = -gt;
  gd[I11] = gx;
  gd[I22] = gx;
  gd[I33] = gx;

  gu[I00] = -1.0/gt;
  gu[I11] = 1.0/gx;
  gu[I22] = 1.0/gx;
  gu[I33] = 1.0/gx;
}
// }}}

// ScrewballMinkowskiMetric {{{
void ScrewballMinkowskiMetric(Real gd[NMETRIC], Real gu[NMETRIC]) {
  for (int i = 0; i < NMETRIC; i++) {
    gd[i] = 0.0;
    gu[i] = 0.0;
  }

  // Screwball coordinates:
  // a = t
  // b = x - y
  // c = x + y
  // d = z + x - 3y
  gd[I00] = -1.0;
  gd[I11] = 9.0/2.0;
  gd[I12] = -2.0;
  gd[I13] = -2.0;
  gd[I22] = 3.0/2.0;
  gd[I23] = 1.0;
  gd[I33] = 1.0;

  gu[I00] = 2.0;
  gu[I11] = 2.0;
  gu[I13] = 4.0;
  gu[I22] = 2.0;
  gu[I23] = -2.0;
  gu[I33] = 11.0;
}
// }}}

// ScrewballSchwarzschildMetric {{{
void ScrewballSchwarzschildMetric(Real gd[NMETRIC], Real gu[NMETRIC]) {
  for (int i = 0; i < NMETRIC; i++) {
    gd[i] = 0.0;
    gu[i] = 0.0;
  }
  Real R = 4.0;
  Real rs = 1.0;
  Real cv = 1.0 - rs/R;
  
  // Eddington-Finkelstein coordinates
  gd[I00] = -cv;
  gd[I01] = 2.0;
  gd[I22] = R*R;
  gd[I33] = R*R;

  gu[I01] = 0.5;
  gu[I11] = cv/4.0;
  gu[I22] = 1.0/(R*R);
  gu[I33] = 1.0/(R*R);
}
// }}}

// ZeroVelocity {{{
void ZeroVelocity(Real prim[NPRIM]) {
  prim[IVX] = 0.0;
  prim[IVY] = 0.0;
  prim[IVZ] = 0.0;
}
// }}}

// StrongVelocity {{{
void StrongVelocity(Real prim[NPRIM]) {
  prim[IVX] = 55.0;
  prim[IVY] = 15.0;
  prim[IVZ] = 30.0;
}
// }}}

// ZeroField {{{
void ZeroField(Real bu[NMAG]) {
  bu[IB1] = 0.0;
  bu[IB2] = 0.0;
  bu[IB3] = 0.0;
}
// }}}

// StrongField {{{
void StrongField(Real bu[NMAG]) {
  bu[IB1] = 40.0;
  bu[IB2] = 60.0;
  bu[IB3] = 70.0;
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

// ParticleFractions {{{
void ParticleFractions(Real prim[NPRIM], int s) {
  switch(s) {
    case 3:
      prim[IYF + 2] = 0.25;
    case 2:
      prim[IYF + 1] = 0.25;
    case 1:
      prim[IYF] = 0.25;
      break;
    case 0:
      break;
  }
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
  Real gd[NMETRIC] = {0.0};
  Real gu[NMETRIC] = {0.0};
  Real mb = eos.GetBaryonMass();
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
  /*ZeroField(bu);
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
                 &ps, prim, cons, bu, gd, gu, tol);*/


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
}
