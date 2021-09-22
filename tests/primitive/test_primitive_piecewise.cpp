//! \file test_primitive_piecewise.cpp
//  \brief Unit tests for a PiecewisePolytrope EOS primitive solver.
#include <iostream>
#include <cmath>

#include <eos.hpp>
#include <ps_types.hpp>
#include <piecewise_polytrope.hpp>
#include <do_nothing.hpp>
#include <primitive_solver.hpp>

#include <testing.hpp>
#include <testfunctions.hpp>
#include <primitive_utility.hpp>

using namespace Primitive;

// Test functions specific to this particular test suite.
// TestConstruction {{{
bool TestConstruction() {
  EOS<PiecewisePolytrope, DoNothing> eos;
  PrimitiveSolver<PiecewisePolytrope, DoNothing> ps{&eos};
  return (&eos == ps.GetEOS());
}
// }}}

// InitializeVariables {{{
void InitializeVariables(EOS<PiecewisePolytrope, DoNothing>& eos, Real prim[NPRIM]) {
  // Initialize the EOS.
  const int N = 3;
  Real gamma_pieces[N] = {1.8, 2.3, 1.9};
  Real density_pieces[N] = {1.0, 1.5, 3.0};
  Real mb_nuc = 1.0;
  //Real rho_nuc = 1.0;
  Real rho_min = 0.5;
  Real P0 = 10.0;

  eos.InitializeFromData(density_pieces, gamma_pieces, rho_min, P0, mb_nuc, N);

  Real mb = eos.GetBaryonMass();
  prim[IDN] = rho_min/2.0;
  prim[ITM] = 10.0;
  prim[IPR] = eos.GetPressure(prim[IDN]/mb, prim[ITM], nullptr);
}
// }}}

void SetDensityToPolytrope(Real prim[NPRIM], int i) {
  switch(i) {
    case 0:
      prim[IDN] = 0.75;
      break;
    case 1:
      prim[IDN] = 1.25;
      break;
    case 2:
      prim[IDN] = 2.0;
      break;
    default:
      prim[IDN] = 0.25;
      break;
  }
}

int main(int argc, char *argv[]) {
  UnitTests tester{"Primitive Solver, Piecewise Polytrope"};

  // Validate that the primitive solver gets constructed as expected.
  tester.RunTest(&TestConstruction, "Construction Test");

  // Construct a simple test case of primitive and conserved variables
  // with no additional species, no magnetic field, zero velocity,
  // and no curvature.
  EOS<PiecewisePolytrope, DoNothing> eos;
  PrimitiveSolver<PiecewisePolytrope, DoNothing> ps{&eos};
  Real tol = 1e-10;
  Real prim[NPRIM] = {0.0};
  Real cons[NCONS] = {0.0};
  Real bu[NMAG] = {0.0};
  Real gd[NSPMETRIC] = {0.0};
  Real gu[NSPMETRIC] = {0.0};
  // Initialize the variables on a Minkowski metric.
  InitializeVariables(eos, prim);
  ParticleFractions(prim, MAX_SPECIES);
  MinkowskiMetric(gd, gu);

  tester.RunTest(&TestConToPrim<PiecewisePolytrope, DoNothing>, 
                 "Conserved to Primitive Test - Flat, Fieldless, Static", 
                 &ps, prim, cons, bu, gd, gu, tol);

  // Test a strong field at zero velocity and curvature.
  StrongField(bu);
  tester.RunTest(&TestConToPrim<PiecewisePolytrope, DoNothing>,
                 "Conserved to Primitive Test - Flat, Static, Strong Field",
                 &ps, prim, cons, bu, gd, gu, tol);

  // Test a strongly relativistic flow with no field or curvature.
  ZeroField(bu);
  StrongVelocity(prim);
  tester.RunTest(&TestConToPrim<PiecewisePolytrope, DoNothing>,
                 "Conserved to Primitive Test - Flat, Fieldless, Relativistic Flow",
                 &ps, prim, cons, bu, gd, gu, tol);

  // Test static flow with no field but strong gravity.
  ZeroVelocity(prim);
  SchwarzschildMetric(gd, gu);
  tester.RunTest(&TestConToPrim<PiecewisePolytrope, DoNothing>,
                 "Conserved to Primitive Test - Fieldless, Static, Strong Gravity",
                 &ps, prim, cons, bu, gd, gu, tol);

  // Test a strongly relativistic flow with a strong field but no curvature.
  StrongField(bu);
  StrongVelocity(prim);
  MinkowskiMetric(gd, gu);
  tester.RunTest(&TestConToPrim<PiecewisePolytrope, DoNothing>,
                 "Conserved to Primitive Test - Flat, Strong Field and Flow",
                 &ps, prim, cons, bu, gd, gu, tol);

  // Test a strongly relativistic flow with strong gravity but no magnetic field.
  ZeroField(bu);
  SchwarzschildMetric(gd, gu);
  tester.RunTest(&TestConToPrim<PiecewisePolytrope, DoNothing>,
                 "Conserved to Primitive Test - Fieldless, Strong Flow and Gravity",
                 &ps, prim, cons, bu, gd, gu, tol);

  // Test a static fluid with strong gravity and a strong magnetic field.
  ZeroVelocity(prim);
  StrongField(bu);
  tester.RunTest(&TestConToPrim<PiecewisePolytrope, DoNothing>,
                 "Conserved to Primitive Test - Static, Strong Field and Gravity",
                 &ps, prim, cons, bu, gd, gu, tol);

  // Test a strongly relativistic flow with strong gravity and a strong magnetic field.
  StrongVelocity(prim);
  tester.RunTest(&TestConToPrim<PiecewisePolytrope, DoNothing>,
                 "Conserved to Primitive Test - Strong Flow, Field, and Gravity",
                 &ps, prim, cons, bu, gd, gu, tol);

  tester.PrintSummary();
  return 0;
}
