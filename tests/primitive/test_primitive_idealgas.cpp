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
  const int n_species = eos.GetNSpecies();
  return (&eos == ps.GetEOS() && ps.GetNSpecies() == n_species);
}
// }}}

// MinkowskiMetric {{{
void MinkowskiMetric(AthenaArray<Real>& gd, AthenaArray<Real>& gu, int nx) {
  gd.ZeroClear();
  gu.ZeroClear();
  
  for (int i = 0; i < nx; i++) {
    gd(I00, i) = -1.0;
    gd(I11, i) = 1.0;
    gd(I22, i) = 1.0;
    gd(I33, i) = 1.0;

    gu(I00, i) = -1.0;
    gu(I11, i) = 1.0;
    gu(I22, i) = 1.0;
    gu(I33, i) = 1.0;
  }
}
// }}}

// SchwarzschildMetric {{{
void SchwarzschildMetric(AthenaArray<Real> gd, AthenaArray<Real> gu, int nx) {
  gd.ZeroClear();
  gu.ZeroClear();
  Real R = 1.0;
  Real rs = 1.0;
  Real hp = 1.0 + rs/(4.0*R);
  Real hm = 1.0 - rs/(4.0*R);
  Real gt = hm*hm/(hp*hp);
  Real gx = hp*hp*hp*hp;

  for (int i = 0; i < nx; i++) {
    gd(I00, i) = -gt;
    gd(I11, i) = gx;
    gd(I22, i) = gx;
    gd(I33, i) = gx;

    gu(I00, i) = -1.0/gt;
    gu(I11, i) = 1.0/gx;
    gu(I22, i) = 1.0/gx;
    gu(I33, i) = 1.0/gx;
  }
}
// }}}

// ScrewballMinkowskiMetric {{{
void ScrewballMinkowskiMetric(AthenaArray<Real> gd, AthenaArray<Real> gu, int nx) {
  gd.ZeroClear();
  gu.ZeroClear();

  // This horribly useless metric corresponds to the equations
  // a = t - x,
  // b = x + t + y,
  // c = y - z,
  // d = z + y + t.
  // It's not useful for anything, it just has a really ugly
  // metric that tests that everything works correctly.
  for (int i = 0; i < nx; i++) {
    /*gd(I01, i) = -0.5;
    gd(I22, i) = 0.5;
    gd(I33, i) = 0.5;

    gu(I01, i) = -2.0;
    gu(I22, i) = 2.0;
    gu(I33, i) = 2.0;*/

    gd(I00, i) = 2.0/9.0;
    gd(I01, i) = 8.0/9.0;
    gd(I02, i) = -4.0/9.0;
    gd(I03, i) = -1.0/9.0;
    gd(I11, i) = 5.0/9.0;
    gd(I12, i) = -7.0/9.0;
    gd(I13, i) = -4.0/9.0;
    gd(I22, i) = 8.0/9.0;
    gd(I23, i) = 2.0/9.0;
    gd(I33, i) = 5.0/9.0;
    
    gu(I00, i) = 1.0;
    gu(I01, i) = 2.0;
    gu(I02, i) = 2.0;
    gu(I03, i) = 1.0;
    gu(I12, i) = 1.0;
    gu(I22, i) = 3.0;
    gu(I33, i) = 2.0;
  }
}
// }}}

// ScrewballSchwarzschildMetric {{{
void ScrewballSchwarzschildMetric(AthenaArray<Real> gd, AthenaArray<Real> gu, int nx) {
  gd.ZeroClear();
  gu.ZeroClear();
  Real R = 4.0;
  Real rs = 1.0;
  Real cv = 1.0 - rs/R;
  
  // Eddington-Finkelstein coordinates
  for (int i = 0; i < nx; i++) {
    gd(I00, i) = -cv;
    gd(I01, i) = 2.0;
    gd(I22, i) = R*R;
    gd(I33, i) = R*R;

    gu(I01, i) = 0.5;
    gu(I11, i) = cv/4.0;
    gu(I22, i) = 1.0/(R*R);
    gu(I33, i) = 1.0/(R*R);
  }
}
// }}}

// ZeroVelocity {{{
void ZeroVelocity(AthenaArray<Real>& prim, int nx, int ny, int nz) {
  for (int k = 0; k < nz; k++) {
    for (int j = 0; j < ny; j++) {
      for (int i = 0; i < nx; i++) {
        prim(IVX, k, j, i) = 0.0;
        prim(IVY, k, j, i) = 0.0;
        prim(IVZ, k, j, i) = 0.0;
      }
    }
  }
}
// }}}

// StrongVelocity {{{
void StrongVelocity(AthenaArray<Real>& prim, int nx, int ny, int nz) {
  for (int k = 0; k < nz; k++) {
    for (int j = 0; j < ny; j++) {
      for (int i = 0; i < nx; i++) {
        prim(IVX, k, j, i) = 15.0;
        prim(IVY, k, j, i) = 5.0;
        prim(IVZ, k, j, i) = 10.0;
      }
    }
  }
}
// }}}

// ZeroField {{{
void ZeroField(AthenaArray<Real>& bu, int nx, int ny, int nz) {
  for (int k = 0; k < nz; k++) {
    for (int j = 0; j < ny; j++) {
      for (int i = 0; i < nx; i++) {
        bu(IB1, k, j, i) = 0.0;
        bu(IB2, k, j, i) = 0.0;
        bu(IB3, k, j, i) = 0.0;
      }
    }
  }
}
// }}}

// StrongField {{{
void StrongField(AthenaArray<Real>& bu, int nx, int ny, int nz) {
  for (int k = 0; k < nz; k++) {
    for (int j = 0; j < ny; j++) {
      for (int i = 0; i < nx; i++) {
        bu(IB1, k, j, i) = 20.0;
        bu(IB2, k, j, i) = 40.0;
        bu(IB3, k, j, i) = 30.0;
      }
    }
  }
}
// }}}

// InitializeVariables {{{
void InitializeVariables(EOS<IdealGas, DoNothing>& eos, AthenaArray<Real>& prim, int nx, int ny, int nz) {
  Real mb = eos.GetBaryonMass();
  for (int k = 0; k < nz; k++) {
    for (int j = 0; j < ny; j++) {
      for (int i = 0; i < nx; i++) {
        prim(IDN, k, j, i) = 10.0;
        prim(ITM, k, j, i) = 5.0;
        prim(IPR, k, j, i) = eos.GetPressure(prim(IDN, k, j, i)/mb, prim(ITM, k, j, i), nullptr);
      }
    }
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
  AthenaArray<Real> prim, cons, bu, gd, gu;
  int nx = 5;
  int ny = nx;
  int nz = nx;
  Real tol = 1e-12;
  prim.NewAthenaArray(7, nx, ny, nz);
  cons.NewAthenaArray(8, nx, ny, nz);
  bu.NewAthenaArray(3, nx, ny, nz);
  gd.NewAthenaArray(NMETRIC, nx);
  gu.NewAthenaArray(NMETRIC, nx);
  Real mb = eos.GetBaryonMass();
  // Initialize the new variables on a Minkowski metric.
  InitializeVariables(eos, prim, nx, ny, nz);
  MinkowskiMetric(gd, gu, nx);

  tester.RunTest(&TestConToPrim<IdealGas, DoNothing>, 
                 "Conserved to Primitive Test - Flat, Fieldless, Static", 
                 &ps, prim, cons, bu, gd, gu, 0, 0, 0, tol);

  // Test a strong field at zero velocity and curvature.
  StrongField(bu, nx, ny, nz);
  tester.RunTest(&TestConToPrim<IdealGas, DoNothing>,
                 "Conserved to Primitive Test - Flat, Static, Strong Field",
                 &ps, prim, cons, bu, gd, gu, 0, 0, 0, tol);

  // Test a strongly relativistic flow with no field or curvature.
  ZeroField(bu, nx, ny, nz);
  StrongVelocity(prim, nx, ny, nz);
  tester.RunTest(&TestConToPrim<IdealGas, DoNothing>,
                 "Conserved to Primitive Test - Flat, Fieldless, Relativistic Flow",
                 &ps, prim, cons, bu, gd, gu, 0, 0, 0, tol);

  // Test static flow with no field but strong gravity.
  ZeroVelocity(prim, nx, ny, nz);
  SchwarzschildMetric(gd, gu, nx);
  tester.RunTest(&TestConToPrim<IdealGas, DoNothing>,
                 "Conserved to Primitive Test - Fieldless, Static, Strong Gravity",
                 &ps, prim, cons, bu, gd, gu, 0, 0, 0, tol);

  // Test a strongly relativistic flow with a strong field but no curvature.
  StrongField(bu, nx, ny, nz);
  StrongVelocity(prim, nx, ny, nz);
  MinkowskiMetric(gd, gu, nx);
  tester.RunTest(&TestConToPrim<IdealGas, DoNothing>,
                 "Conserved to Primitive Test - Flat, Strong Field and Flow",
                 &ps, prim, cons, bu, gd, gu, 0, 0, 0, tol);

  // Test a strongly relativistic flow with strong gravity but no magnetic field.
  ZeroField(bu, nx, ny, nz);
  SchwarzschildMetric(gd, gu, nx);
  tester.RunTest(&TestConToPrim<IdealGas, DoNothing>,
                 "Conserved to Primitive Test - Fieldless, Strong Flow and Gravity",
                 &ps, prim, cons, bu, gd, gu, 0, 0, 0, tol);

  // Test a static fluid with strong gravity and a strong magnetic field.
  ZeroVelocity(prim, nx, ny, nz);
  StrongField(bu, nx, ny, nz);
  tester.RunTest(&TestConToPrim<IdealGas, DoNothing>,
                 "Conserved to Primitive Test - Static, Strong Field and Gravity",
                 &ps, prim, cons, bu, gd, gu, 0, 0, 0, tol);

  // Test a strongly relativistic flow with strong gravity and a strong magnetic field.
  StrongVelocity(prim, nx, ny, nz);
  tester.RunTest(&TestConToPrim<IdealGas, DoNothing>,
                 "Conserved to Primitive Test - Strong Flow, Field, and Gravity",
                 &ps, prim, cons, bu, gd, gu, 0, 0, 0, tol);

  // Now we need to test how the code behaves in a screwball coordinate system.
  // Static flow, no magnetic field.
  ZeroField(bu, nx, ny, nz);
  ZeroVelocity(prim, nx, ny, nz);
  ScrewballMinkowskiMetric(gd, gu, nx);
  tester.RunTest(&TestConToPrim<IdealGas, DoNothing>,
                 "Conserved to Primitive Test - Screwball Minkowski, Static and Fieldless",
                 &ps, prim, cons, bu, gd, gu, 0, 0, 0, tol);

  // Strong flow, no magnetic field.
  StrongVelocity(prim, nx, ny, nz);
  tester.RunTest(&TestConToPrim<IdealGas, DoNothing>,
                 "Conserved to Primitive Test - Screwball Minkowski, Strong Flow and Fieldless",
                 &ps, prim, cons, bu, gd, gu, 0, 0, 0, tol);

  // Strong field, no flow.
  ZeroVelocity(prim, nx, ny, nz);
  StrongField(bu, nx, ny, nz);
  tester.RunTest(&TestConToPrim<IdealGas, DoNothing>,
                 "Conserved to Primitive Test - Screwball Minkowski, Static and Strong Flow",
                 &ps, prim, cons, bu, gd, gu, 0, 0, 0, tol);

  // Strong field and flow
  StrongVelocity(prim, nx, ny, nz);
  tester.RunTest(&TestConToPrim<IdealGas, DoNothing>,
                 "Conserved to Primitive Test - Screwball Minkowski, Strong Flow and Field",
                 &ps, prim, cons, bu, gd, gu, 0, 0, 0, tol);

  
  // Now we can try the same thing for the Schwarzschild metric.
  ZeroField(bu, nx, ny, nz);
  ZeroVelocity(prim, nx, ny, nz);
  ScrewballSchwarzschildMetric(gd, gu, nx);
  tester.RunTest(&TestConToPrim<IdealGas, DoNothing>,
                 "Conserved to Primitive Test - Screwball Schwarzschild, Static and Fieldless",
                 &ps, prim, cons, bu, gd, gu, 0, 0, 0, tol);

  // Strong flow, no magnetic field.
  StrongVelocity(prim, nx, ny, nz);
  tester.RunTest(&TestConToPrim<IdealGas, DoNothing>,
                 "Conserved to Primitive Test - Screwball Schwarzschild, Strong Flow and Fieldless",
                 &ps, prim, cons, bu, gd, gu, 0, 0, 0, tol);

  // Strong field, no flow.
  ZeroVelocity(prim, nx, ny, nz);
  StrongField(bu, nx, ny, nz);
  tester.RunTest(&TestConToPrim<IdealGas, DoNothing>,
                 "Conserved to Primitive Test - Screwball Schwarzschild, Static and Strong Flow",
                 &ps, prim, cons, bu, gd, gu, 0, 0, 0, tol);

  // Strong field and flow
  StrongVelocity(prim, nx, ny, nz);
  tester.RunTest(&TestConToPrim<IdealGas, DoNothing>,
                 "Conserved to Primitive Test - Screwball Schwarzschild, Strong Flow and Field",
                 &ps, prim, cons, bu, gd, gu, 0, 0, 0, tol);


  // Baryon mass consistency tests.
  // Zero field or flow
  MinkowskiMetric(gd, gu, nx);
  ZeroField(bu, nx, ny, nz);
  ZeroVelocity(prim, nx, ny, nz);
  eos.SetBaryonMass(10.3);
  InitializeVariables(eos, prim, nx, ny, nz);
  tester.RunTest(&TestConToPrim<IdealGas, DoNothing>,
                 "Baryon Mass Consistency Test - Static and Fieldless",
                 &ps, prim, cons, bu, gd, gu, 0, 0, 0, tol);

  // Relativistic flow, no field
  StrongVelocity(prim, nx, ny, nz);
  tester.RunTest(&TestConToPrim<IdealGas, DoNothing>,
                 "Baryon Mass Consistency Test - Strong Flow and Fieldless",
                 &ps, prim, cons, bu, gd, gu, 0, 0, 0, tol);
  
  // Strong field, no flow
  ZeroVelocity(prim, nx, ny, nz);
  StrongField(bu, nx, ny, nz);
  tester.RunTest(&TestConToPrim<IdealGas, DoNothing>,
                 "Baryon Mass Consistency Test - Static and Strong Field",
                 &ps, prim, cons, bu, gd, gu, 0, 0, 0, tol);

  // Relativistic flow and strong field
  StrongVelocity(prim, nx, ny, nz);
  tester.RunTest(&TestConToPrim<IdealGas, DoNothing>,
                 "Baryon Mass Consistency Test - Strong Flow and Field",
                 &ps, prim, cons, bu, gd, gu, 0, 0, 0, tol);

  tester.PrintSummary();
}
