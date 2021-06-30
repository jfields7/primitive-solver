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

  // Construct a simple test case of primitive and conserved variables
  // with no additional species, no magnetic field, zero velocity,
  // and no curvature.
  EOS<IdealGas, DoNothing> eos;
  PrimitiveSolver<IdealGas, DoNothing> ps{&eos};
  AthenaArray<Real> prim, cons, bu, gd, gu;
  int nx = 5;
  int ny = nx;
  int nz = nx;
  prim.NewAthenaArray(7, nx, ny, nz);
  cons.NewAthenaArray(8, nx, ny, nz);
  bu.NewAthenaArray(3, nx, ny, nz);
  gd.NewAthenaArray(NMETRIC, nx);
  gu.NewAthenaArray(NMETRIC, nx);
  Real mb = eos.GetBaryonMass();
  // Initialize the new variables.
  for (int k = 0; k < nz; k++) {
    for (int j = 0; j < ny; j++) {
      for (int i = 0; i < nx; i++) {
        prim(IDN, k, j, i) = 10.0;
        prim(ITM, k, j, i) = 5.0;
        prim(IPR, k, j, i) = eos.GetPressure(prim(IDN, k, j, i)/mb, prim(ITM, k, j, i), nullptr);
      }
    }
  }
  for (int i = 0; i < nx; i++) {
    gu(I00, i) = -1.0;
    gd(I00, i) = -1.0;
    gu(I11, i) = 1.0;
    gd(I11, i) = 1.0;
    gu(I22, i) = 1.0;
    gd(I22, i) = 1.0;
    gu(I33, i) = 1.0;
    gd(I33, i) = 1.0;
  }
  tester.RunTest(&TestConToPrim<IdealGas, DoNothing>, 
                 "Conserved to Primitive Test - Flat, Fieldless, Static", 
                 &ps, prim, cons, bu, gd, gu, 0, 0, 0, 1e-10);

  // Test a strong field at zero velocity and curvature.
  for (int k = 0; k < nz; k++) {
    for (int j = 0; j < ny; j++) {
      for (int i = 0; i < nx; i++) {
        bu(IB1, k, j, i) = 10.0;
        bu(IB2, k, j, i) = 10.0;
        bu(IB3, k, j, i) = 10.0;
      }
    }
  }
  tester.RunTest(&TestConToPrim<IdealGas, DoNothing>,
                 "Conserved to Primitive Test - Flat, Static, Strong Field",
                 &ps, prim, cons, bu, gd, gu, 0, 0, 0, 1e-10);

  // Test a strongly relativistic flow with no field or curvature.
  for (int k = 0; k < nz; k++) {
    for (int j = 0; j < ny; j++) {
      for (int i = 0; i < nx; i++) {
        bu(IB1, k, j, i) = 0.0;
        bu(IB2, k, j, i) = 0.0;
        bu(IB3, k, j, i) = 0.0;

        prim(IVX, k, j, i) = 10.0;
        prim(IVY, k, j, i) = 10.0;
        prim(IVZ, k, j, i) = 10.0;
      }
    }
  }
  tester.RunTest(&TestConToPrim<IdealGas, DoNothing>,
                 "Conserved to Primitive Test - Flat, Fieldless, Relativistic Flow",
                 &ps, prim, cons, bu, gd, gu, 0, 0, 0, 1e-10);

  // Test static flow with no field but strong gravity.
  for (int k = 0; k < nz; k++) {
    for (int j = 0; j < ny; j++) {
      for (int i = 0; i < nx; i++) {
        prim(IVX, k, j, i) = 0.0;
        prim(IVY, k, j ,i) = 0.0;
        prim(IVZ, k, j, i) = 0.0;
      }
    }
  }
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

    gu(I00, i) = 1.0/gt;
    gu(I11, i) = 1.0/gx;
    gu(I22, i) = 1.0/gx;
    gu(I33, i) = 1.0/gx;
  }
  tester.RunTest(&TestConToPrim<IdealGas, DoNothing>,
                 "Conserved to Primitive Test - Fieldless, Static, Strong Gravity",
                 &ps, prim, cons, bu, gd, gu, 0, 0, 0, 1e-10);

  tester.PrintSummary();
}
