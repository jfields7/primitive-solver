//! \file test_primitive_unphysical.cpp
//  \brief Unit test for unphysical states in PrimitiveSolver

#include <iostream>
#include <cmath>

#include <eos.hpp>
#include <ps_types.hpp>
#include <idealgas.hpp>
#include <do_nothing.hpp>
#include <primitive_solver.hpp>
#include <geom_math.hpp>
#include <unit_system.hpp>

#include <testing.hpp>
#include <testfunctions.hpp>
#include <primitive_utility.hpp>

using namespace Primitive;

// Take an unphysical input state and compare it to the expected unphysical primitive state.
template<typename EOSPolicy, typename ErrorPolicy>
bool TestUnphysical(PrimitiveSolver<EOSPolicy, ErrorPolicy>* ps, Real cons[NCONS], 
                    Real prim_exp[NPRIM], Real bu[NMAG], Real gd[NMETRIC], Real gu[NMETRIC], 
                    const Real tol) {
  Real prim_actual[NPRIM];
  SolverResult result = ps->ConToPrim(prim_actual, cons, bu, gd, gu);

  if (result.error != Error::SUCCESS) {
    std::cout << "An error occurred during the primitive solve.\n";
    switch(result.error) {
      case Primitive::Error::RHO_TOO_BIG:
        std::cout << "  Rho is too large.\n";
        break;
      case Primitive::Error::RHO_TOO_SMALL:
        std::cout << "  Rho is too small.\n";
        break;
      case Primitive::Error::NANS_IN_CONS:
        std::cout << "  Incoming conservative variables contain NaNs.\n";
        break;
      case Primitive::Error::MAG_TOO_BIG:
        std::cout << "  Incoming magnetic field is too large.\n";
        break;
      case Primitive::Error::BRACKETING_FAILED:
        std::cout << "  Failed to bracket the root correctly.\n";
        break;
      case Primitive::Error::NO_SOLUTION:
        std::cout << "  Failed to find a solution.\n";
        break;
      case Primitive::Error::CONS_FLOOR:
        std::cout << "  Unpermitted conservative flooring occurred.\n";
        break;
      case Primitive::Error::PRIM_FLOOR:
        std::cout << "  Unpermitted primitive flooring occurred.\n";
        break;
      case Primitive::Error::CONS_ADJUSTED:
        std::cout << "  Unpermitted adjustment required in the conserved variables.\n";
        break;
      default:
        std::cout << "  An unknown error occurred.\n";
    }
    return false;
  }

  Real err[NPRIM] = {0.0};
  bool success = true;
  for (int i = 0; i < NPRIM; i++) {
    err[i] = GetError(prim_exp[i], prim_actual[i]);
    // If the error exceeds the tolerance, print it out and fail the test.
    if (err[i] > tol) {
      std::cout << "  prim[" << i << "]\n";
      PrintError(prim_exp[i], prim_actual[i]);
      success = false;
    }
    // If the error is non-finite but the expected value *is*, then there was
    // an issue with the actual value, and we should fail the test.
    if (!std::isfinite(err[i]) && std::fabs(prim_exp[i]) > 0) {
      std::cout << "  Non-finite result in error for prim[" << i << "]!\n";
      success = false;
    }
  }

  if (success == false) {
    Real W = std::sqrt(1.0 + SquareVector(&prim_actual[IVX], gd));
    Real bsq = SquareVector(bu, gd)/(prim_actual[IDN]*W);
    std::cout << "  Other information:\n";
    std::cout << "  Lorentz factor: " << W << "\n";
    std::cout << "  Magnetization: " << bsq << "\n";
  }

  return success;
}

int main(int argc, char *argv[]) {
  (void) argv[argc-1];

  UnitTests tester{"Primitive Solver, Unphysical States"};

  // Initialize all the stuff we'll need, including the EOS, primitive solver,
  // primitive and conserved states, and the metric.
  EOS<IdealGas, DoNothing> eos;
  PrimitiveSolver<IdealGas, DoNothing> ps{&eos};
  Real tol = 1e-10;
  Real prim[NPRIM] = {0.0};
  Real cons[NCONS] = {0.0};
  Real bu[NMAG] = {0.0};
  Real gd[NSPMETRIC] = {0.0};
  Real gu[NSPMETRIC] = {0.0};
  MinkowskiMetric(gd, gu);

  // Set up a state with tau < 0 and S^2 = 0, which should correspond to a negative
  // pressure. Notice that this still obeys the dominant energy condition, as 
  // S^2 = 0 < (D + tau)^2, but we have u < -rho/Gamma, so it results in an unphysical
  // sound speed.
  cons[IDN] = prim[IDN] = 1.0;
  cons[IEN] = -0.5;
  prim[IPR] = cons[IEN]/(1. - eos.GetGamma());
  prim[ITM] = eos.GetTemperatureFromP(prim[IDN], prim[IPR], &prim[IYF]);
  tester.RunTest(&TestUnphysical<IdealGas, DoNothing>,
                 "Unphysical State Test - Negative Tau",
                 &ps, cons, prim, bu, gd, gu, tol);

  // Set up a state with tau < 0 and S^2 = 0 which violates the dominant energy condition.
  cons[IEN] = -2.0;
  prim[IPR] = cons[IEN]/(1. - eos.GetGamma());
  prim[ITM] = eos.GetTemperatureFromP(prim[IDN], prim[IPR], &prim[IYF]);
  tester.RunTest(&TestUnphysical<IdealGas, DoNothing>,
                 "Unphysical State Test - Negative Tau, DEC Violation",
                 &ps, cons, prim, bu, gd, gu, tol);

  // Set up a state with tau = 0 and S^2 > 0 which violates the dominant energy condition.
  cons[IEN] = 0.0;
  cons[IM1] = 2.0;
  cons[IDN] = 1.0;
  prim[IVX] = 0.5;
  prim[IPR] = cons[IM1]/prim[IVX] - prim[IDN];
  prim[ITM] = eos.GetTemperatureFromP(prim[IDN], prim[IPR], &prim[IYF]);
  tester.RunTest(&TestUnphysical<IdealGas, DoNothing>,
                 "Unphysical State Test - Excess Momentum, DEC Violation",
                 &ps, cons, prim, bu, gd, gu, tol);

  tester.PrintSummary();
  return 0;
}
