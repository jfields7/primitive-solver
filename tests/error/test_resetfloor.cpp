//! \test_resetfloor.cpp
//  \brief Unit test for ResetFloor error policy.

#include <iostream>

#include <eos.hpp>
#include <ps_types.hpp>
#include <idealgas.hpp>
#include <reset_floor.hpp>

#include <testing.hpp>
#include <testfunctions.hpp>

// Functions specific to ResetFloor
bool TestConstruction() {
  EOS<IdealGas, ResetFloor> eos;
  Real n_atm = 1e-10;
  Real T_atm = 1.0;
  return (n_atm == eos.GetDensityFloor() && T_atm == eos.GetTemperatureFloor());
}

bool TestPrimitiveFloor(EOS<IdealGas, ResetFloor>* eos, Real n, Real v[3], Real p) {
  Real n_new = n;
  Real v_new[3] = {v[0], v[1], v[2]};
  Real p_new = p;
  Real *Y = nullptr;
  Real T = eos->GetTemperatureFromP(n, p, Y);

  bool result = eos->ApplyPrimitiveFloor(n_new, v_new, p_new, T, Y);

  Real n_atm = eos->GetDensityFloor();
  Real p_atm = eos->GetPressureFloor(Y);

  // Make sure that the test worked when it was supposed to.
  if (!result) {
    if (n < n_atm) {
      std::cout << "  Density was not floored as expected.\n";
      return false;
    }
    if (p < p_atm) {
      std::cout << "  Pressure was not floored as expected.\n";
      return false;
    }
  }
  else {
    if (n >= n_atm && p >= p_atm) {
      std::cout << "  Floor was applied to valid variables.\n";
      std::cout << "  n = " << n << "\n";
      std::cout << "  v = (" << v[0] << "," << v[1] << "," << v[2] <<")\n";
      std::cout << "  p = " << p << "\n";
      std::cout << "  n_new = " << n_new << "\n";
      std::cout << "  v_new = (" << v_new[0] << "," << v_new[1] << "," << v_new[2] <<")\n";
      std::cout << "  p_new = " << p_new << "\n";
      std::cout << "  n_atm = " << n_atm << "\n";
      std::cout << "  p_atm = " << p_atm << "\n";
      return false;
    }
    // If the floor was applied to the density, make sure the variables
    // are all zeroed out the way they should be.
    else if (n < n_atm) {
      if (n_new == n_atm && 
          v_new[0] == 0.0 && v_new[1] == 0.0 && v_new[2] == 0.0 && 
          p_new == p_atm) {
        return true;
      }
      else {
        std::cout << "  Density floor was not applied correctly.\n";
        std::cout << "  n = " << n << "\n";
        std::cout << "  v = (" << v[0] << "," << v[1] << "," << v[2] <<")\n";
        std::cout << "  p = " << p << "\n";
        std::cout << "  n_new = " << n_new << "\n";
        std::cout << "  v_new = (" << v_new[0] << "," << v_new[1] << "," << v_new[2] <<")\n";
        std::cout << "  p_new = " << p_new << "\n";
        std::cout << "  n_atm = " << n_atm << "\n";
        std::cout << "  p_atm = " << p_atm << "\n";
        return false;
      }
    }
    else if (p < p_atm) {
      if (n_new == n &&
          v_new[0] == v[0] && v_new[1] == v[1] && v_new[2] == v[2] &&
          p_new == p_atm) {
        return true;
      }
      else {
        std::cout << "  Pressure floor was not applied correctly.\n";
        std::cout << "  n = " << n << "\n";
        std::cout << "  v = (" << v[0] << "," << v[1] << "," << v[2] <<")\n";
        std::cout << "  p = " << p << "\n";
        std::cout << "  n_new = " << n_new << "\n";
        std::cout << "  v_new = (" << v_new[0] << "," << v_new[1] << "," << v_new[2] <<")\n";
        std::cout << "  p_new = " << p_new << "\n";
        std::cout << "  n_atm = " << n_atm << "\n";
        std::cout << "  p_atm = " << p_atm << "\n";
        return false;
      }
    }
  }
}

int main(int argc, char *argv[]) {
  UnitTests tester("Reset Floor Error Policy");

  // Validate that the constructor works as expected.
  tester.RunTest(&TestConstruction, "Construction Test");

  EOS<IdealGas, ResetFloor> eos;

  // Validate that the primitive variables get floored as expected.
  Real n = 1.0;
  Real v[3] = {0.1, 0.1, 0.1};
  Real p = 1.0;
  tester.RunTest(&TestPrimitiveFloor, "Valid Primitive Floor Test", &eos, n, v, p);
  n = 0.0;
  tester.RunTest(&TestPrimitiveFloor, "Invalid Density Primitive Floor Test", &eos, n, v, p);
  n = 1.0;
  p = 0.0;
  tester.RunTest(&TestPrimitiveFloor, "Invalid Pressure Primitive Floor Test", &eos, n, v, p);

  tester.PrintSummary();
}
