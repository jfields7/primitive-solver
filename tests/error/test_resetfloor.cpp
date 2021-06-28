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

bool TestPrimitiveFloor(EOS<IdealGas, ResetFloor>* eos, Real n, Real v[3], Real T) {
  Real n_new = n;
  Real v_new[3] = {v[0], v[1], v[2]};
  Real *Y = nullptr;
  Real T_new = T;
  Real p_new = eos->GetPressure(n, T, Y);

  bool result = eos->ApplyPrimitiveFloor(n_new, v_new, p_new, T_new, Y);

  Real n_atm = eos->GetDensityFloor();
  Real t_atm = eos->GetTemperatureFloor();

  // Make sure that the test worked when it was supposed to.
  if (!result) {
    if (n < n_atm) {
      std::cout << "  Density was not floored as expected.\n";
      return false;
    }
    if (T< t_atm) {
      std::cout << "  Temperature was not floored as expected.\n";
      return false;
    }
  }
  else {
    if (n >= n_atm && T >= t_atm) {
      std::cout << "  Floor was applied to valid variables.\n";
      std::cout << "  n = " << n << "\n";
      std::cout << "  v = (" << v[0] << "," << v[1] << "," << v[2] <<")\n";
      std::cout << "  T = " << T << "\n";
      std::cout << "  n_new = " << n_new << "\n";
      std::cout << "  v_new = (" << v_new[0] << "," << v_new[1] << "," << v_new[2] <<")\n";
      std::cout << "  T_new = " << T_new << "\n";
      std::cout << "  n_atm = " << n_atm << "\n";
      std::cout << "  t_atm = " << t_atm << "\n";
      return false;
    }
    // If the floor was applied to the density, make sure the variables
    // are all zeroed out the way they should be.
    else if (n < n_atm) {
      if (n_new == n_atm && 
          v_new[0] == 0.0 && v_new[1] == 0.0 && v_new[2] == 0.0 && 
          T_new == t_atm) {
        return true;
      }
      else {
        std::cout << "  Density floor was not applied correctly.\n";
        std::cout << "  n = " << n << "\n";
        std::cout << "  v = (" << v[0] << "," << v[1] << "," << v[2] <<")\n";
        std::cout << "  T = " << T << "\n";
        std::cout << "  n_new = " << n_new << "\n";
        std::cout << "  v_new = (" << v_new[0] << "," << v_new[1] << "," << v_new[2] <<")\n";
        std::cout << "  T_new = " << T_new << "\n";
        std::cout << "  n_atm = " << n_atm << "\n";
        std::cout << "  t_atm = " << t_atm << "\n";
        return false;
      }
    }
    else if (T < t_atm) {
      if (n_new == n &&
          v_new[0] == v[0] && v_new[1] == v[1] && v_new[2] == v[2] &&
          T_new == t_atm) {
        return true;
      }
      else {
        std::cout << "  Temperature floor was not applied correctly.\n";
        std::cout << "  n = " << n << "\n";
        std::cout << "  v = (" << v[0] << "," << v[1] << "," << v[2] <<")\n";
        std::cout << "  T = " << T << "\n";
        std::cout << "  n_new = " << n_new << "\n";
        std::cout << "  v_new = (" << v_new[0] << "," << v_new[1] << "," << v_new[2] <<")\n";
        std::cout << "  T_new = " << T_new << "\n";
        std::cout << "  n_atm = " << n_atm << "\n";
        std::cout << "  t_atm = " << t_atm << "\n";
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
  Real T = 1.0;
  tester.RunTest(&TestPrimitiveFloor, "Valid Primitive Floor Test", &eos, n, v, T);
  n = 0.0;
  tester.RunTest(&TestPrimitiveFloor, "Invalid Density Primitive Floor Test", &eos, n, v, T);
  n = 1.0;
  T = 0.0;
  tester.RunTest(&TestPrimitiveFloor, "Invalid Pressure Primitive Floor Test", &eos, n, v, T);

  tester.PrintSummary();
}
