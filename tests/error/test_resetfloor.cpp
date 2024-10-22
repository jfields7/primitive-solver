//! \test_resetfloor.cpp
//  \brief Unit test for ResetFloor error policy.

#include <iostream>

#include <eos.hpp>
#include <ps_types.hpp>
#include <idealgas.hpp>
#include <reset_floor.hpp>

#include <testing.hpp>
#include <testfunctions.hpp>

using namespace Primitive;

// Functions specific to ResetFloor
bool TestConstruction() {
  EOS<IdealGas, ResetFloor> eos;
  Real n_atm = 1e-10;
  Real T_atm = 1e-10;
  Real n_threshold = 1.0;
  return (n_atm == eos.GetDensityFloor() && T_atm == eos.GetTemperatureFloor() &&
          n_threshold == eos.GetThreshold());
}

void PrintDetails(Real n, Real v[3], Real T, Real n_new, Real v_new[3], 
                  Real T_new, Real n_atm, Real t_atm, Real n_threshold) {
  std::cout << "  n = " << n << "\n";
  std::cout << "  v = (" << v[0] << "," << v[1] << "," << v[2] <<")\n";
  std::cout << "  T = " << T << "\n";
  std::cout << "  n_new = " << n_new << "\n";
  std::cout << "  v_new = (" << v_new[0] << "," << v_new[1] << "," << v_new[2] <<")\n";
  std::cout << "  T_new = " << T_new << "\n";
  std::cout << "  n_atm = " << n_atm << "\n";
  std::cout << "  t_atm = " << t_atm << "\n";
  std::cout << "  n_threshold = " << n_threshold << "\n";
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
  Real n_threshold = eos->GetThreshold();

  // Make sure that the test worked when it was supposed to.
  if (!result) {
    if (n < n_atm*n_threshold) {
      std::cout << "  Density was not floored as expected.\n";
      return false;
    }
    if (T < t_atm) {
      std::cout << "  Temperature was not floored as expected.\n";
      return false;
    }
  }
  else {
    if (n >= n_atm*n_threshold && T >= t_atm) {
      std::cout << "  Floor was applied to valid variables.\n";
      PrintDetails(n, v, T, n_new, v_new, T_new, n_atm, t_atm, n_threshold);
      return false;
    }
    // If the floor was applied to the density, make sure the variables
    // are all zeroed out the way they should be.
    else if (n < n_atm*n_threshold) {
      if (n_new == n_atm && 
          v_new[0] == 0.0 && v_new[1] == 0.0 && v_new[2] == 0.0 && 
          T_new == t_atm) {
        return true;
      }
      else {
        std::cout << "  Density floor was not applied correctly.\n";
        PrintDetails(n, v, T, n_new, v_new, T_new, n_atm, t_atm, n_threshold);
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
        PrintDetails(n, v, T, n_new, v_new, T_new, n_atm, t_atm, n_threshold);
        return false;
      }
    }
  }

  return true;
}

bool TestMagnetizationResponse(EOS<IdealGas, ResetFloor>* eos, Real bsq, Real b_u[3]) {
  Real bsq_old = bsq;
  Real b_u_old[3] = {b_u[0], b_u[1], b_u[2]};

  Error result = eos->DoMagnetizationResponse(bsq, b_u);

  Real bsq_max = eos->GetMaximumMagnetization();

  // Make sure that the test worked when it was supposed to.
  if ((result == Error::CONS_ADJUSTED && bsq != bsq_max) || 
      (result == Error::SUCCESS && bsq > bsq_max)) {
    std::cout << "  Magnetization was not rescaled properly.\n";
    std::cout << "  bsq_old = " << bsq_old << "\n";
    std::cout << "  bsq     = " << bsq << "\n";
    std::cout << "  b_u_old = (" << b_u_old[0] << ", " << b_u_old[1] << ", " << b_u_old[2] << ")\n";
    std::cout << "  b_u     = (" << b_u[0] << ", " << b_u[1] << ", " << b_u[2] << ")\n";
    return false;
  }
  return true;
}

bool TestTemperatureLimits(EOS<IdealGas, ResetFloor>* eos, Real n, Real T) {
  //Real e = eos->GetEnergy(n, T, nullptr);

  //Real e_adjusted = e;
  Real T_adjusted = T;
  eos->ApplyTemperatureLimits(T_adjusted);

  Real min_T = eos->GetMinimumTemperature();
  Real max_T = eos->GetMaximumTemperature();
  if (T < min_T && T_adjusted != min_T) {
    std::cout << "  Small temperature was not rescaled properly.\n";
    std::cout << "  Expected: " << min_T << "\n";
    std::cout << "  Actual: " << T_adjusted << "\n";
    return false;
  }
  else if (T > max_T && T_adjusted != max_T) {
    std::cout << "  Large temperature was not rescaled properly.\n";
    std::cout << "  Expected: " << max_T << "\n";
    std::cout << "  Actual: " << T_adjusted << "\n";
    return false;
  }
  else if (T >= min_T && T <= max_T && T != T_adjusted) {
    std::cout << "  Valid temperature was unexpectedly rescaled.\n";
    std::cout << "  Expected: " << T << "\n";
    std::cout << "  Actual: " << T_adjusted << "\n";
    return false;
  }

  return true;
}

bool TestSpeciesLimits(EOS<IdealGas, ResetFloor>* eos, Real *Y) {
  Real Y_adjusted[MAX_SPECIES] = {0.0};
  int n_species = eos->GetNSpecies();
  for (int i = 0; i < n_species; i++) {
    Y_adjusted[i] = Y[i];
  }
  eos->ApplySpeciesLimits(Y_adjusted);
  for (int i = 0; i < n_species; i++) {
    Real min_Y = eos->GetMinimumSpeciesFraction(i);
    Real max_Y = eos->GetMaximumSpeciesFraction(i);
    // Check that the bounds are correct.
    if (Y[i] < min_Y && Y_adjusted[i] != min_Y) {
      std::cout << "  Small species fraction was not rescaled properly.\n";
      std::cout << "  Expected: " << min_Y << "\n";
      std::cout << "  Actual: " << Y_adjusted[i] << "\n";
      return false;
    }
    else if (Y[i] > max_Y && Y_adjusted[i] != max_Y) {
      std::cout << "  Large species fraction was not rescaled properly.\n";
      std::cout << "  Expected: " << max_Y << "\n";
      std::cout << "  Actual: " << Y_adjusted[i] << "\n";
      return false;
    }
    else if (Y[i] >= max_Y && Y[i] <= min_Y && Y[i] != Y_adjusted[i]) {
      std::cout << "  Valid species fraction unexpectedly rescaled.\n";
      std::cout << "  Expected: " << Y[i] << "\n";
      std::cout << "  Actual: " << Y_adjusted[i] << "\n";
      return false;
    }
  }
  return true;
}

bool TestFailureResponse(EOS<IdealGas, ResetFloor>* eos) {
  Real prim[NPRIM];
  bool result = eos->DoFailureResponse(prim);

  Real n_atm = eos->GetDensityFloor();
  Real t_atm = eos->GetTemperatureFloor();
  Real p_atm = eos->GetPressure(n_atm, t_atm, &prim[IYF]);

  if (!result) {
    std::cout << "  Result from error response was unexpectedly false.\n";
    return false;
  }

  if (prim[IDN] != n_atm) {
    std::cout << "  Density was not reset properly.\n";
    std::cout << "  Expected: " << n_atm << "\n";
    std::cout << "  Actual: " << prim[IDN] << "\n";
    return false;
  }
  if (prim[IPR] != p_atm) {
    std::cout << "  Pressure was not reset properly.\n";
    std::cout << "  Expected: " << p_atm << "\n";
    std::cout << "  Actual: " << prim[IPR] << "\n";
    return false;
  }
  if (prim[ITM] != t_atm) {
    std::cout << "  Temperature was not reset properly.\n";
    std::cout << "  Expected: " << t_atm << "\n";
    std::cout << "  Actual: " << prim[ITM] << "\n";
    return false;
  }
  if (prim[IVX] != 0.0 || prim[IVY] != 0.0 || prim[IVZ] != 0.0) {
    std::cout << "  Velocity was not reset properly.\n";
    std::cout << "  Expected: (0, 0, 0)" << "\n";
    std::cout << "  Actual: " << "(" << prim[IVX] << ", "
                                     << prim[IVY] << ", "
                                     << prim[IVZ] << ")\n";
    return false;
  }
  
  return true;
}

bool TestConservedFloor(EOS<IdealGas,ResetFloor>* eos, Real D, Real Sd[3], Real tau) {
  Real tau_floor = eos->GetTauFloor(D, nullptr, 0.0);
  Real D_floor = eos->GetDensityFloor()*eos->GetBaryonMass();
  Real tau_floor_2 = eos->GetTauFloor(D_floor, nullptr, 0.0);

  Real D_new = D;
  Real Sd_new[3] = {Sd[0], Sd[1], Sd[2]};
  Real tau_new = tau;

  eos->ApplyConservedFloor(D_new, Sd_new, tau_new, nullptr, 0.0);

  // Check if D should have been floored.
  if (D < D_floor) {
    bool success = true;
    if (std::fabs(D_new - D_floor) > 1e-15) {
      std::cout << "  Relativistic density was not reset correctly!\n";
      std::cout << "  Expected: " << D_floor << "\n";
      std::cout << "  Actual: " << D_new << "\n";
      success = false;
    }
    if (Sd_new[0] != 0.0 || Sd_new[1] != 0.0 || Sd_new[2] != 0.0) {
      std::cout << "  Momentum density was not reset correctly!\n";
      std::cout << "  Expected: (0, 0, 0)\n";
      std::cout << "  Actual: " << Sd_new[0] << ", "
                                << Sd_new[1] << ", "
                                << Sd_new[2] << ")\n";
      success = false;
    }
    if (tau_new != tau_floor_2) {
      std::cout << "  Tau was not reset correctly!\n";
      std::cout << "  Expected: " << tau_floor_2 << "\n";
      std::cout << "  Actual: " << tau_new << "\n";
      success = false;
    }
    return success;
  }
  else if (tau < tau_floor) {
    bool success = true;
    if (tau_new != tau_floor) {
      std::cout << "  Tau was not reset correctly!\n";
      std::cout << "  Expected: " << tau_floor << "\n";
      std::cout << "  Actual: " << tau_new << "\n";
      success = false;
    }
    if (D_new != D) {
      std::cout << "  Relativistic density was unexpectedly reset during tau floor!\n";
      success = false;
    }
    if (Sd_new[0] != Sd[0] || Sd_new[1] != Sd[1] || Sd_new[2] != Sd[2]) {
      std::cout << "  Momentum density was unexpectedly reset during tau floor!\n";
      success = false;
    }
    return success;
  }
  else {
    if (D != D_new || tau != tau_new || 
        Sd_new[0] != Sd[0] || Sd_new[1] != Sd[1] || Sd_new[2] != Sd[2]) {
      std::cout << "  Conserved variables unexpectedly reset despite being valid!\n";
      return false;
    }
  }
  return true;
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
  tester.RunTest(&TestPrimitiveFloor, "Invalid Temperature Primitive Floor Test", &eos, n, v, T);
  n = 0.0;
  T = 0.0;
  tester.RunTest(&TestPrimitiveFloor, "All Invalid Primitive Floor Test", &eos, n, v, T);
  eos.SetThreshold(10.0);
  n = eos.GetDensityFloor();
  T = 1.0;
  tester.RunTest(&TestPrimitiveFloor, "Sub-Threshold Primitive Floor Test", &eos, n, v, T);

  // Validate that the conserved variables get floored as expected.
  eos.SetThreshold(1.0);
  Real D = 1.0;
  Real Sd[3] = {0.1, 0.1, 0.1};
  Real tau = 1.0;
  tester.RunTest(&TestConservedFloor, "Valid Conserved Floor Test", &eos, D, Sd, tau);
  D = 0.0;
  tester.RunTest(&TestConservedFloor, "Invalid Density Conserved Floor Test", &eos, D, Sd, tau);
  D = 1.0;
  tau = 0.0;
  tester.RunTest(&TestConservedFloor, "Invalid Tau Conserved Floor Test", &eos, D, Sd, tau);
  D = 0.0;
  tau = 0.0;
  tester.RunTest(&TestConservedFloor, "All Invalid Conserved Floor Test", &eos, D, Sd, tau);

  // Do some magnetization tests.
  eos.SetMaximumMagnetization(300.0);
  Real b_u[3] = {11.0, 11.0, 11.0};
  // We'll just hard-code flat Cartesian coordinates here.
  Real bsq = b_u[0]*b_u[0] + b_u[1]*b_u[1] + b_u[2]*b_u[2];
  tester.RunTest(&TestMagnetizationResponse, "Invalid Magnetization Test", &eos, bsq, b_u);
  // Now we'll make sure that a fully saturated magnetization still passes.
  b_u[0] = 10.0;
  b_u[1] = 10.0;
  b_u[2] = 10.0;
  bsq = b_u[0]*b_u[0] + b_u[1]*b_u[1] + b_u[2]*b_u[2];
  tester.RunTest(&TestMagnetizationResponse, "Saturated Magnetization Test", &eos, bsq, b_u);
  // Check that a valid magnetization doesn't get rescaled.
  b_u[0] = 5.0;
  b_u[1] = 5.0;
  b_u[2] = 5.0;
  bsq = b_u[0]*b_u[0] + b_u[1]*b_u[1] + b_u[2]*b_u[2];
  tester.RunTest(&TestMagnetizationResponse, "Valid Magnetization Test", &eos, bsq, b_u);

  // Make sure that energy is rescaled correctly.
  // Valid energy
  n = 10.0;
  T = 3.0;
  tester.RunTest(&TestTemperatureLimits, "Valid Temperature Test", &eos, n, T);
  // Negative temperature 
  T = -1.0;
  tester.RunTest(&TestTemperatureLimits, "Negative Temperature Test", &eos, n, T);

  // Make sure that species are rescaled properly.
  eos.SetNSpecies(3);
  Real Y[MAX_SPECIES] = {0.0};
  Y[0] = -0.5;
  Y[1] = 0.;
  Y[2] = 0.;
  // Too small species test
  tester.RunTest(&TestSpeciesLimits, "Small Species Test", &eos, Y);
  // Too large species test
  Y[0] = 1.5;
  tester.RunTest(&TestSpeciesLimits, "Large Species Test", &eos, Y);

  // Make sure the primitive variables are reset to floor after failure.
  tester.RunTest(&TestFailureResponse, "Failure Response Test", &eos);

  tester.PrintSummary();
}
