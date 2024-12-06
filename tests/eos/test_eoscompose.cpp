//! \file test_eoscompose.cpp
//  \brief Unit tests for the EOSCompOSE EOS

#include <iostream>
#include <sstream>
#include <cmath>

#include <eos.hpp>
#include <ps_types.hpp>
#include <eos_compose.hpp>
#include <do_nothing.hpp>

#include <testing.hpp>
#include <testfunctions.hpp>

using namespace Primitive;

bool TestConstruction(EOS<EOSCompOSE, DoNothing> * eos, Real tol) {
  try {
    eos->ReadTableFromFile("../data/SFHo.h5");
  } catch (std::exception & e) {
    std::cout << "  Failed to initialize the EOS" << std::endl;
    std::cout << e.what() << std::endl;
    return false;
  }

  Real mb = eos->GetBaryonMass();
  Real mb_want = 939.56535;
  Real err = GetError(mb_want, mb);
  if (err > tol) {
    std::cout << "  Could not read the baryon mass" << std::endl;
    PrintError(mb_want, mb);
    return false;
  }

  Real min_n = eos->GetMinimumDensity();
  Real min_n_want = 1e-5;
  err = GetError(min_n_want, min_n);
  if (err > tol) {
    std::cout << "  Wrong minimum density" << std::endl;
    PrintError(min_n_want, min_n);
    return false;
  }

  Real max_n = eos->GetMaximumDensity();
  Real max_n_want = 1;
  err = GetError(max_n_want, max_n);
  if (err > tol) {
    std::cout << "  Wrong maximum density" << std::endl;
    PrintError(max_n_want, max_n);
    return false;
  }

  Real const * table = eos->GetRawTable();
  Real log_p = table[eos->index(EOSCompOSE::ECLOGP, 3, 2, 5)];
  Real log_p_want = 0.47462794383834883;
  err = GetError(log_p_want, log_p);
  if (err > tol) {
    std::cout << "  Data not read correctly" << std::endl;
    PrintError(log_p_want, log_p);
    return false;
  }

  return true;
}

bool TestPressure(EOS<EOSCompOSE, DoNothing>* peos,
    Real p_want, Real n, Real T, Real* Y, Real tol) {

  Real p = peos->GetPressure(n, T, Y);
  Real err = GetError(p_want, p);
  if (err > tol) {
    std::cout << "  Pressure not interpolated correctly" << std::endl;
    PrintError(p_want, p);
    return false;
  }
  return true;
}

bool TestChemicalPotentials(EOS<EOSCompOSE, DoNothing>* peos,
    Real mub_want, Real muq_want, Real mul_want, Real n, Real T, Real* Y, Real tol) {

  Real mub = peos->GetBaryonChemicalPotential(n, T, Y);
  Real muq = peos->GetChargeChemicalPotential(n, T, Y);
  Real mul = peos->GetElectronLeptonChemicalPotential(n, T, Y);
  
  bool success = true;

  Real err = GetError(mub_want, mub);
  if (err > tol) {
    std::cout << "  Baryon chemical potential not interpolated correctly" << std::endl;
    PrintError(mub_want, mub);
    success = false;
  }

  err = GetError(muq_want, muq);
  if (err > tol) {
    std::cout << "  Charge chemical potential not interpolated correctly" << std::endl;
    PrintError(muq_want, muq);
    success = false;
  }

  err = GetError(mul_want, mul);
  if (err > tol) {
    std::cout << "  Electron-lepton chemical potential not interpolated correctly" << std::endl;
    PrintError(mul_want, mul);
    success = false;
  }

  return success;
}

bool TestNeutrinoEquilibrium(EOS<EOSCompOSE, DoNothing>* peos,
    Real n, Real T, Real* Y, Real tol) {

  Real mul = peos->GetElectronLeptonChemicalPotential(n, T, Y);
  
  bool success = true;

  Real e = 0.0; 
  Real Yl[MAX_SPECIES] = {0.0};

  Real eF  = peos->GetEnergy(n,T,Y);
  Real eta = mul/T;

  Real pi = 3.14159265358979323846;
  Real hc = 1.23984172e3; // [MeV fm]

  Real eta2 = eta*eta;
  Real hc3  = hc*hc*hc;
  Real T3   = T*T*T;
  Real T4   = T3*T;
  Real pi2  = pi*pi;
  Real pi4  = pi2*pi2;

  Real eR = (4*pi/hc3) * T4 * ((7*pi4/20) + 0.5*eta2*(pi2 + 0.5*eta2));
  Real nR = (4*pi/(3*hc3)) * T3 * (eta * (pi2 + eta2));

  e = eF + eR;
  Yl[0] = Y[0] + nR/n;

  Real T_eq;
  Real Y_eq[MAX_SPECIES];

  Real T_guess = std::min(1.5*T,peos->GetMaximumTemperature());
  Real Y_guess[MAX_SPECIES] = {0.0};
  Y_guess[0] = Yl[0];

  bool eq_success = peos->GetBetaEquilibriumTrapped(n, e, Yl, T_eq, Y_eq, T_guess, Y_guess);

  if (!eq_success) {
    std::cout << "  Neutrino equilibrium solve not successful" << std::endl;
  }

  Real err = GetError(T, T_eq);
  if (err > tol) {
    std::cout << "  Temperature found does not match expected" << std::endl;
    PrintError(T, T_eq);
    success = false;
  }

  err = GetError(Y[0], Y_eq[0]);
  if (err > tol) {
    std::cout << "  Electron fraction found does not match expected" << std::endl;
    PrintError(Y[0], Y_eq[0]);
    success = false;
  }

  Real n_nu_eq[3] = {0.0};
  Real e_nu_eq[3] = {0.0};
  peos->GetTrappedNeutrinos(n, T_eq, Y_eq, n_nu_eq, e_nu_eq);

  Real e_eq = peos->GetEnergy(n,T_eq,Y_eq);
  for (int i=0; i<3; ++i){
    e_eq += e_nu_eq[i];
  }
  Real Yle_eq = Y_eq[0] + n_nu_eq[0]/n;

  err = GetError(Yl[0], Yle_eq);
  if (err > tol) {
    std::cout << "  Electron lepton number not conserved in equilibrium" << std::endl;
    PrintError(Yl[0], Yle_eq);
    success = false;
  }

  err = GetError(e, e_eq);
  if (err > tol) {
    std::cout << "  Energy not conserved in equilibrium" << std::endl;
    PrintError(e, e_eq);
    success = false;
  }

  return success;
}

void RunTestSuite(UnitTests& tester, EOS<EOSCompOSE, DoNothing>* peos,
      Real n, Real T, Real* Y, Real tol) {
  std::stringstream ss;

  ss << "Temperature From Energy Test";
  tester.RunTest(&TestTemperatureFromEnergy<EOSCompOSE, DoNothing>,
                 ss.str(),
                 peos, n, T, Y, tol);
  ss.str(std::string());

  ss << "Temperature From Pressure Test";
  tester.RunTest(&TestTemperatureFromPressure<EOSCompOSE, DoNothing>,
                 ss.str(),
                 peos, n, T, Y, tol);
  ss.str(std::string());

  ss << "Enthalpy Test";
  tester.RunTest(&TestEnthalpy<EOSCompOSE, DoNothing>,
                 ss.str(),
                 peos, n, T, Y, tol);
  ss.str(std::string());

  ss << "Specific Energy Test";
  tester.RunTest(&TestSpecificInternalEnergy<EOSCompOSE, DoNothing>,
                 ss.str(),
                 peos, n, T, Y, tol);
  ss.str(std::string());
}

int main(int argc, char ** argv) {
  UnitTests tester("EOSCompOSE Tabulated EOS");
  EOS<EOSCompOSE, DoNothing> eos;

  Real const tol = 1e-12;

  tester.RunTest(&TestConstruction, "Read HDF5 Table", &eos, tol);

  Real np = 0.5;
  Real Tp = 10.0;
  Real Yp[] = {0.3};
  Real p_want = 1.039487191089063e02;
  Real mub_want = 1490.8075594908983;
  Real muq_want = -97.4389607494398;
  Real mul_want = 247.24747752185698;

  tester.RunTest(&TestPressure, "Evaluate Pressure",
    &eos, p_want, np, Tp, Yp, tol);
  tester.RunTest(&TestChemicalPotentials, "Evaluate Chemical Potentials",
    &eos, mub_want, muq_want, mul_want, np, Tp, Yp, tol);
  tester.RunTest(&TestNeutrinoEquilibrium, "Test Neutrino Equilibrium",
    &eos, np, Tp, Yp, tol);

  RunTestSuite(tester, &eos, np, Tp, Yp, tol);

  tester.PrintSummary();

  return 0;
}
