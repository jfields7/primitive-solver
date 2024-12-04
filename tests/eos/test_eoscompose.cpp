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

  Real e = peos->GetEnergy(n,T,Y);
  Real Yl[MAX_SPECIES] = {0.0};
  Yl[0] = Y[0];

  Real T_eq;
  Real Y_eq[MAX_SPECIES];

  bool eq_success = peos->GetBetaEquilibriumTrapped(n, e, Yl, T_eq, Y_eq, T, Y);
  
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

  // std::cout << "  Equilibrium success:" << eq_success << std::endl;
  // std::cout << "  e:" << e << ", Yle" << Yl[0] << std::endl;
  // std::cout << "  T: " << T << " -> " << T_eq << std::endl;
  // std::cout << "  Ye:" << Y[0] << " -> " << Y_eq[0] << std::endl;

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

  RunTestSuite(tester, &eos, np, Tp, Yp, tol);

  tester.PrintSummary();

  return 0;
}
