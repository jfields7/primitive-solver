//! \file test_hybridtable.cpp
//  \brief Unit tests for the HybridTable EOS

#include <iostream>
#include <sstream>
#include <cmath>

#include <eos.hpp>
#include <ps_types.hpp>
#include <hybrid_table.hpp>
#include <do_nothing.hpp>

#include <testing.hpp>
#include <testfunctions.hpp>

using namespace Primitive;

bool TestConstruction(EOS<HybridTable, DoNothing> * eos, Real tol) {
  try {
    eos->ReadTableFromFile("../data/RGSLy4_1D.h5");
    eos->SetThermalGamma(2.0);
  } catch (std::exception & e) {
    std::cout << "  Failed to initialize the EOS" << std::endl;
    std::cout << e.what() << std::endl;
    return false;
  }

  Real mb = eos->GetBaryonMass();
  Real mb_want = 939.5654;
  Real err = GetError(mb_want, mb);
  if (err > tol) {
    std::cout << "  Could not read the baryon mass" << std::endl;
    PrintError(mb_want, mb);
    return false;
  }

  Real thermal_gamma = eos->GetThermalGamma();
  Real thermal_gamma_want = 2.0;
  err = GetError(thermal_gamma_want, thermal_gamma);
  if (err > tol) {
    std::cout << "  Could not read the thermal adiabatic index" << std::endl;
    PrintError(thermal_gamma_want, thermal_gamma);
    return false;
  }

  Real min_n = eos->GetMinimumDensity();
  Real min_n_want = 1e-7;
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
  Real log_p = table[eos->index(HybridTable::ECLOGP, 13)];
  Real log_p_want = 3.041156875973929;
  err = GetError(log_p_want, log_p);
  if (err > tol) {
    std::cout << "  Data not read correctly" << std::endl;
    PrintError(log_p_want, log_p);
    return false;
  }

  return true;
}

bool TestPressure(EOS<HybridTable, DoNothing>* peos,
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

void RunTestSuite(UnitTests& tester, EOS<HybridTable, DoNothing>* peos,
      Real n, Real T, Real* Y, Real tol) {
  std::stringstream ss;

  ss << "Temperature From Energy Test";
  tester.RunTest(&TestTemperatureFromEnergy<HybridTable, DoNothing>,
                 ss.str(),
                 peos, n, T, Y, tol);
  ss.str(std::string());

  ss << "Temperature From Pressure Test";
  tester.RunTest(&TestTemperatureFromPressure<HybridTable, DoNothing>,
                 ss.str(),
                 peos, n, T, Y, tol);
  ss.str(std::string());

  ss << "Enthalpy Test";
  tester.RunTest(&TestEnthalpy<HybridTable, DoNothing>,
                 ss.str(),
                 peos, n, T, Y, tol);
  ss.str(std::string());

  ss << "Specific Energy Test";
  tester.RunTest(&TestSpecificInternalEnergy<HybridTable, DoNothing>,
                 ss.str(),
                 peos, n, T, Y, tol);
  ss.str(std::string());
}

int main(int argc, char ** argv) {
  UnitTests tester("HybridTable Tabulated EOS");
  EOS<HybridTable, DoNothing> eos;

  Real const tol = 1e-12;

  tester.RunTest(&TestConstruction, "Read HDF5 Table", &eos, tol);

  Real np = 0.5;
  Real Tp = 10.0;
  Real Yp[] = {0.0};
  Real p_want = 7.991029653148705e01;
  tester.RunTest(&TestPressure, "Evaluate Pressure",
    &eos, p_want, np, Tp, Yp, tol);

  RunTestSuite(tester, &eos, np, Tp, Yp, tol);

  tester.PrintSummary();

  return 0;
}
