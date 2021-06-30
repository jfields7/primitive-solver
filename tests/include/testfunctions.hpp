#ifndef TEST_FUNCTIONS_HPP
#define TEST_FUNCTIONS_HPP

#include <eos.hpp>
#include <primitive_solver.hpp>
#include <ps_types.hpp>
#include <testing.hpp>

/// Check that the temperature and energy density equations are consistent.
template<typename EOSPolicy, typename ErrorPolicy>
bool TestTemperatureFromEnergy(Primitive::EOS<EOSPolicy, ErrorPolicy>* eos, 
			Real n, Real T, Real *Y, const Real tol) {
  Real e = eos->GetEnergy(n, T, Y);
  Real T2 = eos->GetTemperatureFromE(n, e, Y);

  Real err = GetError(T, T2);
  if (err > tol) {
    PrintError(T, T2);
    return false;
  }

  return true;
}

/// Check that the temperature and pressure equations are consistent.
template<typename EOSPolicy, typename ErrorPolicy>
bool TestTemperatureFromPressure(Primitive::EOS<EOSPolicy, ErrorPolicy>* eos, 
			Real n, Real T, Real *Y, const Real tol) {
  Real p = eos->GetPressure(n, T, Y);
  Real T2 = eos->GetTemperatureFromP(n, p, Y);

  Real err = GetError(T, T2);
  if (err > tol) {
    PrintError(T, T2);
    return false;
  }

  return true;
}

/// Check that the enthalpy is consistent with the enthalpy calculated
/// directly from pressure and energy density.
template<typename EOSPolicy, typename ErrorPolicy>
bool TestEnthalpy(Primitive::EOS<EOSPolicy, ErrorPolicy>* eos, Real n, Real T, Real *Y, const Real tol) {
  Real h = eos->GetEnthalpy(n, T, Y);
  Real p = eos->GetPressure(n, T, Y);
  Real e = eos->GetEnergy(n, T, Y);

  Real expected = (e + p)/n;
  Real err = GetError(expected, h);
  if (err > tol) {
    PrintError(expected, h);
    return false;
  }

  return true;
}

/// Check that the specific energy is consistent with the energy density.
template<typename EOSPolicy, typename ErrorPolicy>
bool TestSpecificEnergy(Primitive::EOS<EOSPolicy, ErrorPolicy>* eos, Real n, Real T, Real *Y, const Real tol) {
  Real eps = eos->GetSpecificEnergy(n, T, Y);
  Real e = eos->GetEnergy(n, T, Y);
  Real mb = eos->GetBaryonMass();

  Real expected = (e/(n*mb) - 1.0);

  Real err = GetError(expected, eps);
  if (err > tol) {
    PrintError(expected, eps);
    return false;
  }

  return true;
}

/// Check that ConToPrim and PrimToCon are consistent.
template<typename EOSPolicy, typename ErrorPolicy>
bool TestConToPrim(Primitive::PrimitiveSolver<EOSPolicy, ErrorPolicy>* ps, AthenaArray<Real> prim, AthenaArray<Real> cons,
    AthenaArray<Real> bu, AthenaArray<Real> gd, AthenaArray<Real> gu, int i, int j, int k, const Real tol) {
  Real rho_old = prim(IDN, k, j, i);
  Real Wvx_old = prim(IVX, k, j, i);
  Real Wvy_old = prim(IVY, k, j, i);
  Real Wvz_old = prim(IVZ, k, j, i);
  Real p_old = prim(IPR, k, j, i);
  Real T_old = prim(ITM, k, j, i);

  ps->PrimToCon(prim, cons, bu, gd, gu, i, j, k);
  ps->ConToPrim(prim, cons, bu, gd, gu, i, j, k);

  Real rho_new = prim(IDN, k, j, i);
  Real Wvx_new = prim(IVX, k, j, i);
  Real Wvy_new = prim(IVY, k, j, i);
  Real Wvz_new = prim(IVZ, k, j, i);
  Real p_new = prim(IPR, k, j, i);
  Real T_new = prim(ITM, k, j, i);

  Real err_rho = GetError(rho_old, rho_new);
  Real err_Wvx = GetError(Wvx_old, Wvx_new);
  Real err_Wvy = GetError(Wvy_old, Wvy_new);
  Real err_Wvz = GetError(Wvz_old, Wvz_new);
  Real err_p = GetError(p_old, p_new);
  Real err_T = GetError(T_old, T_new);

  bool success = true;
  if (err_rho > tol) {
    std::cout << "  rho\n";
    PrintError(rho_old, rho_new);
    success = false;
  }

  if (err_Wvx > tol) {
    std::cout << "  Wvx\n";
    PrintError(Wvx_old, Wvx_new);
    success = false;
  }
  if (err_Wvy > tol) {
    std::cout << "  Wvy\n";
    PrintError(Wvy_old, Wvy_new);
    success = false;
  }
  if (err_Wvz > tol) {
    std::cout << "  Wvz\n";
    PrintError(Wvz_old, Wvz_new);
    success = false;
  }

  if (err_p > tol) {
    std::cout << "  p\n";
    PrintError(p_old, p_new);
    success = false;
  }

  if (err_T > tol) {
    std::cout << "  T\n";
    PrintError(T_old, T_new);
    success = false;
  }

  // Reset the primitive variables to their old values.
  prim(IDN, k, j, i) = rho_old;
  prim(IVX, k, j, i) = Wvx_old;
  prim(IVY, k, j, i) = Wvy_old;
  prim(IVZ, k, j, i) = Wvz_old;
  prim(IPR, k, j, i) = p_old;
  prim(ITM, k, j, i) = T_old;

  return success;
}

#endif
