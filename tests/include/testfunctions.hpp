#ifndef TEST_FUNCTIONS_HPP
#define TEST_FUNCTIONS_HPP
#include <cmath>

#include <eos.hpp>
#include <primitive_solver.hpp>
#include <ps_types.hpp>
#include <ps_error.hpp>
#include <testing.hpp>
#include <geom_math.hpp>

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

  Real expected = (e + p)/(n*eos->GetBaryonMass());
  Real err = GetError(expected, h);
  if (err > tol) {
    PrintError(expected, h);
    return false;
  }

  return true;
}

/// Check that the specific energy is consistent with the energy density.
template<typename EOSPolicy, typename ErrorPolicy>
bool TestSpecificInternalEnergy(Primitive::EOS<EOSPolicy, ErrorPolicy>* eos,
                                Real n, Real T, Real *Y, const Real tol) {
  Real eps = eos->GetSpecificInternalEnergy(n, T, Y);
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
bool TestConToPrim(Primitive::PrimitiveSolver<EOSPolicy, ErrorPolicy>* ps, Real prim[NPRIM], Real cons[NCONS],
    Real bu[NMAG], Real gd[NMETRIC], Real gu[NMETRIC], const Real tol) {
  Real rho_old = prim[IDN];
  Real Wvx_old = prim[IVX];
  Real Wvy_old = prim[IVY];
  Real Wvz_old = prim[IVZ];
  Real p_old = prim[IPR];
  Real T_old = prim[ITM];
  const int n_species = ps->GetEOS()->GetNSpecies();
  Real Y_old[MAX_SPECIES] = {0.0};
  for (int s = 0; s < n_species; s++) {
    Y_old[s] = prim[IYF + s];
  }

  ps->PrimToCon(prim, cons, bu, gd);
  Primitive::SolverResult result = ps->ConToPrim(prim, cons, bu, gd, gu);

  if(result.error != Primitive::Error::SUCCESS) {
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

  Real rho_new = prim[IDN];
  Real Wvx_new = prim[IVX];
  Real Wvy_new = prim[IVY];
  Real Wvz_new = prim[IVZ];
  Real p_new = prim[IPR];
  Real T_new = prim[ITM];
  Real Y_new[MAX_SPECIES] = {0.0};
  for (int s = 0; s < n_species; s++) {
    Y_new[s] = prim[IYF + s];
  }

  Real err_rho = GetError(rho_old, rho_new);
  Real err_Wvx = GetError(Wvx_old, Wvx_new);
  Real err_Wvy = GetError(Wvy_old, Wvy_new);
  Real err_Wvz = GetError(Wvz_old, Wvz_new);
  Real err_p = GetError(p_old, p_new);
  Real err_T = GetError(T_old, T_new);
  Real err_Y[MAX_SPECIES] = {0.0};
  for (int s = 0; s < n_species; s++) {
    err_Y[s] = GetError(Y_old[s], Y_new[s]);
  }

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
  for (int s = 0; s < n_species; s++) {
    if (err_Y[s] > tol) {
      std::cout << "  Y[\n" << s << "]\n";
      PrintError(Y_old[s], Y_new[s]);
      success = false;
    }
  }

  if (success == false) {
    Real uu[3] = {Wvx_old, Wvy_old, Wvz_old};
    Real W = std::sqrt(1.0 + Primitive::SquareVector(uu, gd));
    Real bsq = Primitive::SquareVector(bu, gd)/(rho_old*W);
    std::cout << "  Other information:\n";
    std::cout << "  Lorentz factor: " << W << "\n";
    std::cout << "  Magnetization: " << bsq << "\n";
  }

  // Reset the primitive variables to their old values.
  prim[IDN] = rho_old;
  prim[IVX] = Wvx_old;
  prim[IVY] = Wvy_old;
  prim[IVZ] = Wvz_old;
  prim[IPR] = p_old;
  prim[ITM] = T_old;
  for (int s = 0; s < n_species; s++) {
    prim[IYF + s] = Y_old[s];
  }

  return success;
}

#endif
