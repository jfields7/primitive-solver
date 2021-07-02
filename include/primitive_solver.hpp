#ifndef PRIMITIVE_SOLVER_HPP
#define PRIMITIVE_SOLVER_HPP

//! \file primitive_solver.hpp
//  \brief Declares PrimitiveSolver class.
//
//  PrimitiveSolver contains all the infrastructure for the inversion
//  procedure from conserved to primitive variables in GRMHD. This
//  particular implementation is based on the solver described in
//  Kastaun et al., Phys. Rev. D 103, 023018 (2021).

#include <cmath>
// FIXME: Debug only!
#include <iostream>

#include <numtoolsroot.h>

#include <eos.hpp>
#include <athena_arrays.hpp>
#include <geom_math.hpp>

namespace Primitive {

template<typename EOSPolicy, typename ErrorPolicy>
class PrimitiveSolver {
  private:
    /// A constant pointer to the EOS.
    /// We make this constant because the
    /// possibility of changing the EOS
    /// during implementation seems both
    /// unlikely and dangerous.
    EOS<EOSPolicy, ErrorPolicy> *const peos;
    
    /// The number of separate particle species in the EOS.
    const int n_species;

    /// The minimum enthalpy for the EOS.
    Real min_h;

    //! \brief function for the upper bound of the root
    //
    //  The upper bound is the solution to the function
    //  \f$\mu\sqrt{h_0^2 + \bar{r}^2(\mu)} - 1 = 0\f$
    //
    //  \param[out] f    The value of the root function at mu
    //  \param[out] df   The derivative of the root function
    //  \param[in   mu   The guess for the root
    //  \param[in]  bsq  The square magnitude of the magnetic field
    //  \param[in]  rsq  The square magnitude of the specific momentum S/D
    //  \param[in]  rbsq The square of the product \f$r\cdot b\f$
    //  \param[in]  h_min The minimum enthalpy
    static void UpperRoot(Real &f, Real &df, Real mu, Real bsq, Real rsq, Real rbsq, Real h_min);

    //! \brief master function for the root solve
    //
    //  The root solve is based on the master function
    //  \f$ f(\mu) = \mu - \hat{mu}(\mu).
    //
    //  \param[in]  mu   The guess for the root
    //  \param[in]  D    The relativistic density
    //  \param[in]  q    The specific tau variable tau/D
    //  \param[in]  bsq  The square magnitude of the magnetic field
    //  \param[in]  rsq  The square magnitude of the specific momentum S/D
    //  \param[in]  rbsq The square of the product \f$r\cdot b\f$
    //  \param[in]  Y    The vector of particle fractions
    //  \param[in]  peos The pointer to the EOS
    //  \param[out] n    The resulting estimate for n
    //  \param[out] T    The resulting estimate for T
    //  \param[out] P    The resulting estimate for P
    //  
    //  \return f evaluated at mu for the given constants.
    static Real RootFunction(Real mu, Real D, Real q, Real bsq, Real rsq, Real rbsq, Real *Y,
        EOS<EOSPolicy, ErrorPolicy> *const peos, Real* n, Real* T, Real* P);
  public:
    /// Constructor
    PrimitiveSolver(EOS<EOSPolicy, ErrorPolicy> *eos) : peos(eos), n_species(eos->GetNSpecies()) {
      min_h = (peos->GetMinimumEnthalpy())/(peos->GetBaryonMass());
    }

    /// Destructor
    ~PrimitiveSolver() = default;

    //! \brief Get the primitive variables from the conserved variables.
    //
    //  \param[out]    prim  The array of primitive variables
    //  \param[in,out] cons  The array of conserved variables
    //  \param[in,out] bu    The magnetic field
    //  \param[in]     gd    The full 4x4 metric
    //  \param[in]     gu    The full 4x4 inverse metric
    //  \param[in]     i,j,k The position in the array
    //
    //  \return success or failure
    bool ConToPrim(AthenaArray<Real>& prim, AthenaArray<Real>& cons,
                   AthenaArray<Real>& b, AthenaArray<Real>& gd,
                   AthenaArray<Real>& gu, int i, int j, int k);

    //! \brief Get the conserved variables from the primitive variables.
    //
    //  \param[in]    prim  The array of primitive variables
    //  \param[out]   cons  The array of conserved variables
    //  \param[in]    bu    The magnetic field
    //  \param[in]    gd    The full 4x4 metric
    //  \param[in]    gu    The full 4x4 inverse metric
    //  \param[in]    i,j,k The position in the array
    //
    //  \return success or failure
    bool PrimToCon(AthenaArray<Real>& prim, AthenaArray<Real>& cons,
                   AthenaArray<Real>& bu, AthenaArray<Real>& gd,
                   AthenaArray<Real>& gu, int i, int j, int k);

    /// Get the EOS used by this PrimitiveSolver.
    inline EOS<EOSPolicy, ErrorPolicy> *const GetEOS() const {
      return peos;
    }

    /// Get the number of species this PrimitiveSolver expects.
    inline const int GetNSpecies() const {
      return n_species;
    }
};

// UpperRoot {{{
template<typename EOSPolicy, typename ErrorPolicy>
void PrimitiveSolver<EOSPolicy, ErrorPolicy>::UpperRoot(Real &f, Real &df, Real mu, Real bsq, Real rsq, Real rbsq, Real min_h) {
  const Real x = 1.0/(1.0 + mu*bsq);
  const Real xsq = x*x;
  const Real rbarsq = rsq*xsq + mu*x*(1.0 + x)*rbsq;
  const Real dis = std::sqrt(min_h*min_h + rbarsq);
  const Real dx = -bsq*xsq;
  //const Real drbarsq = rbsq*x*(1.0 + x) + (mu*rbsq + 2.0*(mu*rbsq + rsq)*x)*dx;
  const Real drbarsq = rbsq*xsq + mu*rbsq*dx + x*(rbsq + 2.0*(mu*rbsq + rsq)*dx);
  f = mu*dis - 1.0;
  df = dis + mu*drbarsq/(2.0*dis);
}
// }}}

// RootFunction {{{
template<typename EOSPolicy, typename ErrorPolicy>
Real PrimitiveSolver<EOSPolicy, ErrorPolicy>::RootFunction(Real mu, Real D, Real q, Real bsq, Real rsq, Real rbsq, Real *Y,
      EOS<EOSPolicy, ErrorPolicy> *const peos, Real* n, Real* T, Real* P) {
  // We need to get some utility quantities first.
  const Real x = 1.0/(1.0 + mu*bsq);
  const Real xsq = x*x;
  const Real musq = mu*mu;
  const Real den = 1.0 + mu*bsq;
  const Real mux = mu/den;
  const Real muxsq = mux/den;
  const Real rbarsq = rsq*xsq + mu*x*(1.0 + x)*rbsq;
  // An alternative calculation of rbarsq that may be more accurate.
  //const Real rbarsq = rsq*xsq + (mux + muxsq)*rbsq;
  //const Real rbarsq = x*(rsq*x + mu*(x + 1.0)*rbsq);
  const Real qbar = q - 0.5*bsq - 0.5*musq*xsq*(bsq*rsq - rbsq);
  const Real mb = peos->GetBaryonMass();

  // Now we can estimate the velocity.
  const Real v_max = peos->GetMaxVelocity();
  const Real vhatsq = std::fmin(musq*rbarsq, v_max*v_max);

  // Using the velocity estimate, predict the Lorentz factor.
  //const Real What = 1.0/std::sqrt(1.0 - vhatsq);
  const Real iWhat = std::sqrt(1.0 - vhatsq);

  // Now estimate the number density.
  Real rhohat = D*iWhat;
  Real nhat = rhohat/mb;
  // TODO: Limit nhat to a physical regime.

  // Estimate the energy density.
  Real eoverD = qbar - mu*rbarsq + 1.0;
  Real ehat = D*eoverD;
  // TODO: Limit ehat to a physical regime.

  // Now we can get an estimate of the temperature, and from that, the pressure and enthalpy.
  Real That = peos->GetTemperatureFromE(nhat, ehat, Y);
  Real Phat = peos->GetPressure(nhat, That, Y);
  Real hhat = peos->GetEnthalpy(nhat, That, Y)/mb;

  // Now we can get two different estimates for nu = h/W.
  Real nu_a = hhat*iWhat;
  Real nu_b = eoverD + Phat/D;
  Real nuhat = std::fmax(nu_a, nu_b);

  // Finally, we can get an estimate for muhat.
  Real muhat = 1.0/(nuhat + mu*rbarsq);

  *n = nhat;
  *T = That;
  *P = Phat;

  // FIXME: Debug only!
  /*std::cout << "    D   = " << D << "\n";
  std::cout << "    q   = " << q << "\n";
  std::cout << "    bsq = " << bsq << "\n";
  std::cout << "    rsq = " << rsq << "\n";
  std::cout << "    rbsq = " << rbsq << "\n";*/

  return mu - muhat;
}
// }}}

// ConToPrim {{{
template<typename EOSPolicy, typename ErrorPolicy>
bool PrimitiveSolver<EOSPolicy, ErrorPolicy>::ConToPrim(AthenaArray<Real>& prim,
       AthenaArray<Real>& cons, AthenaArray<Real>& b, AthenaArray<Real>& gd,
       AthenaArray<Real>& gu, int i, int j, int k) {

  // Extract the 3-metric and inverse 3-metric.
  const Real g3d[NSPMETRIC] = {gd(I11, i), gd(I12, i), gd(I13, i),
                               gd(I22, i), gd(I23, i), gd(I33, i)};
  const Real ialphasq = -gu(I00, i);
  const Real alphasq = 1.0/ialphasq; // Lapse squared
  const Real beta_u[3] = {gu(I01, i)*alphasq, gu(I02, i)*alphasq, gu(I03, i)*alphasq}; // Shift vector
  const Real g3u[NSPMETRIC] = {gu(I11, i) + beta_u[0]*beta_u[0]*ialphasq,
                               gu(I12, i) + beta_u[0]*beta_u[1]*ialphasq,
                               gu(I13, i) + beta_u[0]*beta_u[2]*ialphasq,
                               gu(I22, i) + beta_u[1]*beta_u[1]*ialphasq,
                               gu(I23, i) + beta_u[1]*beta_u[2]*ialphasq,
                               gu(I33, i) + beta_u[2]*beta_u[2]*ialphasq};

  // Get the inverse volume element of the 3-metric.
  Real isdetg = 1.0/std::sqrt(GetDeterminant(g3d));

  // Extract the undensitized conserved variables.
  Real D = cons(IDN, k, j, i)*isdetg;
  Real S_d[3] = {cons(IM1, k, j, i)*isdetg, cons(IM2, k, j, i)*isdetg, cons(IM3, k, j, i)*isdetg};
  Real tau = cons(IEN, k, j, i)*isdetg;
  // FIXME: Confirm that the magnetic field is densitized as well.
  //Real B_u[3] = {b(IB1, k, j, i)*isdetg, b(IB2, k, j, i)*isdetg, b(IB3, k, j, i)*isdetg};
  Real B_u[3] = {b(IB1, k, j, i), b(IB2, k, j, i), b(IB3, k, j, i)};
  // Extract the particle fractions.
  Real Y[n_species] = {0.0};
  for (int s = 0; s < n_species; s++) {
    Y[s] = cons(IYD + s, k, j, i)/cons(IDN, k, j, i);
  }

  // TODO
  // If D is below the atmosphere, we need to do whatever
  // the EOSPolicy wants us to do.

  // Calculate some utility quantities.
  Real sqrtD = std::sqrt(D);
  Real b_u[3] = {B_u[0]/sqrtD, B_u[1]/sqrtD, B_u[2]/sqrtD};
  Real r_d[3] = {S_d[0]/D, S_d[1]/D, S_d[2]/D};
  Real r_u[3];
  RaiseForm(r_u, r_d, g3u);
  Real rsqr   = Contract(r_u, r_d);
  Real rb     = Contract(b_u, r_d);
  Real rbsqr  = rb*rb;
  Real bsqr   = SquareVector(b_u, g3d);
  Real q      = tau/D;

  // Make sure there are no NaNs at this point.
  if (!std::isfinite(D) || !std::isfinite(rsqr) || !std::isfinite(q) ||
      !std::isfinite(rbsqr) || !std::isfinite(bsqr)) {
    return false;    
  }
  // We have to check the particle fractions separately.
  for (int s = 0; s < n_species; s++) {
    if (!std::isfinite(Y[s])) {
      return false;
    }
  }

  // TODO
  // Make sure that the magnetic field is physical.
  
  // Bracket the root.
  Real mul = 0.0;
  Real muh = 1.0/min_h;
  // Check if a tighter upper bound exists.
  if(rsqr > min_h*min_h) {
    Real mu = 0.0;
    // We don't need the bound to be that tight, so we reduce
    // the accuracy of the root solve for speed reasons.
    NumTools::Root::tol = 1e-3;
    NumTools::Root::iterations = 10;
    bool result = NumTools::Root::newton_raphson(&UpperRoot, mu, bsqr, rsqr, rbsqr, min_h);
    // Scream if the bracketing failed.
    if (!result) {
      return false;
    }
    else {
      muh = mu;
    }
  }

  // TODO: There's a corner case here that needs to be checked.
  
  // Do the root solve.
  // TODO: This should be done with something like TOMS748 once it's
  // available.
  NumTools::Root::tol = 1e-15;
  NumTools::Root::iterations = 30;
  Real n, P, T, mu;
  bool result = NumTools::Root::false_position(&RootFunction, mul, muh, mu, D, q, bsqr, rsqr, rbsqr, Y, peos, &n, &T, &P);
  if (!result) {
    return false;
  }

  // Retrieve the primitive variables.
  Real rho = n*peos->GetBaryonMass();
  Real rbmu = rb*mu;
  Real W = D/rho;
  Real Wmux = W*mu/(1.0 + mu*bsqr);
  prim(IDN, k, j, i) = rho;
  prim(IPR, k, j, i) = P;
  prim(ITM, k, j, i) = T;
  // Before we retrieve the velocity, we need to raise S.
  Real S_u[3] = {0.0};
  RaiseForm(S_u, S_d, g3u);
  // Now we can get Wv.
  prim(IVX, k, j, i) = Wmux*(r_u[0] + rbmu*b_u[0]);
  prim(IVY, k, j, i) = Wmux*(r_u[1] + rbmu*b_u[1]);
  prim(IVZ, k, j, i) = Wmux*(r_u[2] + rbmu*b_u[2]);

  // TODO: We probably need to check here for some physical violations.

  return true;
}
// }}}

// PrimToCon {{{
template<typename EOSPolicy, typename ErrorPolicy>
bool PrimitiveSolver<EOSPolicy, ErrorPolicy>::PrimToCon(AthenaArray<Real>& prim,
       AthenaArray<Real>& cons, AthenaArray<Real>& bu, AthenaArray<Real>& gd,
       AthenaArray<Real>& gu,int i, int j, int k) {
  // Extract the three metric.
  const Real g3d[NSPMETRIC] = {gd(I11, i), gd(I12, i), gd(I13, i),
                               gd(I22, i), gd(I23, i), gd(I33, i)};
  
  // Get the volume element of the 3-metric.
  Real sdetg = std::sqrt(GetDeterminant(g3d));

  // Extract the primitive variables
  const Real &rho = prim(IDN, k, j, i); // rest-mass density
  const Real Wv_u[3] = {prim(IVX, k, j, i), prim(IVY, k, j, i), prim(IVZ, k, j, i)};
  const Real &p   = prim(IPR, k, j, i); // pressure
  const Real &t   = prim(ITM, k, j, i); // temperature
  const Real B_u[3] = {bu(IB1, k, j, i), bu(IB2, k, j, i), bu(IB3, k, j, i)};
  Real Y[n_species] = {0.0};
  for (int s = 0; s < n_species; s++) {
    Y[s] = prim(IYF + s, k, j, i);
  }

  // Note that Athena passes in Wv, not v.
  // Lower u.
  Real Wv_d[3];
  LowerVector(Wv_d, Wv_u, g3d);
  Real Wvsq = Contract(Wv_u, Wv_d);
  Real Wsq = 1.0 + Wvsq;
  Real W = sqrt(Wsq);
  // Get the 3-velocity.
  Real v_d[3] = {Wv_d[0]/W, Wv_d[1]/W, Wv_d[2]/W};

  // For the magnetic field contribution, we need to find
  // B_i, B^2, and B^i*v_i.
  Real B_d[3];
  LowerVector(B_d, B_u, g3d);
  Real Bsq = Contract(B_u, B_d);
  Real Bv = Contract(B_u, v_d);

  // Some utility quantities that will be helpful.
  const Real mb = peos->GetBaryonMass();

  // Extract the conserved variables
  Real &D = cons(IDN, k, j ,i); // Relativistic density
  Real &Sx = cons(IM1, k, j, i); // Relativistic momentum density (x)
  Real &Sy = cons(IM2, k, j, i); // Relativistic momentum density (y)
  Real &Sz = cons(IM3, k, j, i); // Relativistic momentum density (z)
  Real &tau = cons(IEN, k, j, i); // Relativistic energy - D

  // Set the conserved quantities.
  // Total enthalpy density
  Real H = rho*peos->GetEnthalpy(rho/mb, t, Y)/mb;
  Real HWsq = H*Wsq;
  D = sdetg*rho*W;
  for (int s = 0; s < n_species; s++) {
    cons(IYD + s, k, j, i) = D*Y[s];
  }
  Real HWsqpb = HWsq + Bsq;
  Sx = sdetg*(HWsqpb*v_d[0] - Bv*B_d[0]);
  Sy = sdetg*(HWsqpb*v_d[1] - Bv*B_d[1]);
  Sz = sdetg*(HWsqpb*v_d[2] - Bv*B_d[2]);
  tau = sdetg*(HWsqpb - p - 0.5*(Bv*Bv + Bsq/Wsq)) - D;
}
// }}}

}; // namespace
#endif
