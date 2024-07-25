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
#include <algorithm>

#include <numtools_root.hpp>

#include <eos.hpp>
#include <geom_math.hpp>
#include <ps_error.hpp>

namespace Primitive {

template<typename EOSPolicy, typename ErrorPolicy>
class PrimitiveSolver {
  private:
    // Inner classes defining functors
    // UpperRootFunctor {{{
    class UpperRootFunctor {
      public:
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
      inline void operator()(Real &f, Real &df, Real mu, Real bsq, Real rsq, Real rbsq, Real min_h) {
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
    };
    // }}}

    // MuFromWFunctor {{{
    class MuFromWFunctor {
      public:
      inline void operator()(Real &f, Real &df, Real mu, Real bsq, Real rsq, Real rbsq, Real W) {
        const Real musq = mu*mu;
        const Real x = 1.0/(1.0 + mu*bsq);
        const Real xsq = x*x;
        const Real rbarsq = rsq*xsq + mu*x*(1.0 + x)*rbsq;
        const Real vsq = musq*rbarsq;
        const Real dx = -bsq*xsq;
        //const Real drbarsq = rbsq*x*(1.0 + x) + (mu*rbsq + 2.0*(mu*rbsq + rsq)*x)*dx;
        const Real drbarsq = rbsq*xsq + mu*rbsq*dx + x*(rbsq + 2.0*(mu*rbsq + rsq)*dx);
        //const Real drbarsq = 2.0*rsq*dx + x*(1.0 + x)*rbsq + mu*dx*(1.0 + 2.0*x)*rbsq;
        const Real dvsq = 2.0*mu*rbarsq + musq*drbarsq;
        f = vsq + 1.0/(W*W) - 1.0;
        df = dvsq;
      }
    };
    // }}}

    // RootFunctor {{{
    class RootFunctor {
      public:
      inline Real operator()(Real mu, Real D, Real q, Real bsq, Real rsq, Real rbsq, Real *Y,
          EOS<EOSPolicy, ErrorPolicy> *const peos, Real* n, Real* T, Real* P) {
        // We need to get some utility quantities first.
        const Real x = 1.0/(1.0 + mu*bsq);
        const Real xsq = x*x;
        const Real musq = mu*mu;
        //const Real den = 1.0 + mu*bsq;
        //const Real mux = mu*x;
        //const Real muxsq = mux/den;
        //const Real rbarsq = rsq*xsq + mu*x*(1.0 + x)*rbsq;
        //const Real rbarsq = xsq*(rsq + mu*(2.0 + mu*bsq)*rbsq);
        // An alternative calculation of rbarsq that may be more accurate.
        //const Real rbarsq = rsq*xsq + (mux + muxsq)*rbsq;
        const Real rbarsq = x*(rsq*x + mu*(x + 1.0)*rbsq);
        //const Real qbar = q - 0.5*bsq - 0.5*musq*xsq*(bsq*rsq - rbsq);
        const Real qbar = q - 0.5*bsq - 0.5*musq*xsq*std::fma(bsq, rsq, -rbsq);
        const Real mb = peos->GetBaryonMass();

        // Now we can estimate the velocity.
        //const Real v_max = peos->GetMaxVelocity();
        const Real h_min = peos->GetMinimumEnthalpy();
        const Real vsq_max = std::min(rsq/(h_min*h_min + rsq), 
                                      peos->GetMaxVelocity()*peos->GetMaxVelocity());
        const Real vhatsq = std::min(musq*rbarsq, vsq_max);

        // Using the velocity estimate, predict the Lorentz factor.
        // NOTE: for extreme velocities, this alternative form of W may be more accurate:
        // Wsq = 1/(eps*(2 - eps)) = 1/(eps*(1 + v)), where eps = 1 - v.
        //const Real What = 1.0/std::sqrt(1.0 - vhatsq);
        const Real iWhat = std::sqrt(1.0 - vhatsq);

        // Now estimate the number density.
        Real rhohat = D*iWhat;
        Real nhat = rhohat/mb;
        peos->ApplyDensityLimits(nhat);

        // Estimate the energy density.
        Real eoverD = qbar - mu*rbarsq + 1.0;
        Real ehat = D*eoverD;
        peos->ApplyEnergyLimits(ehat, nhat, Y);
        //eoverD = ehat/D;

        // Now we can get an estimate of the temperature, and from that, the pressure and enthalpy.
        Real That = peos->GetTemperatureFromE(nhat, ehat, Y);
        peos->ApplyTemperatureLimits(That);
        //ehat = peos->GetEnergy(nhat, That, Y);
        Real Phat = peos->GetPressure(nhat, That, Y);
        Real hhat = peos->GetEnthalpy(nhat, That, Y);

        // Now we can get two different estimates for nu = h/W.
        Real nu_a = hhat*iWhat;
        //Real ahat = Phat / ehat;
        Real nu_b = eoverD + Phat/D;
        //Real nu_b = (1.0 + ahat)*eoverD;
        //Real nu_b = (1.0 + ahat)*eoverD;
        Real nuhat = std::max(nu_a, nu_b);

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
        std::cout << "    rbsq = " << rbsq << "\n"*/

        return mu - muhat;


      }
    };
    // }}}
  private:
    /// A constant pointer to the EOS.
    /// We make this constant because the
    /// possibility of changing the EOS
    /// during implementation seems both
    /// unlikely and dangerous.
    EOS<EOSPolicy, ErrorPolicy> *const peos;

    /// The root solver.
    NumTools::Root root;
    UpperRootFunctor UpperRoot;
    MuFromWFunctor MuFromW;
    RootFunctor RootFunction;
    
    //! \brief Check and handle the corner case for rho being too small or large.
    //
    //  Using the minimum and maximum values of rho along with some physical
    //  limitations on the velocity using S, we can predict if rho is going
    //  to violate constraints set by the EOS on how big or small it can get.
    //  We can also use these constraints to tighten the bounds on mu.
    //  
    //  \param[in,out] mul   The lower bound for mu
    //  \param[in,out] muh   The upper bound for mu
    //  \param[in]     D     The relativistic density
    //  \param[in]     bsq   The square magnitude of the magnetic field
    //  \param[in]     rsq   The square magnitude of the specific momentum S/D
    //  \param[in]     rbsq  The square of the product \f$r\cdot b\f$
    //  \param[in]     h_min The minimum enthalpy
    //
    //  \return an Error code, usually RHO_TOO_BIG, RHO_TOO_SMALL, or SUCCESS
    Error CheckDensityValid(Real& mul, Real& muh, Real D, Real bsq, Real rsq, Real rbsq, Real h_min);
  public:
    Real tol;

    /// Constructor
    PrimitiveSolver(EOS<EOSPolicy, ErrorPolicy> *eos) : peos(eos) {
      //root = NumTools::Root();
      tol = 1e-15;
      root.iterations = 30;
    }

    /// Destructor
    ~PrimitiveSolver() = default;

    //! \brief Get the primitive variables from the conserved variables.
    //
    //  \param[out]    prim  The array of primitive variables
    //  \param[in,out] cons  The array of conserved variables
    //  \param[in,out] bu    The magnetic field
    //  \param[in]     g3d   The 3x3 spatial metric
    //  \param[in]     g3u   The 3x3 inverse spatial metric
    //
    //  \return information about the solve
    SolverResult ConToPrim(Real prim[NPRIM], Real cons[NCONS], Real b[NMAG], 
                           Real g3d[NSPMETRIC], Real g3u[NSPMETRIC]);

    //! \brief Get the conserved variables from the primitive variables.
    //
    //  \param[in]    prim  The array of primitive variables
    //  \param[out]   cons  The array of conserved variables
    //  \param[in]    bu    The magnetic field
    //  \param[in]    g3d   The 3x3 spatial metric
    //
    //  \return an error code
    Error PrimToCon(Real prim[NPRIM], Real cons[NCONS], Real b[NMAG], 
                   Real g3d[NSPMETRIC]);

    /// Get the EOS used by this PrimitiveSolver.
    inline EOS<EOSPolicy, ErrorPolicy> *const GetEOS() const {
      return peos;
    }

    /// Get the root solver used by this PrimitiveSolver.
    inline NumTools::Root& GetRootSolver() {
      return root;
    }

    //! \brief Do failure response and adjust conserved variables if necessary.
    //
    //  Note that in ConToPrim, the error policy dictates whether or not we
    //  should adjust the conserved variables if the primitive variables are
    //  floored. That may appear to be the case here, too, because of the bool
    //  returned by DoFailureResponse. However, DoFailureResponse simply tells
    //  us whether or not the primitives were adjusted in the first place, not
    //  whether or not we should adjust the conserved variables. Thus, if the
    //  primitive variables are modified as part of the error response, the
    //  conserved variables are *always* altered. The reasoning here is that
    //  flooring after a primitive solver indicates a physical state which is
    //  just slightly out of bounds. Depending on how the conserved variables
    //  are to be used afterward, the user may not find it necessary to rescale
    //  them. On the other hand, a failure mode generally indicates that the
    //  state itself is unphysical, so any modification to the primitives
    //  also requires a modification to the conserved variables.
    //
    //  \param[in,out] prim  The array of primitive variables
    //  \param[in,out] cons  The array of conserved variables
    //  \param[in,out] bu    The magnetic field
    //  \param[in]     g3d   The 3x3 spatial metric
    void HandleFailure(Real prim[NPRIM], Real cons[NCONS], Real bu[NMAG],
                       Real g3d[NSPMETRIC]) {
      bool result = peos->DoFailureResponse(prim);
      if (result) {
        PrimToCon(prim, cons, bu, g3d);
      }
    }
};

// CheckDensityValid {{{
template<typename EOSPolicy, typename ErrorPolicy>
inline Error PrimitiveSolver<EOSPolicy, ErrorPolicy>::CheckDensityValid(Real& mul, Real& muh, Real D, 
      Real bsq, Real rsq, Real rbsq, Real h_min) {
  // There are a few things considered:
  // 1. If D > rho_max, we need to make sure that W isn't too large.
  //    W_max can be estimated by considering the zero-field limit
  //    of S^2/D^2 if h = h_min.
  //    - If W is larger than W_max, then rho is just too big.
  //    - Otherwise, we can bound mu by using W to do a root solve
  //      for mu.
  // 2. If D < W_max*rho_min, then we need to make sure W isn't less
  //    than 1.
  //    - If W is less than 1, it means rho is actually smaller than
  //      rho_min.
  //    - Otherwise, we can bound mu by using W to do a root solve
  //      for mu.
  Real W_max = std::sqrt(1.0 + rsq/(h_min*h_min));
  Real rho_max = peos->GetMaximumDensity()*peos->GetBaryonMass();
  Real rho_min = peos->GetMinimumDensity()*peos->GetBaryonMass();
  if (D > rho_max) {
    Real W = D/rho_max;
    Real f, df;
    MuFromW(f, df, muh, bsq, rsq, rbsq, W);
    if (f <= 0) {
      // W is not physical, so rho must be larger than rho_max.
      return Error::RHO_TOO_BIG;
    }
    else {
      MuFromW(f, df, mul, bsq, rsq, rbsq, W);
      if (f < 0) {
        Real mu;
        Real mulc = mul;
        Real muhc = muh;
        // We can tighten up the bounds for mul.
        // The derivative is zero at mu = 0, so we perturb it slightly.
        /*if (mu <= root.tol) {
          mu += root.tol;
        }*/
        bool result = root.NewtonSafe(MuFromW, mulc, muhc, mu, 1e-10, bsq, rsq, rbsq, W);
        if (!result) {
          return Error::BRACKETING_FAILED;
        }
        mul = (mu > mul) ? mu : mul;
      }
    }
  }
  if (D < W_max*rho_min) {
    Real W = D/rho_min;
    Real f, df;
    MuFromW(f, df, mul, bsq, rsq, rbsq, W);
    if (f >= 0) {
      // W is not physical, so rho must be smaller than rho_min.
      return Error::RHO_TOO_SMALL;
    }
    else {
      MuFromW(f, df, muh, bsq, rsq, rbsq, W);
      if (f > 0) {
        Real mu = muh;
        Real mulc = mul;
        Real muhc = muh;
        // We can tighten up the bounds for muh.
        bool result = root.NewtonSafe(MuFromW, mulc, muhc, mu, 1e-10, bsq, rsq, rbsq, W);
        if (!result) {
          return Error::BRACKETING_FAILED;
        }
        muh = (mu < muh) ? mu : muh;
      }
    }
  }
  return Error::SUCCESS;
}
// }}}

// ConToPrim {{{
template<typename EOSPolicy, typename ErrorPolicy>
inline SolverResult PrimitiveSolver<EOSPolicy, ErrorPolicy>::ConToPrim(Real prim[NPRIM], 
      Real cons[NCONS], Real b[NMAG], Real g3d[NSPMETRIC], Real g3u[NSPMETRIC]) {

  SolverResult solver_result{Error::SUCCESS, 0, false, false, false};

  // Extract the undensitized conserved variables.
  Real D = cons[IDN];
  Real S_d[3] = {cons[IM1], cons[IM2], cons[IM3]};
  Real tau = cons[IEN];
  Real B_u[3] = {b[IB1], b[IB2], b[IB3]};
  // Extract the particle fractions.
  const int n_species = peos->GetNSpecies();
  Real Y[MAX_SPECIES] = {0.0};
  for (int s = 0; s < n_species; s++) {
    Y[s] = cons[IYD + s]/cons[IDN];
  }
  // Apply limits to Y to ensure a physical state
  bool Y_adjusted = peos->ApplySpeciesLimits(Y);

  // Check the conserved variables for consistency and do whatever
  // the EOSPolicy wants us to.
  bool floored = peos->ApplyConservedFloor(D, S_d, tau, Y, SquareVector(B_u, g3d));
  solver_result.cons_floor = floored;
  if (floored && peos->IsConservedFlooringFailure()) {
    HandleFailure(prim, cons, b, g3d);
    solver_result.error = Error::CONS_FLOOR;
    return solver_result;
  }
  // If a floor is applied or Y is adjusted, we need to propagate the changes back to
  // DYe.
  if (floored || Y_adjusted) {
    for (int s = 0; s < n_species; s++) {
      cons[IYD + s] = D*Y[s];
    }
  }

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
    HandleFailure(prim, cons, b, g3d);
    solver_result.error = Error::NANS_IN_CONS;
    return solver_result;
  }
  // We have to check the particle fractions separately.
  for (int s = 0; s < n_species; s++) {
    if (!std::isfinite(Y[s])) {
      HandleFailure(prim, cons, b, g3d);
      solver_result.error = Error::NANS_IN_CONS;
      return solver_result;
    }
  }

  // Make sure that the magnetic field is physical.
  Error error = peos->DoMagnetizationResponse(bsqr, b_u);
  if (error == Error::MAG_TOO_BIG) {
    HandleFailure(prim, cons, b, g3d);
    solver_result.error = Error::MAG_TOO_BIG;
    return solver_result;
  }
  else if (error == Error::CONS_ADJUSTED) {
    solver_result.cons_adjusted = true;
    // If b_u is rescaled, we also need to adjust D, which means we'll
    // have to adjust all our other rescalings, too.
    Real Bsq = SquareVector(B_u, g3d);
    D = Bsq/bsqr;
    r_d[0] = S_d[0]/D; r_d[1] = S_d[1]/D; r_d[2] = S_d[2]/D;
    RaiseForm(r_u, r_d, g3d);
    // We need to recalculate rb if b_u is rescaled.
    rb = Contract(b_u, r_d);
    rbsqr = rb*rb;
    q = tau/D;
    rsqr = Contract(r_d, r_u);
  }
  
  // Bracket the root.
  Real min_h = peos->GetMinimumEnthalpy();
  Real mul = 0.0;
  Real muh = 1.0/min_h;
  // Check if a tighter upper bound exists.
  if (rsqr > min_h*min_h) {
    Real mu = 0.0;
    // We don't need the bound to be that tight, so we reduce
    // the accuracy of the root solve for speed reasons.
    Real mulc = mul;
    Real mulh = muh;
    bool result = root.NewtonSafe(UpperRoot, mulc, mulh, mu, 1e-10,
                                  bsqr, rsqr, rbsqr, min_h);
    // Scream if the bracketing failed.
    if (!result) {
      HandleFailure(prim, cons, b, g3d);
      solver_result.error = Error::BRACKETING_FAILED;
      return solver_result;
    }
    else {
      // To avoid problems with the case where the root and the upper bound collide,
      // we will perturb the bound slightly upward.
      // TODO: Is there a more rigorous way of treating this?
      muh = mu*(1. + 1e-10);
    }
  }

  // Check the corner case where the density is outside the permitted
  // bounds according to the ErrorPolicy.
  error = CheckDensityValid(mul, muh, D, bsqr, rsqr, rbsqr, min_h);
  // TODO: This is probably something that should be handled by the ErrorPolicy.
  if (error != Error::SUCCESS) {
    HandleFailure(prim, cons, b, g3d);
    solver_result.error = error;
    return solver_result;
  }

  
  // Do the root solve.
  // TODO: This should be done with something like TOMS748 once it's
  // available.
  Real n, P, T, mu;
  NumTools::Root::RootResult result = root.FalsePosition(RootFunction, mul, muh, mu, tol, D, q, bsqr, rsqr, rbsqr, Y, peos, &n, &T, &P);
  //NumTools::Root::RootResult result = root.ITP(RootFunction, mul, muh, mu, D, q, bsqr, rsqr, rbsqr, Y, peos, &n, &T, &P);
  solver_result.iterations = result.iterations;
  if (!result.success) {
    // It may be the case that the result isn't great, but it's still valid. In this case,
    // we just warn the user about convergence being poor.
    if (result.iterations == root.iterations && 
        result.err < peos->GetFailureTolerance()) {
      solver_result.error = Error::SLOW_CONVERGENCE;
    }
    else {
      HandleFailure(prim, cons, b, g3d);
      if (result.bracketed) {
        solver_result.error = Error::NO_SOLUTION;
      } else {
        solver_result.error = Error::BRACKETING_FAILED;
      }
      return solver_result;
    }
  }

  // Retrieve the primitive variables.
  Real rho = n*peos->GetBaryonMass();
  Real rbmu = rb*mu;
  Real W = D/rho;
  Real Wmux = W*mu/(1.0 + mu*bsqr);
  // Before we retrieve the velocity, we need to raise S.
  Real S_u[3] = {0.0};
  RaiseForm(S_u, S_d, g3u);
  // Now we can get Wv.
  Real Wv_u[3] = {0.0};
  Wv_u[0] = Wmux*(r_u[0] + rbmu*b_u[0]);
  Wv_u[1] = Wmux*(r_u[1] + rbmu*b_u[1]);
  Wv_u[2] = Wmux*(r_u[2] + rbmu*b_u[2]);
  
  // Apply the flooring policy to the primitive variables.
  floored = peos->ApplyPrimitiveFloor(n, Wv_u, P, T, Y);
  solver_result.prim_floor = floored;
  if (floored && peos->IsPrimitiveFlooringFailure()) {
    HandleFailure(prim, cons, b, g3d);
    solver_result.error = Error::PRIM_FLOOR;
    return solver_result;
  }
  solver_result.cons_adjusted = solver_result.cons_adjusted || floored ||
                                solver_result.cons_floor || Y_adjusted;

  prim[IDN] = n;
  prim[IPR] = P;
  prim[ITM] = T;
  prim[IVX] = Wv_u[0];
  prim[IVY] = Wv_u[1];
  prim[IVZ] = Wv_u[2];
  for (int s = 0; s < n_species; s++) {
    prim[IYF + s] = Y[s];
  }

  // If we floored the primitive variables, we should check
  // if the EOS wants us to adjust the conserved variables back
  // in bounds. If that's the case, then we'll do it.
  if (solver_result.cons_adjusted && peos->KeepPrimAndConConsistent()) {
    PrimToCon(prim, cons, b, g3d);
  }
  else {
    solver_result.cons_adjusted = false;
  }

  return solver_result;
}
// }}}

// PrimToCon {{{
template<typename EOSPolicy, typename ErrorPolicy>
inline Error PrimitiveSolver<EOSPolicy, ErrorPolicy>::PrimToCon(Real prim[NPRIM], Real cons[NCONS],
      Real bu[NMAG], Real g3d[NMETRIC]) {
  // Extract the primitive variables
  const Real &n = prim[IDN]; // number density
  const Real Wv_u[3] = {prim[IVX], prim[IVY], prim[IVZ]};
  const Real &p   = prim[IPR]; // pressure
  const Real &t   = prim[ITM]; // temperature
  const Real B_u[3] = {bu[IB1], bu[IB2], bu[IB3]};
  const int n_species = peos->GetNSpecies();
  Real Y[MAX_SPECIES] = {0.0};
  for (int s = 0; s < n_species; s++) {
    Y[s] = prim[IYF + s];
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
  Real &D = cons[IDN]; // Relativistic density
  Real &Sx = cons[IM1]; // Relativistic momentum density (x)
  Real &Sy = cons[IM2]; // Relativistic momentum density (y)
  Real &Sz = cons[IM3]; // Relativistic momentum density (z)
  Real &tau = cons[IEN]; // Relativistic energy - D

  // Set the conserved quantities.
  // Total enthalpy density
  Real H = n*peos->GetEnthalpy(n, t, Y)*mb;
  Real HWsq = H*Wsq;
  D = n*mb*W;
  for (int s = 0; s < n_species; s++) {
    cons[IYD + s]= D*Y[s];
  }
  Real HWsqpb = HWsq + Bsq;
  Sx = (HWsqpb*v_d[0] - Bv*B_d[0]);
  Sy = (HWsqpb*v_d[1] - Bv*B_d[1]);
  Sz = (HWsqpb*v_d[2] - Bv*B_d[2]);
  tau = (HWsqpb - p - 0.5*(Bv*Bv + Bsq/Wsq)) - D;

  return Error::SUCCESS;
}
// }}}

} // namespace
#endif
