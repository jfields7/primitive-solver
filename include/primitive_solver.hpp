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
    const int n_species;
  public:
    /// Constructor
    PrimitiveSolver(EOS<EOSPolicy, ErrorPolicy> *eos) : peos(eos), n_species(eos->GetNSpecies()) {
      
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

// ConToPrim {{{
template<typename EOSPolicy, typename ErrorPolicy>
bool PrimitiveSolver<EOSPolicy, ErrorPolicy>::ConToPrim(AthenaArray<Real>& prim,
       AthenaArray<Real>& cons, AthenaArray<Real>& b, AthenaArray<Real>& gd,
       AthenaArray<Real>& gu, int i, int j, int k) {

  // Extract the 3-metric and inverse 3-metric.
  const Real g3d[NSPMETRIC] = {gd(I11, i), gd(I12, i), gd(I13, i),
                               gd(I22, i), gd(I23, i), gd(I33, i)};
  const Real ialphasq = -gu(I00, i);
  const Real alphasq = -1.0/ialphasq; // Lapse squared
  const Real beta_u[3] = {gu(I01, i)*alphasq, gu(I02, i)*alphasq, gu(I03, i)*alphasq}; // Shift squared
  const Real g3u[NSPMETRIC] = {gu(I11, i) + beta_u[0]*beta_u[0]*ialphasq,
                               gu(I12, i) + beta_u[0]*beta_u[1]*ialphasq,
                               gu(I13, i) + beta_u[0]*beta_u[2]*ialphasq,
                               gu(I22, i) + beta_u[2]*beta_u[2]*ialphasq,
                               gu(I23, i) + beta_u[2]*beta_u[3]*ialphasq,
                               gu(I33, i) + beta_u[3]*beta_u[3]*ialphasq};

  // Get the inverse volume element of the 3-metric.
  Real isdetg = 1.0/std::sqrt(GetDeterminant(g3d));

  // Extract the undensitized conserved variables.
  Real D = cons(IDN, k, j, i)*isdetg;
  Real S_d[3] = {cons(IM1, k, j, i)*isdetg, cons(IM2, k, j, i)*isdetg, cons(IM3, k, j, i)*isdetg};
  Real tau = cons(IEN, k, j, i)*isdetg;
  // FIXME: Confirm that the magnetic field is densitized as well.
  Real B_u[3] = {b(IB1, k, j, i)*isdetg, b(IB2, k, j, i)*isdetg, b(IB3, k, j, i)*isdetg};
  // Extract the particle fractions.
  Real Y[n_species] = {0.0};
  for (int s = 0; s < n_species; s++) {
    Y[s] = cons(IYD + s, k, j, i)/cons(IDN, k, j, i);
  }

  // If D is below the atmosphere, we need to do whatever
  // the EOSPolicy wants us to do.

  // Calculate some utility quantities.
  const Real b_u[3] = {B_u[0]/D, B_u[1]/D, B_u[2]/D};
  const Real r_d[3] = {S_d[0]/D, S_d[1]/D, S_d[2]/D};
  const Real r_u[3];
  RaiseForm(r_u, r_d, g3u);
  const Real rsqr   = Contract(r_u, r_d);
  const Real rb     = Contract(b_u, r_d);
  const Real rbsqr  = rb*rb;
  const Real bsqr   = SquareVector(b_u, g3d);
  const Real q      = tau/D;

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

  // Make sure that the magnetic field is physical.

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
  Sy = sdetg*(HWsqpb*v_d[2] - Bv*B_d[2]);
  tau = sdetg*(HWsqpb - p - 0.5*(Bv*Bv + Bsq/Wsq));
}
// }}}

}; // namespace
#endif
