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

template<typename EOSPolicy, typename ErrorPolicy>
bool PrimitiveSolver<EOSPolicy, ErrorPolicy>::ConToPrim(AthenaArray<Real>& prim,
       AthenaArray<Real>& cons, AthenaArray<Real>& b, AthenaArray<Real>& gd,
       AthenaArray<Real>& gu, int i, int j, int k) {
  

  // Extract the undensitized conserved variables.
  /*Real D = cons(IDN, k, j, i)*isdetg;
  Real Sd[3] = {cons(IM1, k, j, i)*isdetg, cons(IM2, k, j, i)*isdetg, cons(IM3, k, j, i)*isdetg};
  Real tau = cons(IEN, k, j, i)*isdetg;
  // FIXME: Confirm that the magnetic field is densitized as well.
  Real Bu[3] = {b(IB1, k, j, i)*isdetg, b(IB2, k, j, i)*isdetg, b(IB3, k, j, i)*isdetg};

  // Make sure nothing is a NaN at this point.
  if(!isfinite(D) || !isfinite(Sd[0]) || !isfinite(Sd[1]) || 
     !isfinite(Sd[2]) || !isfinite(Bu[0]) || !isfinite(Bu[1]) ||
     !isfinite(Bu[2]) || !isfinite(Bu[2])) {
    return false;
  }*/
  return true;
}

// PrimToCon {{{
template<typename EOSPolicy, typename ErrorPolicy>
bool PrimitiveSolver<EOSPolicy, ErrorPolicy>::PrimToCon(AthenaArray<Real>& prim,
       AthenaArray<Real>& cons, AthenaArray<Real>& bu, AthenaArray<Real>& gd,
       AthenaArray<Real>& gu,int i, int j, int k) {
  // Extract the three metric.
  Real g3d[NSPMETRIC] = {gd(I11, i), gd(I12, i), gd(I13, i),
                        gd(I22, i), gd(I23, i), gd(I33, i)};
  Real alpha = std::sqrt(-1.0/gu(I00, i));
  
  // Get the volume element of the 3-metric.
  Real sdetg = std::sqrt(g3d[S11]*g3d[S22]*g3d[S33] + 2.0*g3d[S12]*g3d[S13]*g3d[S23] - 
              (g3d[S11]*g3d[S23]*g3d[S23] + g3d[S12]*g3d[S12]*g3d[S33] + g3d[S13]*g3d[S13]*g3d[S22]));

  // Extract the primitive variables
  const Real &rho = prim(IDN, k, j, i); // rest-mass density
  const Real &uu1 = prim(IVX, k, j, i); // Wvx
  const Real &uu2 = prim(IVY, k, j, i); // Wvy
  const Real &uu3 = prim(IVZ, k, j, i); // Wvz
  const Real &p   = prim(IPR, k, j, i); // pressure
  const Real &t   = prim(ITM, k, j, i); // temperature
  const Real &bb1 = bu(IB1, k, j, i);   // B^x
  const Real &bb2 = bu(IB2, k, j, i);   // B^y
  const Real &bb3 = bu(IB3, k, j, i);   // B^z
  Real Y[n_species] = {0.0};
  for (int s = 0; s < n_species; s++) {
    Y[s] = prim(IYF + s, k, j, i);
  }

  // Note that Athena passes in Wv, not v.
  // Lower u.
  Real uu_1 = g3d[S11]*uu1 + g3d[S12]*uu2 + g3d[S13]*uu3;
  Real uu_2 = g3d[S12]*uu1 + g3d[S22]*uu2 + g3d[S23]*uu3;
  Real uu_3 = g3d[S13]*uu1 + g3d[S23]*uu2 + g3d[S33]*uu3;
  Real usq = uu_1*uu1 + uu_2*uu2 + uu_3*uu3;
  Real Wsq = 1.0 + usq;
  Real W = sqrt(Wsq);
  // Get the 3-velocity.
  Real v_1 = uu_1/W;
  Real v_2 = uu_2/W;
  Real v_3 = uu_3/W;

  // For the magnetic field contribution, we need to find
  // B_i, B^2, and B^i*v_i.
  Real bb_1 = g3d[S11]*bb1 + g3d[S12]*bb2 + g3d[S13]*bb3;
  Real bb_2 = g3d[S12]*bb1 + g3d[S22]*bb2 + g3d[S23]*bb3;
  Real bb_3 = g3d[S13]*bb1 + g3d[S23]*bb2 + g3d[S33]*bb3;
  Real Bsq = bb1*bb_1 + bb2*bb_2 + bb3*bb_3;
  Real Bv = bb1*v_1 + bb2*v_2 + bb3*v_3;

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
  Sx = sdetg*(HWsqpb*v_1 - Bv*bb_1);
  Sy = sdetg*(HWsqpb*v_2 - Bv*bb_2);
  Sy = sdetg*(HWsqpb*v_3 - Bv*bb_3);
  tau = sdetg*(HWsqpb - p - 0.5*(Bv*Bv + Bsq/Wsq));
}
// }}}

#endif
