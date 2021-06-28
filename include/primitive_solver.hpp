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
};

template<typename EOSPolicy, typename ErrorPolicy>
bool PrimitiveSolver<EOSPolicy, ErrorPolicy>::ConToPrim(AthenaArray<Real>& prim,
       AthenaArray<Real>& cons, AthenaArray<Real>& b, AthenaArray<Real>& gd,
       AthenaArray<Real>& gu, int i, int j, int k) {
  

  // Extract the undensitized conserved variables.
  Real D = cons(IDN, k, j, i)*isdetg;
  Real Sd[3] = {cons(IM1, k, j, i)*isdetg, cons(IM2, k, j, i)*isdetg, cons(IM3, k, j, i)*isdetg};
  Real tau = cons(IEN, k, j, i)*isdetg;
  // FIXME: Confirm that the magnetic field is densitized as well.
  Real Bu[3] = {b(IB1, k, j, i)*isdetg, b(IB2, k, j, i)*isdetg, b(IB3, k, j, i)*isdetg};

  // Make sure nothing is a NaN at this point.
  if(!isfinite(D) || !isfinite(Sd[0]) || !isfinite(Sd[1]) || 
     !isfinite(Sd[2]) || !isfinite(Bu[0]) || !isfinite(Bu[1]) ||
     !isfinite(Bu[2]) || !isfinite(Bu[2])) {
    return false;
  }

}

template<typename EOSPolicy, typename ErrorPolicy>
bool PrimitiveSolver<EOSPolicy, ErrorPolicy>::PrimToCon(AthenaArray<Real>& prim,
       AthenaArray<Real>& cons, AthenaArray<Real>& bu, AthenaArray<Real>& gd,
       AthenaArray<Real>& gu, AthenaArray<Real>& alpha, AthenaArray<Real>& beta,
       int i, int j, int k) {
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

  // Calculate the Lorentz factor. Note: Athena passes in Wv, not v.
  Real usq = g3d[S11]*u1*u1 + g3d[S22]*u2*u2 + g3d[S33]*u3*u3 + 2.0*((g3d[S12]*u2 + g3d[S13]*u3)*u1 + g3d[S23]*u2*u3);
  Real W = std::sqrt(1.0 + usq);

  // Get the four-velocity.
  Real u0 = W/alpha;
  Real u1 = uu1 - W * alpha * gu(I01, i);
  Real u2 = uu2 - W * alpha * gu(I02, i);
  Real u3 = uu3 - W * alpha * gu(I03, i);
  // Lower the four-velocity.
  Real u_0 = gd(I00,i)*u0 + gd(I01,i)*u1 + gd(I02,i)*u2 + gd(I03, i)*u3;
  Real u_1 = gd(I01,i)*u0 + gd(I11,i)*u1 + gd(I12,i)*u2 + gd(I13, i)*u3;
  Real u_2 = gd(I02,i)*u0 + gd(I12,i)*u1 + gd(I22,i)*u2 + gd(I23, i)*u3;
  Real u_3 = gd(I03,i)*u0 + gd(I13,i)*u1 + gd(I23,i)*u2 + gd(I33, i)*u3;

  // Get the four-magnetic field.
  Real b0 = gd(I01,i)*u0*bb1 + gd(I02,i)*u0*bb2 + gd(I03,i)*u0*bb3
          + gd(I11,i)*u1*bb1 + gd(I12,i)*u1*bb2 + gd(I13,i)*u1*bb3
          + gd(I12,i)*u2*bb1 + gd(I22,i)*u2*bb2 + gd(I23,i)*u2*bb3
          + gd(I13,i)*u3*bb1 + gd(I23,i)*u3*bb2 + gd(I33,i)*u3*bb3;
  Real b1 = (bb1 + b0*u1)/u0;
  Real b2 = (bb2 + b0*u2)/u0;
  Real b3 = (bb3 + b0*u3)/u0;
  // Lower the magnetic field.
  Real b_0 = gd(I00,i)*b0 + gd(I01,i)*b1 + gd(I02,i)*b2 + gd(I03,I)*b3;
  Real b_1 = gd(I01,i)*b0 + gd(I11,i)*b1 + gd(I12,i)*b2 + gd(I13,I)*b3;
  Real b_2 = gd(I02,i)*b0 + gd(I12,i)*b1 + gd(I22,i)*b2 + gd(I23,I)*b3;
  Real b_3 = gd(I03,i)*b0 + gd(I13,i)*b1 + gd(I23,i)*b2 + gd(I33,I)*b3;
  // Get the square of the field.
  Real b_sq = b0*b_0 + b1*b_1 + b2*b_2 + b3*b_3;

  // Some utility quantities that will be helpful.
  const Real mb = eos->GetBaryonMass();
  Real Wsq = W*W;

  // Extract the conserved variables
  Real &D = cons(IDN, k, j ,i); // Relativistic density
  Real &Sx = cons(IM1, k, j, i); // Relativistic momentum density (x)
  Real &Sy = cons(IM2, k, j, i); // Relativistic momentum density (y)
  Real &Sz = cons(IM3, k, j, i); // Relativistic momentum density (z)
  Real &tau = cons(IEN, k, j, i); // Relativistic energy - D

  // Set the conserved quantities.
  // Total enthalpy density
  Real H = rho*eos->GetEnthalpy(rho/mb, t, Y)/mb;
  Real HW = H*W;
  D = sdetg*rho*W;
  for (int s = 0; s < n_species; s++) {
    cons(IYD + s, k, j, i) = D*Y[s];
  }
  Sx = sdetg*HW*u_1;
  Sy = sdetg*HW*u_2;
  Sz = sdetg*HW*u_3;
  tau = sdetg*(HW*W - p - D);
}


#endif
