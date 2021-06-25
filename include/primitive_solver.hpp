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
  public:
    /// Constructor
    PrimitiveSolver(EOS<EOSPolicy, ErrorPolicy> *eos) : peos(eos) {
      
    }

    /// Destructor
    ~PrimitiveSolver() = default;

    //! \brief Get the primitive variables from the conserved variables.
    //
    //  \param[out]    prim  The array of primitive variables
    //  \param[in,out] cons  The array of conserved variables
    //  \param[in,out] bu    The magnetic field
    //  \param[in]     gd    The metric
    //  \param[in]     gu    The inverse metric
    //  \param[in]     i,j,k The position in the array
    //
    //  \return success or failure
    bool ConToPrim(AthenaArray<Real>& prim, AthenaArray<Real>& cons,
                   AthenaArray<Real>& b, AthenaArray<Real>& gu
                   int i, int j, int k);

    //! \brief Get the conserved variables from the primitive variables.
    //
    //  \param[in]    prim  The array of primitive variables
    //  \param[out]   cons  The array of conserved variables
    //  \param[in]    bu    The magnetic field
    //  \param[in]    detg  The determinant of the metric
    //  \param[in]    i,j,k The position in the array
};

template<typename EOSPolicy, typename ErrorPolicy>
bool PrimitiveSolver<EOSPolicy, ErrorPolicy>::ConToPrim(AthenaArray<Real>& prim,
       AthenaArray<Real>& cons, AthenaArray<Real>& b, AthenaArray<Real>& gd,
       AthenaArray<Real>& gu, int i, int j, int k){
  

  // Let's extract the undensitized conserved variables.
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


#endif
