#ifndef PRIMITIVE_SOLVER_HPP
#define PRIMITIVE_SOLVER_HPP

//! \file primitive_solver.hpp
//  \brief Declares PrimitiveSolver class.
//
//  PrimitiveSolver contains all the infrastructure for the inversion
//  procedure from conserved to primitive variables in GRMHD. This
//  particular implementation is based on the solver described in
//  Kastaun et al., Phys. Rev. D 103, 023018 (2021).

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

    
};

#endif
