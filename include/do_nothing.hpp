#ifndef DO_NOTHING_HPP
#define DO_NOTHING_HPP
//! \file do_nothing.hpp
//  \brief Describes an error policy that basically does nothing.
//
//  Since this policy does nothing, it should really only be used for testing
//  purposes.

#include <ps_types.hpp>

namespace Primitive {

class DoNothing {
  protected:
    Real n_atm;
    Real T_atm;

    bool PrimitiveFloor(Real& n, Real v[3], Real& t) {return false;}
    bool ConservedFloor(Real& D, Real Sd[3], Real& tau, Real tau_floor) {return false;}
};

}; // namespace

#endif
