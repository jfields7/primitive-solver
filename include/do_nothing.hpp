#ifndef DO_NOTHING_HPP
#define DO_NOTHING_HPP
//! \file do_nothing.hpp
//  \brief Describes an error policy that basically does nothing.
//
//  Since this policy does nothing, it should really only be used for testing
//  purposes.

#include <ps_types.hpp>
#include <error_policy_interface.hpp>

namespace Primitive {

class DoNothing : public ErrorPolicyInterface {
  protected:
    DoNothing() {
      fail_conserved_floor = false;
      fail_primitive_floor = false;
      adjust_conserved = false;
    }

    bool PrimitiveFloor(Real& n, Real v[3], Real& t) {return false;}
    bool ConservedFloor(Real& D, Real Sd[3], Real& tau, Real tau_floor) {return false;}
};

}; // namespace

#endif
