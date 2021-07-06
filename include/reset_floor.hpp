#ifndef RESET_FLOOR_HPP
#define RESET_FLOOR_HPP
//! \file reset_floor.hpp
//  \brief Describes an error floor that simply resets nonphysical values.
//
//  If the density or pressure fall below the atmosphere, they get floored.
//  We impose similar limits for D and tau. If the density is floored,
//  the velocity is zeroed out and the pressure is also reset to the floor.
//  If the pressure is floored, all other quantities are ignored.

#include <ps_types.hpp>
#include <error_policy_interface.hpp>

namespace Primitive {

class ResetFloor : public ErrorPolicyInterface {
  protected:
    /// Constructor
    ResetFloor();

    /// Floor for primitive variables
    bool PrimitiveFloor(Real& n, Real v[3], Real& t);

    /// Floor for conserved variables
    bool ConservedFloor(Real& D, Real Sd[3], Real& tau, Real tau_floor);

  public:
    /// Set the failure mode for conserved flooring
    inline void SetConservedFloorFailure(bool failure) {
      fail_conserved_floor = failure;
    }

    /// Set the failure mode for primitive flooring
    inline void SetPrimitiveFloorFailure(bool failure) {
      fail_primitive_floor = failure;
    }

    /// Set whether or not it's okay to adjust the conserved variables.
    inline void SetAdjustConserved(bool adjust) {
      adjust_conserved = adjust;
    }
};

}; // namespace

#endif
