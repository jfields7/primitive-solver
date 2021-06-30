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

namespace Primitive {

class ResetFloor {
  protected:
    Real n_atm;
    Real T_atm;
    Real v_max;

    /// Constructor
    ResetFloor();

    /// Floor for primitive variables
    bool PrimitiveFloor(Real& n, Real v[3], Real& t);

    /// Floor for conserved variables
    bool ConservedFloor(Real& D, Real Sd[3], Real& tau, Real tau_floor);
};

}; // namespace

#endif
