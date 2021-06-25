//! \file reset_floor.cpp
//  \brief Implementation of the ResetFloor policy

#include <reset_floor.hpp>

/// Constructor
ResetFloor::ResetFloor() {
}

/// Floor for the primitive variables
bool ResetFloor::PrimitiveFloor(Real& n, Real v[3], Real& t) {
  if (n < n_atm) {
    n = n_atm;
    v[0] = 0.0;
    v[1] = 0.0;
    v[2] = 0.0;
    t = T_atm;
    return true;
  }
  else if (t < T_atm) {
    t = T_atm;
    return true;
  }
  return false;
}

/// Floor for the conserved variables
/// FIXME: This currently sets tau to the pressure atmosphere,
/// but it needs to be replaced with something which is
/// thermodynamically consistent.
bool ResetFloor::ConservedFloor(Real& D, Real Sd[3], Real& tau, Real tau_floor) {
  if (D < n_atm) {
    D = n_atm;
    Sd[0] = 0.0;
    Sd[1] = 0.0;
    Sd[2] = 0.0;
    tau = tau_floor;
    return true;
  }
  else if (tau < tau_floor) {
    tau = tau_floor;
    return true;
  }
  return false;
}
