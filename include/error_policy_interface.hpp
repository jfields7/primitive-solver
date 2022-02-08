#ifndef ERROR_POLICY_INTERFACE_HPP
#define ERROR_POLICY_INTERFACE_HPP
//! \file error_policy_interface.hpp
//  \brief Defines a class that provides all the basic members
//         needed by an ErrorPolicy.
//
//  It cannot be instantiated and, in fact, has no purpose in
//  being instantiated. It literally just provides member
//  variables for an ErrorPolicy;

namespace Primitive {

class ErrorPolicyInterface {
  protected:
    ErrorPolicyInterface() = default;
    ~ErrorPolicyInterface() = default;

    Real n_atm;
    Real n_threshold;
    Real p_atm;
    Real v_max;
    Real max_bsq;
    bool fail_conserved_floor;
    bool fail_primitive_floor;
    bool adjust_conserved;
};

} // namespace

#endif
