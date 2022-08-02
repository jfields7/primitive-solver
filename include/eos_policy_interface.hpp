#ifndef EOS_POLICY_INTERFACE_HPP
#define EOS_POLICY_INTERFACE_HPP

//! \file eos_policy_interface.hpp
//  \brief Defines a class that provides all the common class
//         member variables needed by an EOSPolicy.

#include <ps_types.hpp>

namespace Primitive {

class EOSPolicyInterface {
  protected:
    EOSPolicyInterface() = default;
    ~EOSPolicyInterface() = default;

    /// Number of particle species
    int n_species;
    /// Baryon mass
    Real mb;
    /// maximum number density
    Real max_n;
    /// minimum number density
    Real min_n;
    /// maximum temperature
    Real max_T;
    /// minimum temperature
    Real min_T;
};

} // namespace

#endif
