#ifndef IDEAL_GAS_HPP
#define IDEAL_GAS_HPP

//! \file idealgas.hpp
//  \brief Defines an ideal gas equation of state.

#include <ps_types.hpp>
#include <eos_policy_interface.hpp>

namespace Primitive{

class IdealGas : public EOSPolicyInterface{
  protected:
    /// Adiabatic index
    Real gamma;
    Real gammam1;

    /// Constructor
    IdealGas();

    /// Calculate the temperature using the ideal gas law.
    Real TemperatureFromE(Real n, Real e, Real *Y);

    /// Calculate the temperature using the ideal gas law.
    Real TemperatureFromP(Real n, Real p, Real *Y);

    /// Calculate the energy density using the ideal gas law.
    Real Energy(Real n, Real T, Real *Y);

    /// Calculate the pressure using the ideal gas law.
    Real Pressure(Real n, Real T, Real *Y);

    /// Calculate the entropy per baryon using the ideal gas law.
    [[ noreturn ]]
    Real Entropy(Real n, Real T, Real *Y);

    /// Calculate the enthalpy per baryon using the ideal gas law.
    Real Enthalpy(Real n, Real T, Real *Y);

    /// Get the minimum enthalpy per baryon according to the ideal gas law.
    Real MinimumEnthalpy();

    /// Calculate the sound speed for an ideal gas.
    Real SoundSpeed(Real n, Real T, Real *Y);

    /// Calculate the energy per baryon (NOT energy per mass!)
    Real SpecificEnergy(Real n, Real T, Real *Y);

  public:
    /// Set the adiabatic index for the ideal gas. 
    /// The range \f$1 < \gamma < 1\f$ is imposed. The lower
    /// constraint ensures that enthalpy is finite, and the upper
    /// bound keeps the sound speed causal.
    inline void SetGamma(Real g) {
      gamma = (g <= 1.0) ? 1.00001 : ((g >= 2.0) ? 2.00001 : g);
      gammam1 = gamma - 1.0;
    }

    /// Get the adiabatic index.
    inline Real GetGamma() const {
      return gamma;
    }

    /// Set the baryon mass
    inline void SetBaryonMass(Real m) {
      mb = m;
    }

    /// Set the number of species. Throw an exception if
    /// the number of species is invalid.
    void SetNSpecies(int n);
};

}; // namespace

#endif
