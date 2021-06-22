#ifndef EOS_HPP
#define EOS_HPP

//! \file eos.hpp
//  \brief Defines an equation of state.
//
//  EOS is effectively an interface that describes how to create an
//  equation of state. It must be instantiated with an object implementing
//  the following functions:
//    Real Temperature(Real n, Real e, Real *Y)
//    Real Energy(Real n, Real T, Real *Y)
//    Real Pressure(Real n, Real T, Real *Y)
//    Real Entropy(Real n, Real T, Real *Y)
//    Real Enthalpy(Real n, Real e, Real T)
//  And it must also have the following member variables:
//    const int NSpecies

#include <ps_types.hpp>

template <typename EOSPolicy>
class EOS : private EOSPolicy {
  private:
    // Member functions
    using EOSPolicy::Temperature;
    using EOSPolicy::Energy;
    using EOSPolicy::Pressure;
    using EOSPolicy::Entropy;
    using EOSPolicy::Enthalpy;

    // Member variables
    // The number of particle species used by the EOS.
    using EOSPolicy::NSpecies;
    
  public:
    //! \fn Real GetTemperature(Real n, Real e, Real *Y)
    //  \brief Calculate the temperature from number density, energy density, and
    //         particle fractions.
    //
    //  \param[in] n  The number density
    //  \param[in] e  The energy density
    //  \param[in] Y  An array of particle fractions, expected to be of size NSpecies.
    //  \returns The temperature according to the EOS.
    inline Real GetTemperature(Real n, Real e, Real *Y){
      return Temperature(n, p, Y);
    }

    //! \fn Real GetEnergy(Real n, Real T, Real *Y)
    //  \brief Get the energy density from the numbe density, temperature, and
    //         particle fractions.
    //
    //  \param[in] n  The number density
    //  \param[in] T  The temperature
    //  \param[in] Y  An array of size NSpecies of the particle fractions.
    //  \returns The energy density according to the EOS.



    //! \fn int GetNSpecies() const
    //  \brief Get the number of particle species in this EOS.
    inline int GetNSpecies() const {
      return NSpecies;
    }
}

#endif
