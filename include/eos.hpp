#ifndef EOS_HPP
#define EOS_HPP

//! \file eos.hpp
//  \brief Defines an equation of state.
//
//  EOS is effectively an interface that describes how to create an
//  equation of state. It must be instantiated with an object implementing
//  the following protected functions:
//    Real Temperature(Real n, Real e, Real *Y)
//    Real TemperatureFromP(Real n, Real p, Real *Y)
//    Real Energy(Real n, Real T, Real *Y)
//    Real Pressure(Real n, Real T, Real *Y)
//    Real Entropy(Real n, Real T, Real *Y)
//    Real Enthalpy(Real n, Real T, Real *Y)
//    Real SoundSpeed(Real n, Real T, Real *Y)
//    Real SpecificEnergy(Real n, Real T, Real *Y)
//  And it must also have the following protected member variables:
//    const int NSpecies

#include <ps_types.hpp>

template <typename EOSPolicy>
class EOS : public EOSPolicy {
  private:
    // Member functions
    using EOSPolicy::Temperature;
    using EOSPolicy::TemperatureFromP;
    using EOSPolicy::Energy;
    using EOSPolicy::Pressure;
    using EOSPolicy::Entropy;
    using EOSPolicy::Enthalpy;
    using EOSPolicy::SoundSpeed;
    using EOSPolicy::SpecificEnergy;

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
    //  \return The temperature according to the EOS.
    inline Real GetTemperature(Real n, Real e, Real *Y) {
      return Temperature(n, e, Y);
    }

    //! \fn Real GetTemperatureFromP(Real n, Real p, Real *Y)
    //  \brief Calculate the temperature from number density, pressure, and
    //         particle fractions.
    //  \param[in] n  The number density
    //  \param[in] p  The pressure
    //  \param[in] Y  An array of particle fractions, expected to be of size NSpecies.
    //  \return The temperature according to the EOS.
    inline Real GetTemperatureFromP(Real n, Real p, Real *Y) {
      return TemperatureFromP(n, p, Y);
    }

    //! \fn Real GetEnergy(Real n, Real T, Real *Y)
    //  \brief Get the energy density from the number density, temperature, and
    //         particle fractions.
    //
    //  \param[in] n  The number density
    //  \param[in] T  The temperature
    //  \param[in] Y  An array of size NSpecies of the particle fractions.
    //  \return The energy density according to the EOS.
    inline Real GetEnergy(Real n, Real T, Real *Y) {
      return Energy(n, T, Y);
    }

    //! \fn Real GetPressure(Real n, Real T, Real *Y)
    //  \brief Get the pressure from the number density, temperature, and
    //         particle fractions.
    //
    //  \param[in] n  The number density
    //  \param[in] T  The temperature
    //  \param[in] Y  An array of size NSpecies of the particle fractions.
    //  \return The pressure according to the EOS.
    inline Real GetPressure(Real n, Real T, Real *Y) {
      return Pressure(n, T, Y);
    }

    //! \fn Real GetEntropy(Real n, Real T, Real *Y)
    //  \brief Get the entropy per baryon from the number density, temperature,
    //         and particle fractions.
    //
    //  \param[in] n  The number density
    //  \param[in] T  The temperature
    //  \param[in] Y  An array of size NSpecies of the particle fractions.
    //  \return The entropy per baryon for this EOS.
    inline Real GetEntropy(Real n, Real T, Real *Y) {
      return Entropy(n, T, Y);
    }

    //! \fn Real GetEnthalpy(Real n, Real T, Real *Y)
    //  \brief Get the enthalpy per baryon from the number density, temperature,
    //         and particle fractions.
    //
    //  \param[in] n  The number density
    //  \param[in] T  The temperature
    //  \param[in] Y  An array of size NSpecies of the particle fractions.
    //  \return The enthalpy per baryon for this EOS.
    inline Real GetEnthalpy(Real n, Real T, Real *Y) {
      return Enthalpy(n, T, Y);
    }

    //! \fn Real GetSoundSpeed(Real n, Real T, Real *Y)
    //  \brief Get the sound speed from the number density, temperature, and 
    //         particle fractions.
    //
    //  \param[in] n  The number density
    //  \param[in] T  The temperature
    //  \param[in] Y  An array of size NSpecies of the particle fractions.
    //  \return The sound speed for this EOS.
    inline Real GetSoundSpeed(Real n, Real T, Real *Y) {
      return SoundSpeed(n, T, Y);
    }

    //! \fn Real GetSpecificEnergy(Real n, Real T, Real *Y)
    //  \brief Get the energy per baryon from the number density, temperature,
    //         and particle fractions.
    //
    //  \param[in] n  The number density
    //  \param[in] T  The temperature
    //  \param[in] Y  An array of size NSpecies of the particle fractions.
    //  \return The energy per baryon for the EOS.
    inline Real GetSpecificEnergy(Real n, Real T, Real *Y) {
      return SpecificEnergy(n, T, Y);
    }

    //! \fn int GetNSpecies() const
    //  \brief Get the number of particle species in this EOS.
    inline int GetNSpecies() const {
      return NSpecies;
    }
};

#endif
