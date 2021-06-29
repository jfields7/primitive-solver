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
//    const int n_species
//    Real mb
//
//  It must also take an object for an error policy that implements the
//  following functions:
//    bool PrimitiveFloor(Real& n, Real& v[3], Real& p)
//    bool ConservedFloor(Real& D, Real& Sd[3], Real& tau, Real& Bu[3])
//  And the following protected variables:
//    Real n_atm
//    Real T_atm

#include <ps_types.hpp>

namespace Primitive {

template <typename EOSPolicy, typename ErrorPolicy>
class EOS : public EOSPolicy, public ErrorPolicy {
  private:
    // EOSPolicy member functions
    using EOSPolicy::Temperature;
    using EOSPolicy::TemperatureFromP;
    using EOSPolicy::Energy;
    using EOSPolicy::Pressure;
    using EOSPolicy::Entropy;
    using EOSPolicy::Enthalpy;
    using EOSPolicy::SoundSpeed;
    using EOSPolicy::SpecificEnergy;

    // EOSPolicy member variables
    // The number of particle species used by the EOS.
    using EOSPolicy::n_species;
    // The baryon mass
    using EOSPolicy::mb;

    // ErrorPolicy member functions
    using ErrorPolicy::PrimitiveFloor;
    using ErrorPolicy::ConservedFloor;

    // ErrorPolicy member variables
    using ErrorPolicy::n_atm;
    using ErrorPolicy::T_atm;

  public:
    //! \fn EOS()
    //  \brief Constructor for the EOS. It sets a default value for the floor.
    //
    //  n_atm gets fixed to 1e-10, and T_atm is set to 1.0.
    EOS() {
      n_atm = 1e-10;
      T_atm = 1.0;
    }

    //! \fn Real GetTemperature(Real n, Real e, Real *Y)
    //  \brief Calculate the temperature from number density, energy density, and
    //         particle fractions.
    //
    //  \param[in] n  The number density
    //  \param[in] e  The energy density
    //  \param[in] Y  An array of particle fractions, expected to be of size n_species.
    //  \return The temperature according to the EOS.
    inline Real GetTemperature(Real n, Real e, Real *Y) {
      return Temperature(n, e, Y);
    }

    //! \fn Real GetTemperatureFromP(Real n, Real p, Real *Y)
    //  \brief Calculate the temperature from number density, pressure, and
    //         particle fractions.
    //  \param[in] n  The number density
    //  \param[in] p  The pressure
    //  \param[in] Y  An array of particle fractions, expected to be of size n_species.
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
    //  \param[in] Y  An array of size n_species of the particle fractions.
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
    //  \param[in] Y  An array of size n_species of the particle fractions.
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
    //  \param[in] Y  An array of size n_species of the particle fractions.
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
    //  \param[in] Y  An array of size n_species of the particle fractions.
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
    //  \param[in] Y  An array of size n_species of the particle fractions.
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
    //  \param[in] Y  An array of size n_species of the particle fractions.
    //  \return The energy per baryon for the EOS.
    inline Real GetSpecificEnergy(Real n, Real T, Real *Y) {
      return SpecificEnergy(n, T, Y);
    }

    //! \fn int Getn_species() const
    //  \brief Get the number of particle species in this EOS.
    inline int GetNSpecies() const {
      return n_species;
    }

    //! \fn Real GetBaryonMass() const
    //  \brief Get the baryon mass used by this EOS.
    inline Real GetBaryonMass() const {
      return mb;
    }

    //! \fn void ApplyPrimitiveFloor(Real& n, Real& vu[3], Real& p, Real& T)
    //  \brief Apply the floor to the primitive variables.
    //
    //  \param[in,out] n  The number density
    //  \param[in,out] vu The velocity vector (contravariant)
    //  \param[in,out] p  The pressure
    //  \param[out]    T  The temperature
    //  \param[in]     Y  An array of size n_species of the particle fractions.
    //
    //  \return true if the primitives were adjusted, false otherwise.
    bool ApplyPrimitiveFloor(Real& n, Real v[3], Real& p, Real& T, Real *Y) {
      bool result = PrimitiveFloor(n, v, T);
      if (result) {
        p = Pressure(n, T, Y);
      }
      return result;
    }

    //! \fn void ApplyConservedFloor(Real& D, Real& Sd[3], Real& tau, Real& Bu[3])
    //  \brief Apply the floor to the conserved variables.
    //
    //  \param[in,out] D   The relativistic number density
    //  \param[in,out] Sd  The momentum density vector (covariant)
    //  \param[in,out] tau The tau variable (relativistic energy - D)
    //  \param[in,out] Bu  The magnetic field vector (contravariant)
    //  \param[in]     Y   An array of size_species of the particle fractions.
    //
    //  \return true if the conserved variables were adjusted, false otherwise.
    bool ApplyConservedFloor(Real& D, Real Sd[3], Real& tau, Real *Y) {
      return ConservedFloor(D, Sd, tau, GetTauFloor(Y));
    }

    //! \fn Real GetDensityFloor() const
    //  \brief Get the density floor used by the EOS ErrorPolicy.
    inline Real GetDensityFloor() const {
      return n_atm;
    }

    //! \fn Real GetPressureFloor(Real *Y) const
    //  \brief Get the pressure floor based on the current particle
    //         composition.
    //
    //  \param[in] Y A n_species-sized array of particle fractions.
    inline Real GetPressureFloor(Real *Y) {
      return GetPressure(n_atm, T_atm, Y);
    }

    //! \fn Real GetTemperatureFloor() const
    //  \brief Get the temperature floor used by the EOS ErrorPolicy.
    inline Real GetTemperatureFloor() const {
      return T_atm;
    }

    //! \fn Real GetTauFloor() const
    //  \brief Get the tau floor used by the EOS ErrorPolicy based
    //         on the current particle composition.
    //
    //  \param[in] Y A n_species-sized array of particle fractions.
    inline Real GetTauFloor(Real *Y) {
      return GetEnergy(n_atm, T_atm, Y) - mb*n_atm;
    }

    //! \fn Real SetDensityFloor(Real floor)
    //  \brief Set the density floor used by the EOS ErrorPolicy.
    //         Also adjusts the pressure and tau floor to be consistent.
    inline void SetDensityFloor(Real floor) {
      n_atm = (floor >= 0.0) ? floor : 0.0;
    }

    //! \fn Real SetPressureFloor(Real floor)
    //  \brief Set the pressure floor used by the EOS ErrorPolicy.
    //         Also adjusts the tau floor to be consistent.
    inline void SetTemperatureFloor(Real floor) {
      T_atm = (floor >= 0.0) ? floor : 0.0;
    }
};

}; // namespace

#endif
