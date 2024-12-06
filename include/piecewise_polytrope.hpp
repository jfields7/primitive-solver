#ifndef PIECEWISE_POLYTROPE_H
#define PIECEWISE_POLYTROPE_H

//! \file piecewise_polytrope.hpp
//  \brief Defines a piecewise-polytropic equation of state.
//
//  Each individual piece satisfies the form
//  \f$P_\textrm{cold} = P_i \frac{\rho}{\rho_i}^{\gamma_i}\f$,
//  for some density \f$\rho > \rho_i\f$. There is an additional
//  finite-temperature portion added on top using the ideal gas
//  law:
//  \f$P_\textrm{therm} = nk_B T\f$

#include <cassert>

#include <ps_types.hpp>
#include <eos_policy_interface.hpp>
#include <unit_system.hpp>

namespace Primitive {

class PiecewisePolytrope : public EOSPolicyInterface {
  protected:
    /// Number of polytropes in the EOS
    int n_pieces;
    
    /// Parameters for the EOS
    Real *density_pieces;
    Real *gamma_pieces;
    Real *pressure_pieces;
    Real *eps_pieces;
    Real gamma_thermal;
    bool initialized;

    /// Allocate memory for the different EOS pieces.
    void AllocateMemory();
  protected:
    /// Constructor
    PiecewisePolytrope();

    /// Destructor
    ~PiecewisePolytrope();

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

    /// Calculate the internal energy per mass.
    Real SpecificInternalEnergy(Real n, Real T, Real *Y);

    /// Calculate the baryon chemical potential
    [[ noreturn ]]
    Real BaryonChemicalPotential(Real n, Real T, Real *Y);

    /// Calculate the charge chemical potential
    [[ noreturn ]]
    Real ChargeChemicalPotential(Real n, Real T, Real *Y);

    /// Calculate the electron-lepton chemical potential
    [[ noreturn ]]
    Real ElectronLeptonChemicalPotential(Real n, Real T, Real *Y);

    [[ noreturn ]]
    int BetaEquilibriumTrapped(Real n, Real e, Real *Yl, Real &T_eq, Real *Y_eq, Real T_guess, Real *Y_guess);

    [[ noreturn ]]
    void TrappedNeutrinos(Real n, Real T, Real *Y, Real n_nu[3], Real e_nu[3]);

    /// Calculate the minimum pressure at a given density and composition
    Real MinimumPressure(Real n, Real *Y);

    /// Calculate the maximum pressure at a given density and composition
    Real MaximumPressure(Real n, Real *Y);

    /// Calculate the minimum energy at a given density and composition
    Real MinimumEnergy(Real n, Real *Y);

    /// Calculate the maximum energy at a given density and composition
    Real MaximumEnergy(Real n, Real *Y);

  public:
    /// Load the EOS parameters from a file.
    bool ReadParametersFromFile(std::string fname);

    //! \brief Initialize PiecewisePolytrope from data.
    //  
    //  \param[in] densities The dividing densities
    //  \param[in] gammas    The adiabatic index for each polytrope
    //  \param[in] P0        The pressure at the first polytrope division
    //  \param[in] m         The baryon mass
    //  \param[in] n         The number of pieces in the EOS
    bool InitializeFromData(Real *densities, Real *gammas, 
                            Real P0, Real m, int n);

    /// Check if the EOS has been initialized properly.
    inline bool IsInitialized() const {
      return initialized;
    }

    /// Find out how many polytropes are in the EOS.
    inline int GetNPieces() const {
      return n_pieces;
    }

    /// Get the adiabatic constant for a particular density.
    inline Real GetGamma(Real n) const {
      return gamma_pieces[FindPiece(n)];
    }

    /// Set the adiabatic constant for the thermal part.
    inline void SetThermalGamma(Real g) {
      assert(g > 1.0);
      gamma_thermal = g;
    }

    /// Get the adiabatic constant for the thermal part.
    inline Real GetThermalGamma() const {
      return gamma_thermal;
    }

    /// Find the index of the piece that the density aligns with.
    int FindPiece(Real n) const;

    /// Polytropic Energy Density
    Real GetColdEnergy(Real n, int p);

    /// Polytropic Pressure
    Real GetColdPressure(Real n, int p);

    /// Set the number of species. Throw an exception if
    /// the number of species is invalid.
    void SetNSpecies(int n);

    /// Set the EOS unit system
    inline void SetEOSUnitSystem(UnitSystem* units) {
      eos_units = units;
    }
};

} // namespace

#endif
