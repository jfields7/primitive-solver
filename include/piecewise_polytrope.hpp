#ifndef PIECEWISE_POLYTROPE_H
#define PIECEWISE_POLYTROPE_H

//! \file piecewise_polytrope.hpp
//  \brief Defines a piecewise-polytropic equation of state.
//
//  Each individual piece satisfies the form
//  \f$P = \Kappa_i \rho^{\gamma_i}\f$,
//  for some density \f$\rho > \rho_i\f$. The temperature is
//  defined via the ideal gas law:
//  \f$P = nk_B T\f$.

#include <ps_types.hpp>
#include <eos_policy_interface.hpp>

namespace Primitive {

class PiecewisePolytrope : EOSPolicyInterface {
  private:
    /// Number of polytropes in the EOS
    int n_pieces;
    
    /// Parameters for the EOS
    Real *density_pieces;
    Real *a_pieces;
    Real *gamma_pieces;
    bool initialized;

    /// Allocate memory for the different EOS pieces.
    void AllocateMemory();

    /// Find the index of the piece that the density aligns with.
    int FindPiece(Real n) const;
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

    /// Calculate the energy per baryon (NOT energy per mass!)
    Real SpecificEnergy(Real n, Real T, Real *Y);

    /// Calculate the maximum allowed energy at a given density.
    Real MaximumEnergyAtDensity(Real n);

    /// Calculate the minimum allowed energy at a given density.
    Real MinimumEnergyAtDensity(Real n);

  public:
    /// Load the EOS parameters from a file.
    bool ReadParametersFromFile(std::string fname);

    /// Initialize PiecewisePolytrope from data.
    bool InitializeFromData(Real *densities, Real *kappas, Real *gammas, Real m, int n);

    /// Check if the EOS has been initialized properly.
    inline bool IsInitialized() const {
      return initialized;
    }
};

}; // namespace

#endif
