#ifndef PIECEWISE_POLYTROPE_H
#define PIECEWISE_POLYTROPE_H

//! \file piecewise_polytrope.hpp
//  \brief Defines a piecewise-polytropic equation of state.

#include <ps_types.hpp>

namespace Primitive {

class PiecewisePolytrope {
  private:
    /// Number of polytropes in the EOS
    int n_pieces;
    
    /// Parameters for the EOS
    Real *density_pieces;
    Real *kappa_pieces;
    Real *a_pieces;
    Real *gamma_pieces;
    bool initialized;

    /// Allocate memory for the different EOS pieces.
    void AllocateMemory();

    /// Find the index of the piece that the density aligns with.
    int FindPiece(Real n) const;
  protected:
    /// An ideal gas doesn't really care about the particle species.
    const int n_species = 0;

    /// Baryon mass
    Real mb;

    /// Maximum density
    Real max_rho;
    /// Minimum density
    Real min_rho;

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

  public:
    /// Load the EOS parameters from a file.
    bool ReadParametersFromFile(std::string fname);

    /// Check if the EOS has been initialized properly.
    inline bool IsInitialized() const {
      return initialized;
    }
};

}; // namespace

#endif
