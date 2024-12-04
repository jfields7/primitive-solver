#ifndef EOS_COMPOSE_H
#define EOS_COMPOSE_H

//! \file eos_compose.hpp
//  \brief Defines EOSTable, which stores information from a tabulated
//         equation of state in CompOSE format.
//
//  Tables should be generated using
//  <a href="https://bitbucket.org/dradice/pycompose">PyCompOSE</a>

///  \warning This code assumes the table to be uniformly spaced in
///           log nb, log t, and yq

#include <cstddef>
#include <string>

#include <ps_types.hpp>
#include <eos_policy_interface.hpp>

namespace Primitive {

class EOSCompOSE : public EOSPolicyInterface {
  public:
    enum TableVariables {
      ECLOGP  = 0,  //! log (pressure / 1 MeV fm^-3)
      ECENT   = 1,  //! entropy per baryon [kb]
      ECMUB   = 2,  //! baryon chemical potential [MeV]
      ECMUQ   = 3,  //! charge chemical potential [MeV]
      ECMUL   = 4,  //! lepton chemical potential [MeV]
      ECLOGE  = 5,  //! log (total energy density / 1 MeV fm^-3)
      ECCS    = 6,  //! sound speed [c]
      ECNVARS = 7
    };

  protected:
    /// Constructor
    EOSCompOSE();

    /// Destructor
    ~EOSCompOSE();

    /// Temperature from energy density
    Real TemperatureFromE(Real n, Real e, Real *Y);

    /// Calculate the temperature using.
    Real TemperatureFromP(Real n, Real p, Real *Y);

    /// Calculate the energy density using.
    Real Energy(Real n, Real T, Real *Y);

    /// Calculate the pressure using.
    Real Pressure(Real n, Real T, Real *Y);

    /// Calculate the entropy per baryon using.
    Real Entropy(Real n, Real T, Real *Y);

    /// Calculate the enthalpy per baryon using.
    Real Enthalpy(Real n, Real T, Real *Y);

    /// Calculate the sound speed.
    Real SoundSpeed(Real n, Real T, Real *Y);

    /// Calculate the specific internal energy per unit mass
    Real SpecificInternalEnergy(Real n, Real T, Real *Y);

    /// Calculate the baryon chemical potential
    Real BaryonChemicalPotential(Real n, Real T, Real *Y);

    /// Calculate the charge chemical potential
    Real ChargeChemicalPotential(Real n, Real T, Real *Y);

    /// Calculate the electron-lepton chemical potential
    Real ElectronLeptonChemicalPotential(Real n, Real T, Real *Y);

    /// Calculate hot (neutrino trapped) beta equilibrium T_eq and Y_eq given n, e, and Yl
    int BetaEquilibriumTrapped(Real n, Real e, Real *Yl, Real &T_eq, Real *Y_eq, Real T_guess, Real *Y_guess);

    /// Get the minimum enthalpy per baryon.
    Real MinimumEnthalpy();

    /// Get the minimum pressure at a given density and composition
    Real MinimumPressure(Real n, Real *Y);

    /// Get the maximum pressure at a given density and composition
    Real MaximumPressure(Real n, Real *Y);

    /// Get the minimum energy at a given density and composition
    Real MinimumEnergy(Real n, Real *Y);

    /// Get the maximum energy at a given density and composition
    Real MaximumEnergy(Real n, Real *Y);

  public:
    /// Reads the table file.
    void ReadTableFromFile(std::string fname);

    /// Get the raw number density
    Real const * GetRawLogNumberDensity() const {
      return m_log_nb;
    }
    Real const * GetRawYq() const {
      return m_yq;
    }
    /// Get the raw number density
    Real const * GetRawLogTemperature() const {
      return m_log_t;
    }
    /// Get the raw table data
    Real const * GetRawTable() const {
      return m_table;
    }

    // Indexing used to access the data
    inline ptrdiff_t index(int iv, int in, int iy, int it) const {
      return it + m_nt*(iy + m_ny*(in + m_nn*iv));
    }

    /// Check if the EOS has been initialized properly.
    inline bool IsInitialized() const {
      return m_initialized;
    }

    /// Set the number of species. Throw an exception if
    /// the number of species is invalid.
    void SetNSpecies(int n);

  private:
    /// Low level function, not intended for outside use
    Real temperature_from_var(int vi, Real var, Real n, Real Yq) const;
    /// Low level evaluation function, not intended for outside use
    Real eval_at_nty(int vi, Real n, Real T, Real Yq) const;
    /// Low level evaluation function, not intended for outside use
    Real eval_at_lnty(int vi, Real ln, Real lT, Real Yq) const;

    /// Functions for neutrino equilibrium
    int trapped_equilibrium_2DNR(Real n, Real e, Real Yle, Real x0[2], Real x1[2]);
    void func_eq_weak(Real n, Real e, Real Yle, Real x1[2], Real y[2]);
    Real error_func_eq_weak(Real Yle, Real e, Real y[2]);
    int jacobi_eq_weak(Real n, Real e, Real Yle, Real x1[2], Real J[2][2]);
    int eta_e_gradient(Real n, Real T, Real *Y, Real eta, Real &detadt, Real &detadye, Real &dedt, Real &dedye);
    void inv_jacobi(Real det, Real J[2][2], Real invJ[2][2]);

    /// Evaluate interpolation weight for density
    void weight_idx_ln(Real *w0, Real *w1, int *in, Real log_n) const;
    /// Evaluate interpolation weight for composition
    void weight_idx_yq(Real *w0, Real *w1, int *iy, Real yq) const;
    /// Evaluate interpolation weight for temperature
    void weight_idx_lt(Real *w0, Real *w1, int *it, Real log_t) const;

  private:
    // Inverse of table spacing
    Real m_id_log_nb, m_id_log_t, m_id_yq;
    // Table size
    int m_nn, m_nt, m_ny;
    // Minimum enthalpy per baryon
    Real m_min_h;

    // Table storage, care should be made to store these data on the GPU later
    Real * m_log_nb;
    Real * m_log_t;
    Real * m_yq;
    Real * m_table;

    bool m_initialized;

  private:
  // Constants for neutrino calculations
#define WR_SQR(x) ((x)*(x))
#define WR_CUBE(x) ((x)*(x)*(x))
#define WR_QUAD(x) ((x)*(x)*(x)*(x))

    const Real hc_mevfm = 1.23984172e3;           // hc    [MeV fm] (not reduced)
    const Real pi       = 3.14159265358979323846; // pi    [-]
    const Real pi2      = WR_SQR(pi);             // pi**2 [-]

    const Real pref1 = 4.0/3.0*pi/WR_CUBE(hc_mevfm); // 4/3 *pi/(hc)**3 [1/MeV^3/fm^3]
    const Real pref2 = 4.0*pi/WR_CUBE(hc_mevfm);     // 4*pi/(hc)**3    [1/MeV^3 fm^3]

    const Real cnst1 = 7.0*WR_QUAD(pi)/20.0;  // 7*pi**4/20  [-]
    const Real cnst2 = 7.0*WR_QUAD(pi)/5.0;   // 7*pi**4/5   [-]
    const Real cnst3 = 7.0*WR_QUAD(pi)/15.0;  // 7*pi**4/15  [-]
    const Real cnst4 = 14.0*WR_QUAD(pi)/15.0; // 14*pi**4/15 [-]
    const Real cnst5 = 7.0*WR_QUAD(pi)/60.0;  // 7*pi**4/60  [-]
    const Real cnst6 = 7.0*WR_QUAD(pi)/30.0;  // 7*pi**4/30  [-]

#undef WR_SQR
#undef WR_CUBE
#undef WR_QUAD
};

} // namespace Primitive

#endif
