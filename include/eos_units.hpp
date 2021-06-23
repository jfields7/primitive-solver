#ifndef EOS_UNITS_HPP
#define EOS_UNITS_HPP

/*!
 * \file eos_units.hpp
 * \brief contains unit definitions and conversion for the EOS solver.
 */

namespace EOSUnits{
  // Constants

  // Speed of light
  const double c_in_cgs = 2.99792458e10; // cgs units -- cm/s
  const double c_in_si  = 2.99792458e8;  // SI units -- m/s
  const double c_in_nat = 1.0; // Natural units
  const double c_in_code = 1.0; // Code units

  // Gravitational constant
  const double G_in_cgs = 6.67430e-8;  // cgs units -- cm^3 g^-1 s^-2
  const double G_in_si  = 6.67430e-11; // SI units -- m^3 kg^-1 s^-2
  const double G_in_nat = 1.0; // Natural units
  const double G_in_code = 1.0; // Code units

  // Boltzmann constant
  const double kb_in_cgs = 1.380649e-16; // cgs units -- erg K^-1
  const double kb_in_si  = 1.380649e-23; // SI units -- J K^-1
  const double kb_in_mev = 8.617333262145e1; // MeV/K
  const double kb_in_nat = 1.0; // natural units
  const double kb_in_code = 1.0; // Code units


  //Units

  // Length unit (fm)
  const double length_in_cgs = 1e-13; // cgs units -- cm
  const double length_in_si  = 1e-15; // SI units -- m
  const double length_in_nuc = 1.0;   // nuclear units
  const double length_in_code = 1.0; // Code units
  
  // Number density unit (fm^-3)
  const double nb_in_cgs = 1e39; // cgs units -- cm^-3
  const double nb_in_si  = 1e45; // SI units -- m^-3
  const double nb_in_nuc = 1.0;  // nuclear units
  const double nb_in_code = 1.0; // Code units

  // Energy units (MeV)
  const double energy_in_cgs  = 1.602176634e-6;  // cgs units -- erg
  const double energy_in_si   = 1.602176634e-13; // SI units -- J
  const double energy_in_natm = 1.323833314e-57; // Natural units -- m
  const double energy_in_nuc  = 1.0; // Nuclear units
  const double energy_in_code = 1.0; // Code units

  // Baryon mass (amu)
  /*const double baryon_in_cgs = 1.66053906660e-24; // cgs units -- g
  const double baryon_in_si  = 1.66053906660e-27; // SI units -- kg
  const double baryon_in_mev = 931.49410242; // MeV/c^2
  const double baryon_in_nuc = 1.0;          // nuclear units*/
  // Baryon mass, if it has to be defined, should be defined by the EOS.

  // Temperature (MeV)

  // Chemical potential (MeV)
};

#endif
