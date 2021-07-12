#ifndef EOS_UNITS_HPP
#define EOS_UNITS_HPP

//! \file eos_units.hpp
//  \brief contains unit definitions and conversion for the EOS solver.
//  
//  Each unit system is defined as its own struct inside the EOSUnits namespace.
//  TODO: Check that these conversions are correct.
//  TODO: Should I use the constants in CompOSE instead?

namespace EOSUnits{

//! CGS Units
//
//  Units with uncertainty are defined using the 2018 CODATA values.
struct CGS {
  static Real const c    = 2.99792458e10; //! cm/s
  static Real const G    = 6.67430e-8; //! cm^3 g^-1 s^-2
  static Real const kb   = 1.380649e-16; //! erg K^-1
  static Real const Msun = 1.98847e33; //! g
  static Real const MeV  = 1.602176634e-6; //! erg

  static Real const length  = 1.0;
  static Real const density = 1.0;
  static Real const mass = 1.0;
  static Real const energy = 1.0;
  static Real const pressure = 1.0;
  static Real const temperature = 1.0;
};

//! Geometric units
struct Geometric {
  static Real const c    = 1.0;
  static Real const G    = 1.0;
  static Real const kb   = 1.0;
  static Real const Msun = 1.4766696910334395;  //! km
  static Real const MeV  = 1.32383331356638e-60; //! km

  // Geometric units for converting from their cgs equivalent.
  static Real const length      = 1.0e-5; //! km/cm
  static Real const density     = 1.0e15; //! cm^3/km^3
  static Real const mass        = 7.426160269118665e-34; //! erg*G/c^2*km/cm
  static Real const energy      = 8.26271763969804e-55; //! cm*c^4/G/erg
  static Real const pressure    = 8.26271763969804e-40; //! c^4/G/cm^2/erg
  static Real const temperature = 1.140791284653146e-70; //! cm*c^4/G*kb/erg
};

//! Nuclear units
//
//  TODO: Especially double-check these.
struct Nuclear {
  static Real const c    = 1.0;
  static Real const G    = 1.323833313566383e-42; //! fm
  static Real const kb   = 1.0;
  static Real const Msun = 1.115449865100704e60; //! MeV
  static Real const MeV  = 1.0;

  // Nuclear units for converting from their cgs equivalent.
  static Real const length      = 1e13; //! fm/cm
  static Real const density     = 1e-39; //! cm^3/fm^3
  static Real const mass        = 5.609588603804452e26; //! MeV/(erg/c^2)
  static Real const energy      = 6.241509074460763e5;  //! MeV/erg
  static Real const pressure    = 6.241509074460763e34; //! (MeV/fm^3)/(erg/cm^3)
  static Real const temperature = 8.61733326214518e-11; //! MeV/(erg/kb)
};

}; // namespace

#endif
