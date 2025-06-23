//! \test_units.cpp
//  \brief Unit tests for the various unit systems.

#include <iostream>
#include <cmath>

#include <unit_system.hpp>

#include <testing.hpp>
#include <testfunctions.hpp>

using namespace Primitive;

#define SQR(x) ((x)*(x))

#define CUBE(x) ((x)*(x)*(x))

bool ConversionCorrect(Real expected, Real actual, std::string name) {
  Real error = GetError(expected, actual);
  if (error > 1e-15) {
    std::cout << "Failed unit conversion for " << name << "\n";
    PrintError(expected, actual);
    return false;
  }

  return true;
}

bool TestGeometricKilometerUnits() {
  UnitSystem& geom = GeometricKilometer;

  // In geometric units, we should expect the following:
  // c = 1, G = 1, kb = 1, 1 cm = 1e-5.
  Real c = CGS.c * CGS.VelocityConversion(geom);
  bool result = ConversionCorrect(1.0, c, "c");

  Real G = CGS.G * CGS.LengthConversion(geom)*SQR(CGS.VelocityConversion(geom))/CGS.MassConversion(geom);
  result = result && ConversionCorrect(1.0, G, "G");

  Real kb = CGS.kb * CGS.EnergyConversion(geom)/CGS.TemperatureConversion(geom);
  result = result && ConversionCorrect(1.0, kb, "kb");

  Real cm = CGS.length * CGS.LengthConversion(geom);
  result = result && ConversionCorrect(1e-5, cm, "cm");

  return result;
}

bool TestGeometricSolarUnits() {
  UnitSystem& geom = GeometricSolar;

  // In geometric units, we should expect the following:
  // c = 1, G = 1, 1 Msol = 1
  Real c = CGS.c * CGS.VelocityConversion(geom);
  bool result = ConversionCorrect(1.0, c, "c");

  Real G = CGS.G * CGS.LengthConversion(geom)*SQR(CGS.VelocityConversion(geom))/CGS.MassConversion(geom);
  result = result && ConversionCorrect(1.0, G, "G");

  //Real kb = CGS.kb * CGS.EnergyConversion(geom)/CGS.TemperatureConversion(geom);
  //result = result && ConversionCorrect(1.0, kb, "kb");

  Real Msun = CGS.Msun * CGS.MassConversion(geom);
  result = result && ConversionCorrect(1.0, Msun, "Msun");

  return result;
}

bool TestNuclearUnits() {
  // In nuclear units, we should expect the following:
  // c = 1, kb = 1, MeV = 1, 1 cm = 1e13
  Real c = CGS.c * CGS.VelocityConversion(Nuclear);
  bool result = ConversionCorrect(1.0, c, "c");

  Real kb = CGS.kb * CGS.EnergyConversion(Nuclear)/CGS.TemperatureConversion(Nuclear);
  result = result && ConversionCorrect(1.0, kb, "kb");

  Real MeV = CGS.MeV * CGS.EnergyConversion(Nuclear);
  result = result && ConversionCorrect(1.0, MeV, "MeV");

  Real cm = CGS.length * CGS.LengthConversion(Nuclear);
  result = result && ConversionCorrect(1.0e13, cm, "cm");

  return result;
}

// Check for self-consistency of a unit system
bool TestConsistency(UnitSystem* units) {
  
  // Density
  bool result = ConversionCorrect(units->density, 1.0/CUBE(units->length), "density");

  // Pressure
  result = result && ConversionCorrect(units->pressure, units->energy*units->density, "pressure");

  return result;
}

int main(int argc, char* argv[]) {
  (void)argv[argc-1];

  UnitTests tester{"Unit Systems"};

  tester.RunTest(&TestGeometricKilometerUnits, "Geometric (Kilometer) Units -- Geometrization");
  tester.RunTest(&TestConsistency, "Geometric (Kilometer) Units -- Consistency", &GeometricKilometer);

  tester.RunTest(&TestGeometricSolarUnits, "Geometric (Solar) Units -- Geometrization");
  tester.RunTest(&TestConsistency, "Geometric (Solar) Units -- Consistency", &GeometricSolar);

  tester.RunTest(&TestNuclearUnits, "Nuclear Units -- Geometrization");
  tester.RunTest(&TestConsistency, "Nuclear Units -- Consistency", &Nuclear);

  tester.PrintSummary();

  return 0;
}
