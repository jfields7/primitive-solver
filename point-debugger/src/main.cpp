// C++ headers
#include <iostream>
#include <string>
#include <sstream>
#include <type_traits>
#include <cmath>

// PrimitiveSolver headers
#include <primitive_solver.hpp>
#include <idealgas.hpp>
#include <reset_floor.hpp>
#include <do_nothing.hpp>
#include <piecewise_polytrope.hpp>
#include <eos_compose.hpp>
#include <eos.hpp>
#include <ps_types.hpp>
#include <ps_error.hpp>
#include <geom_math.hpp>
#include <unit_system.hpp>

// Utility libraries
#include <command_parser.hpp>
#include <paramreader.hpp>

// Function prototypes
bool TestPoint(ParamReader& params);
template<class EOSPolicy> bool RunWithEOS(ParamReader& params);
template<class EOSPolicy,class ErrorPolicy> bool RunWithEOSAndError(ParamReader& params);

template<class ErrorPolicy>
void LoadEOSOptions(Primitive::EOS<Primitive::IdealGas,ErrorPolicy>&, ParamReader&);
template<class ErrorPolicy>
void LoadEOSOptions(Primitive::EOS<Primitive::PiecewisePolytrope,ErrorPolicy>&,
               ParamReader&);
template<class ErrorPolicy>
void LoadEOSOptions(Primitive::EOS<Primitive::EOSCompOSE,ErrorPolicy>&, ParamReader&);

template<class EOSPolicy>
void LoadErrorOptions(Primitive::EOS<EOSPolicy,Primitive::ResetFloor>&, ParamReader&);
template<class EOSPolicy>
void LoadErrorOptions(Primitive::EOS<EOSPolicy,Primitive::DoNothing>&, ParamReader&);

std::string CodeToString(Primitive::Error error) {
  switch(error) {
    case Primitive::Error::SUCCESS:
      return "SUCCESS";
    case Primitive::Error::RHO_TOO_BIG:
      return "RHO_TOO_BIG";
    case Primitive::Error::RHO_TOO_SMALL:
      return "RHO_TOO_SMALL";
    case Primitive::Error::NANS_IN_CONS:
      return "NANS_IN_CONS";
    case Primitive::Error::MAG_TOO_BIG:
      return "MAG_TOO_BIG";
    case Primitive::Error::BRACKETING_FAILED:
      return "BRACKETING_FAILED";
    case Primitive::Error::NO_SOLUTION:
      return "NO_SOLUTION";
    case Primitive::Error::CONS_FLOOR:
      return "CONS_FLOOR";
    case Primitive::Error::PRIM_FLOOR:
      return "PRIM_FLOOR";
    case Primitive::Error::CONS_ADJUSTED:
      return "CONS_ADJUSTED";
    case Primitive::Error::SLOW_CONVERGENCE:
      return "SLOW_CONVERGENCE";
  }
  return "UNKNOWN";
}

bool TestPoint(ParamReader& params) {
  // Find the requested EOS policy, check that it's valid, and load the correct
  // test with the specific EOS policy.
  std::string eos = params.readAsString(std::string("EOS"), "eos_policy");
  if (eos.compare("IdealGas") == 0) {
    return RunWithEOS<Primitive::IdealGas>(params);
  } else if (eos.compare("PiecewisePolytrope") == 0) {
    return RunWithEOS<Primitive::PiecewisePolytrope>(params);
  } else if (eos.compare("EOSCompOSE") == 0) {
    return RunWithEOS<Primitive::EOSCompOSE>(params);
  } else if (eos.compare("NULL") == 0) {
    std::cout << "Error: expected parameter EOS/eos_policy.\n"
              << "  Permitted options are:\n"
              << "    IdealGas\n"
              << "    PiecewisePolytrope\n"
              << "    EOSCompOSE\n";
  } else {
    std::cout << "Error: Unknown EOS policy: " << eos << "\n"
              << "  Permitted options are:\n"
              << "    IdealGas\n"
              << "    PiecewisePolytrope\n"
              << "    EOSCompOSE\n";
  }

  return false;
}

template<class EOSPolicy>
bool RunWithEOS(ParamReader& params) {
  // Find the requested error policy, check that it's valid, and load the correct
  // test with specific EOS and error policies.
  std::string error = params.readAsString(std::string("Error"),"error_policy");
  if (error.compare("ResetFloor") == 0) {
    return RunWithEOSAndError<EOSPolicy,Primitive::ResetFloor>(params);
  } else if (error.compare("DoNothing") == 0) {
    return RunWithEOSAndError<EOSPolicy,Primitive::DoNothing>(params);
  } else if (error.compare("NULL") == 0) {
    std::cout << "Error: expected parameter Error/error_policy.\n";
    std::cout << "  Permitted options are:\n"
              << "    ResetFloor\n"
              << "    DoNothing\n";
  } else {
    std::cout << "Error: Unknown error policy: " << error << "\n";
    std::cout << "  Permitted options are:\n"
              << "    ResetFloor\n"
              << "    DoNothing\n";
  }

  return false;
}

template<class EOSPolicy, class ErrorPolicy>
bool RunWithEOSAndError(ParamReader& params) { 
  // Create the correct PrimitiveSolver object.
  Primitive::EOS<EOSPolicy, ErrorPolicy> eos;
  Primitive::PrimitiveSolver<EOSPolicy, ErrorPolicy> ps{&eos};

  // Load options for this combination of EOS and Error Policy.
  // (This would be *much* easier with C++17)
  LoadEOSOptions(eos, params);
  LoadErrorOptions(eos, params);
  std::string units = params.readAsString("EOS", "code_units");
  if (units.compare("CGS") == 0) {
    eos.SetCodeUnitSystem(&Primitive::CGS);
  } else if (units.compare("GeometricKilometer") == 0) {
    eos.SetCodeUnitSystem(&Primitive::GeometricKilometer);
  } else if (units.compare("GeometricSolar") == 0) {
    eos.SetCodeUnitSystem(&Primitive::GeometricSolar);
  } else if (units.compare("Nuclear") == 0) {
    eos.SetCodeUnitSystem(&Primitive::Nuclear);
  } else if (units.compare("MKS") == 0) {
    eos.SetCodeUnitSystem(&Primitive::MKS);
  } else if (units.compare("NULL") == 0) {
    eos.SetCodeUnitSystem(&Primitive::Nuclear);
  } else {
    std::cout << "Error: Unknown unit system: " << units << "\n";
    std::cout << "  Permitted options are:\n"
              << "    CGS\n"
              << "    GeometricKilometer\n"
              << "    GeometricSolar\n"
              << "    Nuclear\n"
              << "    MKS\n";
    return false;
  }
  Real dfloor = params.readAsDouble("Error", "dfloor");
  if (dfloor <= 0.) {
    dfloor = 1e-10;
  }
  eos.SetDensityFloor(dfloor/eos.GetBaryonMass());
  Real tfloor = params.readAsDouble("Error", "tfloor");
  if (tfloor <= 0.) {
    tfloor = 1e-10;
  }
  eos.SetTemperatureFloor(tfloor);
  Real dthreshold = params.readAsDouble("Error", "dthreshold");
  if (dthreshold <= 0.) {
    dthreshold = 1.0;
  }
  eos.SetThreshold(dthreshold);
  Real vmax = params.readAsDouble("Error", "vmax");
  if (vmax <= 0.) {
    vmax = 1. - 1e-15;
  }
  eos.SetMaxVelocity(vmax);
  Real bsqmax = params.readAsDouble("Error", "bsqmax");
  if (bsqmax <= 0.) {
    bsqmax = 1e6;
  }
  eos.SetMaximumMagnetization(bsqmax);
  Real fail_tol = params.readAsDouble("Error", "fail_tol");
  if (fail_tol <= 0.) {
    fail_tol = 1e-5;
  }
  eos.SetFailureTolerance(fail_tol);
  for (int i = 0; i < eos.GetNSpecies(); i++) {
    std::stringstream ss;
    ss << "y" << (i+1) << "floor";
    Real yfloor = params.readAsDouble("Error", ss.str());
    if (yfloor <= 0.) {
      yfloor = 0.;
    }
    eos.SetSpeciesAtmosphere(yfloor, i);
  }

  // Load PrimitiveSolver options.
  int iterations = params.readAsInt("PrimitiveSolver", "iterations");
  if (iterations <= 0) {
    iterations = 30;
  }
  ps.GetRootSolver().iterations = iterations;
  Real tol = params.readAsDouble("PrimitiveSolver", "tol");
  if (tol <= 0) {
    tol = 1e-15;
  }
  ps.tol = tol;

  // Load conserved variables
  Real cons[NCONS] = {0.0};
  Real cons_old[NCONS] = {0.0};
  Real bu[NMAG] = {0.0};
  Real g3d[NSPMETRIC] = {0.0};

  cons[IDN] = cons_old[IDN] = params.readAsDouble("State", "D");
  cons[IM1] = cons_old[IM1] = params.readAsDouble("State", "Sx");
  cons[IM2] = cons_old[IM2] = params.readAsDouble("State", "Sy");
  cons[IM3] = cons_old[IM3] = params.readAsDouble("State", "Sz");
  cons[IEN] = cons_old[IEN] = params.readAsDouble("State", "tau");
  for (int i = 0; i < eos.GetNSpecies(); i++) {
    std::stringstream ss;
    ss << "Dy" << (i+1);
    cons[IYD+i] = cons_old[IYD+i] = params.readAsDouble("State", ss.str());
  }

  bu[IB1] = params.readAsDouble("State", "Bx");
  bu[IB2] = params.readAsDouble("State", "By");
  bu[IB3] = params.readAsDouble("State", "Bz");

  g3d[S11] = params.readAsDouble("State", "gxx");
  g3d[S12] = params.readAsDouble("State", "gxy");
  g3d[S13] = params.readAsDouble("State", "gxz");
  g3d[S22] = params.readAsDouble("State", "gyy");
  g3d[S23] = params.readAsDouble("State", "gyz");
  g3d[S33] = params.readAsDouble("State", "gzz");

  Real detg = Primitive::GetDeterminant(g3d);

  Real g3u[NSPMETRIC];
  // Construct the inverse.
  Real idetg = 1.0/detg;
  g3u[S11] = (g3d[S22]*g3d[S33] - g3d[S23]*g3d[S23])*idetg;
  g3u[S12] = (g3d[S13]*g3d[S23] - g3d[S12]*g3d[S33])*idetg;
  g3u[S13] = (g3d[S12]*g3d[S23] - g3d[S13]*g3d[S22])*idetg;
  g3u[S22] = (g3d[S11]*g3d[S33] - g3d[S13]*g3d[S13])*idetg;
  g3u[S23] = (g3d[S12]*g3d[S13] - g3d[S11]*g3d[S23])*idetg;
  g3u[S33] = (g3d[S11]*g3d[S22] - g3d[S12]*g3d[S12])*idetg;

  Real prim[NPRIM];
  Primitive::SolverResult result = ps.ConToPrim(prim, cons, bu, g3d, g3u);

  // Print out the results of the solve.
  std::cout << "PrimitiveSolver results: \n"
            << "  error code: " << CodeToString(result.error) << "\n"
            << "  iterations: " << result.iterations << "\n"
            << "  cons floor applied: " << result.cons_floor << "\n"
            << "  prim floor applied: " << result.prim_floor << "\n"
            << "  cons adjusted: " << result.cons_adjusted << "\n\n";
  std::cout << "Input conserved: \n"
            << "  D   = " << cons_old[IDN] << "\n"
            << "  Sx  = " << cons_old[IM1] << "\n"
            << "  Sy  = " << cons_old[IM2] << "\n"
            << "  Sz  = " << cons_old[IM3] << "\n"
            << "  tau = " << cons_old[IEN] << "\n"
            << "  Bx  = " << bu[IB1] << "\n"
            << "  By  = " << bu[IB2] << "\n"
            << "  Bz  = " << bu[IB3] << "\n\n";
  std::cout << "Input metric: \n"
            << "  gxx = " << g3d[S11] << "\n"
            << "  gxy = " << g3d[S12] << "\n"
            << "  gxz = " << g3d[S13] << "\n"
            << "  gyy = " << g3d[S22] << "\n"
            << "  gyz = " << g3d[S23] << "\n"
            << "  gzz = " << g3d[S33] << "\n\n";
  if (result.error == Primitive::Error::SUCCESS) {
    std::cout << "Calculated primitives: \n"
              << "  n  = " << prim[IDN] << "\n"
              << "  vx = " << prim[IVX] << "\n"
              << "  vy = " << prim[IVY] << "\n"
              << "  vz = " << prim[IVZ] << "\n"
              << "  T  = " << prim[ITM] << "\n"
              << "  P  = " << prim[IPR] << "\n";
    for (int i = 0; i < eos.GetNSpecies(); i++) {
      std::cout << "  Y" << (i+1) << " = " << prim[IYF + i] << "\n";
    }
    std::cout << "\n";
  }

  // For a sanity check, convert the primitive variables back to conserved variables.
  Real cons_new[NCONS];
  ps.PrimToCon(prim, cons_new, bu, g3d);
  Real err[NCONS];
  int nhyd = NHYDRO - MAX_SPECIES + eos.GetNSpecies();
  for (int i = 0; i < nhyd; i++) {
    if (std::fabs(cons_old[i]) > 0) {
      err[i] = std::fabs((cons_new[i] - cons_old[i])/cons_old[i]);
    }
    else {
      err[i] = std::fabs(cons_new[i] - cons_old[i]);
    }
  }
  // Print out the errors
  std::cout << "Errors (relative unless it should be zero):\n"
            << "  D:   " << err[IDN] << "\n"
            << "  Sx:  " << err[IM1] << "\n"
            << "  Sy:  " << err[IM2] << "\n"
            << "  Sz:  " << err[IM3] << "\n"
            << "  tau: " << err[IEN] << "\n";
  for (int i = 0; i < eos.GetNSpecies(); i++) {
    std::cout << "  Dy" << i << ": " << err[IYD+i] << "\n";
  }
  std::cout << "\n";

  return true;
}

template<class ErrorPolicy>
void LoadEOSOptions(Primitive::EOS<Primitive::IdealGas,ErrorPolicy>& eos,
    ParamReader& params) {
  std::cout << "Loading parameters for IdealGas...\n";
  Real gamma = params.readAsDouble("EOS", "gamma");
  if (gamma <= 0) {
    gamma = 5./3.;
  }
  eos.SetGamma(gamma);
}

template<class ErrorPolicy>
void LoadEOSOptions(Primitive::EOS<Primitive::PiecewisePolytrope,ErrorPolicy>& eos,
    ParamReader& params) {
  std::cout << "Loading parameters for PiecewisePolytrope...\n";
  int n_pieces = params.readAsInt("EOS", "npieces");
  Real *densities = new Real[n_pieces+1];
  Real *gammas = new Real[n_pieces+1];

  Real poly_rmd = Primitive::MKS.MassDensityConversion(*(eos.GetEOSUnitSystem()))*
                  params.readAsDouble("EOS", "pwp_poly_rmd");
  for (int i = 0; i < n_pieces; i++) {
    std::stringstream ss;
    ss << "pwp_density_pieces_" << i;
    densities[i] = params.readAsDouble("EOS", ss.str())*
                   Primitive::MKS.MassDensityConversion(*(eos.GetEOSUnitSystem()));
  }
  for (int i = 0; i < n_pieces; i++) {
    std::stringstream ss;
    ss << "pwp_gamma_pieces_" << i;
    gammas[i] = params.readAsDouble("EOS", ss.str());
  }

  Real P0 = densities[1]*std::pow(densities[1]/poly_rmd, gammas[0] - 1.0);

  Real gamma_thermal = params.readAsDouble("EOS", "gamma_thermal");
  eos.InitializeFromData(densities, gammas, P0, 1.0, n_pieces);
  eos.SetThermalGamma(gamma_thermal);

  delete densities;
  delete gammas;
}

template<class ErrorPolicy>
void LoadEOSOptions(Primitive::EOS<Primitive::EOSCompOSE,ErrorPolicy>& eos,
    ParamReader& params) {
  std::cout << "Loading parameters for EOSCompOSE...\n";
  std::string fname = params.readAsString("EOS", "table");
  eos.ReadTableFromFile(fname);
}

template<class EOSPolicy>
void LoadErrorOptions(Primitive::EOS<EOSPolicy,Primitive::DoNothing>& eos,
    ParamReader& params) {
  return;
}

template<class EOSPolicy>
void LoadErrorOptions(Primitive::EOS<EOSPolicy,Primitive::ResetFloor>& eos,
    ParamReader& params) {
  return;
}

int main(int argc, char* argv[]) {
  // Parse command line arguments.
  CommandParser parser;
  parser.AddString("input", true, "", true);
  
  CommandParser::Error error = parser.ParseArguments(argc, argv);

  // If there's an error, print out what it is.
  if (error.code != CommandParser::ErrorCode::SUCCESS) {
    std::cout << "Encountered the following error while reading parameters:\n";
    std::cout << "  Argument position: " << error.position << "\n";
    std::cout << "  Argument name: " << error.arg << "\n";
    std::cout << "  Error: " << CommandParser::GetErrorCodeName(error.code) << "\n";
    return -1;
  }

  // Read in the parameter file.
  std::cout << "Reading parameters...\n";
  ParamReader params;
  ParamReader::ParamResult result;
  std::string filename = parser.GetString("input");
  result = params.readFile(filename);
  if (result != ParamReader::SUCCESS) {
    std::cout << "There was an error reading the file " << filename << ".\n";
    return -1;
  }

  bool success = TestPoint(params);

  if (!success) {
    return -1;
  }

  return 0;
}
