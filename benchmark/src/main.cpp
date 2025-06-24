//! \main.cpp
//! \brief Run the benchmarks and save the output data.
#include <cmath>

#include <benchmark.hpp>
#include <benchmark_magW.hpp>
#include <command_parser.hpp>
#include <paramreader.hpp>

#include <primitive_solver.hpp>
#include <eos.hpp>
#include <idealgas.hpp>
#include <piecewise_polytrope.hpp>
#include <eos_compose.hpp>
#include <do_nothing.hpp>
#include <reset_floor.hpp>

using namespace Primitive;

void ReadDataRange(ParamReader& params, DataRange& range, std::string section) {
  range.size = (unsigned int) params.readAsInt(section, "size");
  range.min = (Real) params.readAsDouble(section, "min");
  range.max = (Real) params.readAsDouble(section, "max");
  if (params.hasParameter(section, "log")) {
    range.log = (bool) params.readAsInt(section, "log");
    if (range.log) {
      range.min = std::log10(range.min);
      range.max = std::log10(range.max);
    }
  } else {
    range.log = false;
  }
}

int main(int argc, char *argv[]) {
  CommandParser parser;
  parser.AddString("input", true, "", true);
  parser.AddBoolean("save");

  // Tests
  CommandParser::Error error = parser.ParseArguments(argc, argv);

  // If there's an error, print out what it is.
  if (error.code != CommandParser::ErrorCode::SUCCESS) {
    std::cout << "Encountered the following error while reading parameters:\n"
              << "  Argument position: " << error.position << "\n"
              << "  Argument name: " << error.arg << "\n"
              << "  Error: " << CommandParser::GetErrorCodeName(error.code) << "\n";
  }

  // Read in the parameter file.
  std::cout << "Reading parameters...\n";
  ParamReader params;
  ParamReader::ParamResult result;
  std::string filename = parser.GetString("input");
  result = params.readFile(filename);
  if (result != ParamReader::SUCCESS) {
    std::cout << "There was an error reading the file " << filename << ".\n";
  }

  // Collect general information about the run
  std::string name = params.readAsString("Info", "name");
  std::string name_Wbeta = params.readAsString("Info", "name_Wbeta");
  Real tol = params.readAsDouble("Info", "tol");
  int max_iters = params.readAsInt("Info", "max_iters");
  Real ibeta = params.readAsDouble("Info", "ibeta");
  Real W = params.readAsDouble("Info", "W");
  Real Ye = params.readAsDouble("Info", "Ye");
  Real n = params.readAsDouble("Info", "n");
  Real T = params.readAsDouble("Info", "T");

  // Populate the data ranges using the parameter file.
  std::cout << "Populating data ranges...\n";
  DataRange n_range, T_range, W_range, ibeta_range;
  ReadDataRange(params, n_range, "Density");
  ReadDataRange(params, T_range, "Temperature");
  ReadDataRange(params, W_range, "Lorentz");
  ReadDataRange(params, ibeta_range, "IBeta");

  // Now we can construct and run the benchmark.
  bool save = parser.GetBoolean("save");
  Benchmark benchmark{n_range, T_range, ibeta, W, Ye, name, save};
  BenchmarkMagW benchmark_magW{ibeta_range, W_range, n, T, Ye, name_Wbeta, save};

  // EOS information; set up the primitive solver and run the benchmark at the same time.
  std::string eos = params.readAsString("EOS", "eos");
  if (eos == "ideal") {
    EOS<IdealGas, DoNothing> eos;
    PrimitiveSolver<IdealGas, DoNothing> ps{&eos};
    ps.tol = tol;
    ps.GetRootSolver().iterations = max_iters;

    benchmark.RunBenchmark(&ps, save);
    benchmark_magW.RunBenchmark(&ps, save);
  } else if (eos == "compose") {
    EOS<EOSCompOSE, ResetFloor> eos;
    eos.ReadTableFromFile(params.readAsString("EOS", "table"));
    eos.SetDensityFloor(eos.GetMinimumDensity());
    eos.SetTemperatureFloor(eos.GetMinimumTemperature());

    PrimitiveSolver<EOSCompOSE, ResetFloor> ps{&eos};
    ps.tol = tol;
    ps.GetRootSolver().iterations = max_iters;

    benchmark.RunBenchmark(&ps, save);
    benchmark_magW.RunBenchmark(&ps, save);
  }

  if (save) {
    benchmark.SaveBenchmark();
    benchmark_magW.SaveBenchmark();
  }

  return 0;
}
