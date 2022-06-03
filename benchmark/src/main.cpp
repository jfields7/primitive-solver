//! \main.cpp
//  \brief Run the benchmarks and save the output data.

#include <benchmark.hpp>
#include <command_parser.hpp>
#include <paramreader.hpp>

#include <primitive_solver.hpp>
#include <eos.hpp>
#include <idealgas.hpp>
#include <do_nothing.hpp>

using namespace Primitive;

void ReadDataRange(ParamReader& params, DataRange& range, std::string section) {
  range.size = (unsigned int) params.readAsInt(section, "size");
  range.min = (Real) params.readAsDouble(section, "min");
  range.max = (Real) params.readAsDouble(section, "max");
}

int main(int argc, char *argv[]) {
  // Parse command line arguments.
  CommandParser parser;
  parser.AddString("input", true, "", true);
  parser.AddBoolean("save");

  // Tests
  CommandParser::Error error = parser.ParseArguments(argc, argv);

  // If there's an error, print out what it is.
  if (error.code != CommandParser::ErrorCode::SUCCESS) {
    std::cout << "Encountered the following error while reading parameters:\n";
    std::cout << "  Argument position: " << error.position << "\n";
    std::cout << "  Argument name: " << error.arg << "\n";
    std::cout << "  Error: " << CommandParser::GetErrorCodeName(error.code) << "\n";
    return 0;
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

  // Populate the data ranges using the parameter file.
  std::cout << "Populating data ranges...\n";
  DataRange n_data, T_data, ux_data, uy_data, uz_data, Bx_data, By_data, Bz_data;
  ReadDataRange(params, n_data, "Density");
  ReadDataRange(params, T_data, "Temperature");
  ReadDataRange(params, ux_data, "VelX");
  ReadDataRange(params, uy_data, "VelY");
  ReadDataRange(params, uz_data, "VelZ");
  ReadDataRange(params, Bx_data, "MagX");
  ReadDataRange(params, By_data, "MagY");
  ReadDataRange(params, Bz_data, "MagZ");

  // Set up the primitive solver.
  EOS<IdealGas, DoNothing> eos;
  PrimitiveSolver<IdealGas, DoNothing> ps{&eos};

  // Now we can construct and run the benchmark.
  bool save = parser.GetBoolean("save");
  Benchmark benchmark{n_data, T_data, ux_data, uy_data, uz_data, 
                      Bx_data, By_data, Bz_data, std::string("Basic"), save};
  std::cout << "Running benchmark...\n";
  if (!save) {
    benchmark.RunBenchmark(&ps, false);
  }
  else {
    benchmark.RunBenchmark(&ps, true);
    std::cout << "Saving benchmark...\n";
    benchmark.SaveBenchmark();
    std::cout << "Benchmark saved!\n";
  }

  return 0;
}
