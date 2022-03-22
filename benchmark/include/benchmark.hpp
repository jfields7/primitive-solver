#ifndef BENCHMARK_HPP
#define BENCHMARK_HPP
//! \file benchmark.cpp
//  \brief Declare a class that defines a benchmark taken at
//  a range of initial densities, temperatures, and velocities.

#include <iostream>

#include <ps_types.hpp>
#include <string>
#include <primitive_solver.hpp>
#include <eos.hpp>
#include <geom_math.hpp>

// A utility struct for storing data ranges.
struct DataRange {
  unsigned int size;
  Real min;
  Real max;
};

class Benchmark {
  private:
    // Size of data sets
    unsigned int nn;
    unsigned int nT;
    unsigned int nux;
    unsigned int nuy;
    unsigned int nuz;
    unsigned int nBx;
    unsigned int nBy;
    unsigned int nBz;

    // Data storage
    Real *n; // Density
    Real *T; // Temperature
    Real *ux; // Velocity in the x direction
    Real *uy; // Velocity in the y direction
    Real *uz; // Velocity in the z direction
    Real *Bx; // Magnetic field in the x direction
    Real *By; // Magnetic field in the y direction
    Real *Bz; // Magnetic field in the z direction
    Real *iterations; // Iterations per run
    Real *n_errors; // Error in n over all runs
    Real *T_errors; // Error in n over all runs
    Real *ux_errors; // Error in n over all runs
    Real *uy_errors; // Error in n over all runs
    Real *uz_errors; // Error in n over all runs
    Real *magnetization; // Magnetization per run
    Real *thermal; // Thermal variable per run

    // Metric variables
    Real gd[NMETRIC];
    Real gu[NMETRIC];

    // Name of the benchmark for data output.
    std::string name;

    // Initialize a data array from storage.
    Real* InitializeFromDataRange(DataRange& range);

    inline unsigned int GetIndex(unsigned int in, unsigned int it,
                                 unsigned int iux, unsigned int iuy,
                                 unsigned int iuz, unsigned int ibx,
                                 unsigned int iby, unsigned int ibz) {
      return in + nn*(it + nT*(iux + nux*(iuy + nuy*(iuz + 
             nuz*(ibx + nBx*(iby + nBy*ibz))))));
    }

    Real GetError(Real expected, Real actual);

    void WriteRange(std::ofstream& file, unsigned int size, Real* data);

  public:
    Benchmark(DataRange& range_n, DataRange& range_T,
              DataRange& range_ux, DataRange& range_uy, DataRange& range_uz, 
              DataRange& range_Bx, DataRange& range_By, DataRange& range_Bz,
              std::string benchmark_name);

    ~Benchmark();

    // Run the benchmark for a given primitive solver.
    template<typename EOSPolicy, typename ErrorPolicy>
    void RunBenchmark(Primitive::PrimitiveSolver<EOSPolicy, ErrorPolicy>* ps) {
      // Traverse the grid of values using the mother of all nested loops.
      const NumTools::Root& root = ps->GetRootSolver();
      for (unsigned int ibz = 0; ibz < nBz; ibz++) {
        for (unsigned int iby = 0; iby < nBy; iby++) {
          for (unsigned int ibx = 0; ibx < nBx; ibx++) {
            for (unsigned int iuz = 0; iuz < nuz; iuz++) {
              for (unsigned int iuy = 0; iuy < nuy; iuy++) {
                for (unsigned int iux = 0; iux < nux; iux++) {
                  for (unsigned int it = 0; it < nT; it++) {
                    for (unsigned int in = 0; in < nn; in++) {
                      unsigned int idx = GetIndex(in, it, iux, iuy, iuz, ibx, iby, ibz);

                      //std::cout << "Point: " << idx << "\n";

                      // Extract the primitive variables for this particular solution.
                      Real prim[NPRIM];
                      Real primold[NPRIM];
                      Real bu[NMAG];
                      Real cons[NCONS];

                      prim[IDN] = primold[IDN] = n[in];
                      prim[ITM] = primold[ITM] = T[it];
                      prim[IPR] = primold[IPR] = ps->GetEOS()->GetPressure(prim[IDN], prim[ITM], nullptr);
                      prim[IVX] = primold[IVX] = ux[iux];
                      prim[IVY] = primold[IVY] = uy[iuy];
                      prim[IVZ] = primold[IVZ] = uz[iuz];
                      const int n_species = ps->GetEOS()->GetNSpecies();
                      // Right now there are no EOSes that use particle fractions.
                      for (int s = 0; s < n_species; s++) {
                        prim[IYF + s] = primold[IYF + s] = 0.0;
                      }

                      bu[IB1] = Bx[ibx];
                      bu[IB2] = By[iby];
                      bu[IB3] = Bz[ibz];

                      ps->PrimToCon(prim, cons, bu, gd);
                      //Primitive::Error result = ps->ConToPrim(prim, cons, bu, gd, gu);
                      ps->ConToPrim(prim, cons, bu, gd, gu);

                      // Save the errors from this calculation.
                      n_errors[idx] = GetError(primold[IDN], prim[IDN]);
                      T_errors[idx] = GetError(primold[ITM], prim[ITM]);
                      ux_errors[idx] = GetError(primold[IVX], prim[IVX]);
                      uy_errors[idx] = GetError(primold[IVY], prim[IVY]);
                      uz_errors[idx] = GetError(primold[IVZ], prim[IVZ]);

                      // Write the auxiliary quantities
                      // Calculate the thermal variable.
                      thermal[idx] = primold[IPR]/(primold[IDN]*ps->GetEOS()->GetBaryonMass());
                      // To get the magnetization, we first need to find B^2.
                      // Extract the 3-metric from the full spacetime metric.
                      Real g3d[NSPMETRIC] = {gd[I11], gd[I12], gd[I13],
                                                      gd[I22], gd[I23],
                                                               gd[I33]};
                      // Square the magnetic field.
                      Real Bsq = Primitive::SquareVector(bu, g3d);
                      magnetization[idx] = Bsq/cons[IDN];


                      // Get the number of iterations used by the primitive solver.
                      iterations[idx] = (Real) root.last_count;
                    }
                  }
                }
              }
            }
          }
        }
      }
    }

    // Output the benchmark results to a .bmrk binary file,
    // along with the data ranges in a .csv text file.
    void SaveBenchmark();

    // Set the metric for the benchmark. Default is Minkowski space.
    void SetMetric(Real g_d[NMETRIC], Real g_u[NMETRIC]);
};

#endif
