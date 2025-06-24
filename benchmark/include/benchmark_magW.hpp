#ifndef BENCHMARK_MAGW_HPP_
#define BENCHMARK_MAGW_HPP_

#include <iostream>
#include <random>

#include <ps_types.hpp>
#include <string>
#include <primitive_solver.hpp>
#include <eos.hpp>
#include <geom_math.hpp>

struct DataRange;

class BenchmarkMagW {
 private:
  // Size of data sets
  unsigned int nbeta;
  unsigned int nW;

  // Data storage
  Real n; // Density
  Real T; // Temperature
  int *iterations; // Iterations per run
  Real *error; // Averaged error in result
  bool *success; // Success for each test
  Real *ibeta; // Inverse plasma beta (b^2/2p)
  Real *W; // Lorentz factor
  Real Ye; // Electron fraction

  // Name of the benchmark for data output.
  std::string name;

  // Initialize a data array with the specified ranges
  Real* InitializeFromDataRange(DataRange& range);

  inline Real GetError(Real expected, Real actual) {
    return std::fabs((actual - expected)/expected);
  }
 public:
  BenchmarkMagW(DataRange& range_ibeta, DataRange& range_W,
            Real t_n, Real t_T, Real t_Ye,
            std::string benchmark_name, bool save);

  ~BenchmarkMagW();

  // Run the benchmark for a given primitive solver.
  template<typename EOSPolicy, typename ErrorPolicy, class Boolean>
  void RunBenchmark(Primitive::PrimitiveSolver<EOSPolicy, ErrorPolicy>* ps,
                    const Boolean save) {
    // Assume a flat metric for simplicity.
    Real gd[NSPMETRIC] = {1.0, 0.0, 0.0, 1.0, 0.0, 1.0};
    Real gu[NSPMETRIC] = {1.0, 0.0, 0.0, 1.0, 0.0, 1.0};

    // Set up the random number generator for the velocity
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-1.0, 1.0);

    // Traverse the grid values
    for (unsigned int im = 0; im < nbeta; im++) {
      for (unsigned int iW = 0; iW < nW; iW++) {
        unsigned int idx = iW + nW*im;

        // Extract the primitive variables for this particular solution.
        Real prim[NPRIM];
        Real primold[NPRIM];
        Real bu[NMAG];
        Real cons[NCONS];

        prim[IDN] = primold[IDN] = n;
        prim[ITM] = primold[ITM] = T;
        prim[IYF] = primold[IYF] = Ye;
        prim[IPR] = primold[IPR] = ps->GetEOS()->GetPressure(prim[IDN], prim[ITM],
                                                             &prim[IYF]);
        // The velocity normalization based on our Lorentz factor.
        Real Wv = std::sqrt(W[iW]*W[iW] - 1.0);
        // Compute a random direction, which we then normalize.
        Real dir[3] = {dis(gen), dis(gen), dis(gen)};
        Real norm = std::sqrt(dir[0]*dir[0] + dir[1]*dir[1] + dir[2]*dir[2]);
        dir[0] /= norm;
        dir[1] /= norm;
        dir[2] /= norm;

        // Now compute the velocity.
        prim[IVX] = primold[IVX] = dir[0]*Wv;
        prim[IVY] = primold[IVY] = dir[1]*Wv;
        prim[IVZ] = primold[IVZ] = dir[2]*Wv;

        // Align the magnetic field with the velocity and rescale based on what beta
        // is. Since v || B, it follows that B^2 = b^2, so this is simple.
        Real b = std::sqrt(2.0*prim[IPR]*ibeta[im]);
        bu[IB1] = dir[0]*b;
        bu[IB2] = dir[1]*b;
        bu[IB3] = dir[2]*b;

        ps->PrimToCon(prim, cons, bu, gd);
        Primitive::SolverResult result = ps->ConToPrim(prim, cons, bu, gd, gu);

        if (save) {
          Real nerr = GetError(primold[IDN], prim[IDN]);
          Real vxerr = GetError(primold[IVX]/W[iW], prim[IVX]/W[iW]);
          Real vyerr = GetError(primold[IVY]/W[iW], prim[IVY]/W[iW]);
          Real vzerr = GetError(primold[IVZ]/W[iW], prim[IVZ]/W[iW]);
          Real epserr = GetError(
            ps->GetEOS()->GetSpecificInternalEnergy(primold[IDN],primold[ITM],&prim[IYF]),
            ps->GetEOS()->GetSpecificInternalEnergy(prim[IDN],primold[ITM],&prim[IYF]));
          
          error[idx] = 0.2*(nerr + vxerr + vyerr + vzerr + epserr);
          iterations[idx] = result.iterations;
          success[idx] = (nerr <= 10.0*ps->tol) & (vxerr <= 10.0*ps->tol) &
                         (vyerr <= 10.0*ps->tol) & (vzerr <= 10.0*ps->tol) &
                         (epserr <= 10.0*ps->tol);
        }
      }
    }
  }

  // Ouptut the benchmark results to a .dat ASCII file
  void SaveBenchmark();
};

#endif
