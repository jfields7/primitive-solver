//! \file benchmark.cpp
//! \brief Implementation file for BenchmarkMagW class

#include <benchmark.hpp>
#include <benchmark_magW.hpp>
#include <stdexcept>
#include <cmath>
#include <cstdio>

BenchmarkMagW::BenchmarkMagW(DataRange& range_ibeta, DataRange& range_W,
                     Real t_n, Real t_T, Real t_Ye,
                     std::string benchmark_name, bool save=true) {
  // Initialize the data ranges
  ibeta = InitializeFromDataRange(range_ibeta);
  W = InitializeFromDataRange(range_W);
  for (unsigned int i = 0; i < range_W.size; i++) {
    W[i] += 1.0;
  }

  // Initialize the sizes of the input fields
  nbeta = range_ibeta.size;
  nW = range_W.size;

  n = t_n;
  T = t_T;
  Ye = t_Ye;

  name = benchmark_name;

  // Initialize the output fields.
  unsigned int size;
  if (save) {
    size = nbeta*nW;
  } else {
    size = 1;
  }
  iterations = new int[size];
  error      = new Real[size];
  success    = new bool[size];
}

BenchmarkMagW::~BenchmarkMagW() {
  // Clear all the memory we just allocated
  delete[] iterations;
  delete[] error;
  delete[] success;

  delete[] ibeta;
  delete[] W;
}

Real* BenchmarkMagW::InitializeFromDataRange(DataRange& range) {
  if (range.size == 0) {
    throw std::domain_error("Size must be greater than zero.");
  }

  Real *data = new Real[range.size];
  
  // Calculate the increment between points.
  Real dx = 0.0;
  if (range.size > 1) {
    dx = (range.max - range.min)/(range.size - 1);
  }

  for (unsigned int i = 0; i < range.size; i++) {
    data[i] = range.min + i*dx;
  }

  if (range.log) {
    for (unsigned int i = 0; i < range.size; i++) {
      data[i] = std::pow(10.0, data[i]);
    }
  }

  return data;
}

void BenchmarkMagW::SaveBenchmark() {
  // Open the file
  std::FILE* file = std::fopen(name.c_str(), "w");

  // The columns for the .dat file are as follows:
  // # 1 n : 2 T : 3 iterations : 4 error : 5 success
  std::fprintf(file, "# 1 W : 2 ibeta : 3 iterations : 4 error : 5 success\n");

  // Write out each column
  for (unsigned int im = 0; im < nbeta; im++) {
    for (unsigned int iW = 0; iW < nW; iW++) {
      unsigned int idx = iW + nW*im;

      unsigned int succ = success[idx] ? 1 : 0;

      std::fprintf(file, "%20.15g  %20.15g  %u  %20.15g  %u\n",
                   W[iW], ibeta[im], iterations[idx], error[idx], succ);
    }
  }

  std::fclose(file);
}
