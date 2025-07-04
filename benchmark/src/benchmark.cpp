//! \file benchmark.cpp
//! \brief Implementation file for Benchmark class

#include <benchmark.hpp>
#include <stdexcept>
#include <cmath>
#include <cstdio>

Benchmark::Benchmark(DataRange& range_n, DataRange& range_T,
                     Real t_ibeta, Real t_W, Real t_Ye,
                     std::string benchmark_name, bool save=true) {
  // Initialize the data ranges
  n = InitializeFromDataRange(range_n);
  T = InitializeFromDataRange(range_T);

  // Initialize the sizes of the input fields
  nn = range_n.size;
  nT = range_T.size;

  ibeta = t_ibeta;
  W = t_W;
  Ye = t_Ye;

  name = benchmark_name;

  // Initialize the output fields.
  unsigned int size;
  if (save) {
    size = nn*nT;
  } else {
    size = 1;
  }
  iterations = new int[size];
  error      = new Real[size];
  success    = new bool[size];
}

Benchmark::~Benchmark() {
  // Clear all the memory we just allocated
  delete[] iterations;
  delete[] error;
  delete[] success;

  delete[] n;
  delete[] T;
}

Real* Benchmark::InitializeFromDataRange(DataRange& range) {
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

void Benchmark::SaveBenchmark() {
  // Open the file
  std::FILE* file = std::fopen(name.c_str(), "w");

  // The columns for the .dat file are as follows:
  // # 1 n : 2 T : 3 iterations : 4 error : 5 success
  std::fprintf(file, "# 1 n : 2 T : 3 iterations : 4 error : 5 success\n");

  // Write out each column
  for (unsigned int it = 0; it < nT; it++) {
    for (unsigned int in = 0; in < nn; in++) {
      unsigned int idx = in + nn*it;

      unsigned int succ = success[idx] ? 1 : 0;

      std::fprintf(file, "%20.15g  %20.15g  %u  %20.15g  %u\n",
                   n[in], T[it], iterations[idx], error[idx], succ);
    }
  }

  std::fclose(file);
}
