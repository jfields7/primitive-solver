//! \file benchmark.cpp
//  \brief Implementation file for Benchmark class

#include <benchmark.hpp>
#include <stdexcept>
#include <cmath>
#include <fstream>
#include <sstream>

Benchmark::Benchmark(DataRange& range_n, DataRange& range_T,
                          DataRange& range_ux, DataRange& range_uy, DataRange& range_uz,
                          DataRange& range_Bx, DataRange& range_By, DataRange& range_Bz,
                          std::string benchmark_name, bool save=true) {
  // Initialize all the data ranges.
  n = InitializeFromDataRange(range_n);
  T = InitializeFromDataRange(range_T);
  ux = InitializeFromDataRange(range_ux);
  uy = InitializeFromDataRange(range_uy);
  uz = InitializeFromDataRange(range_uz);
  Bx = InitializeFromDataRange(range_Bx);
  By = InitializeFromDataRange(range_By);
  Bz = InitializeFromDataRange(range_Bz);

  // Initialize the sizes of the input fields.
  nn = range_n.size;
  nT = range_T.size;
  nux = range_ux.size;
  nuy = range_uy.size;
  nuz = range_uz.size;
  nBx = range_Bx.size;
  nBy = range_By.size;
  nBz = range_Bz.size;

  // Initialize the output fields.
  unsigned int size;
  if (save) {
    size = nn*nT*nux*nuy*nuz*nBx*nBy*nBz;
  }
  else {
    size = 1;
  }
  iterations = new Real[size];
  n_errors = new Real[size];
  T_errors = new Real[size];
  ux_errors = new Real[size];
  uy_errors = new Real[size];
  uz_errors = new Real[size];
  magnetization = new Real[size];
  thermal = new Real[size];

  // Initialize the metric to Minkowski space.
  gd[I00] = gu[I00] = -1.0;
  gd[I01] = gu[I01] = 0.0;
  gd[I02] = gu[I02] = 0.0;
  gd[I03] = gu[I03] = 0.0;
  gd[I11] = gu[I11] = 1.0;
  gd[I12] = gu[I12] = 0.0;
  gd[I13] = gu[I13] = 0.0;
  gd[I22] = gu[I22] = 1.0;
  gd[I23] = gu[I23] = 0.0;
  gd[I33] = gu[I33] = 1.0;

  name = benchmark_name;
}

Benchmark::~Benchmark() {
  // Clear all the memory that we just allocated.
  delete[] iterations;
  delete[] n_errors;
  delete[] T_errors;
  delete[] ux_errors;
  delete[] uy_errors;
  delete[] uz_errors;

  delete[] n;
  delete[] T;
  delete[] ux;
  delete[] uy;
  delete[] uz;
  delete[] Bx;
  delete[] By;
  delete[] Bz;
}

Real* Benchmark::InitializeFromDataRange(DataRange& range) {
  Real *data = new Real[range.size];

  if (range.size == 0) {
    throw std::domain_error("Size must be greater than zero.");
  }

  // Calculate the increment between points.
  Real dx = 0.0;
  if (range.size > 1) {
    dx = (range.max - range.min)/(range.size - 1);
  }

  for (unsigned int i = 0; i < range.size; i++) {
    data[i] = range.min + i*dx;
  }

  return data;
}

Real Benchmark::GetError(Real expected, Real actual) {
  return std::fabs((expected - actual)/expected);
}

void Benchmark::SetMetric(Real g_d[NMETRIC], Real g_u[NMETRIC]) {
  for (unsigned int i = 0; i < NMETRIC; i++) {
    gd[i] = g_d[i];
    gu[i] = g_u[i];
  }
}

void Benchmark::WriteRange(std::ofstream& file, unsigned int size, Real* data) {
  std::stringstream ss;
  ss << size << ", " << n[0] << ", " << n[size-1] << "\n";
  file << ss.str();
}

void Benchmark::SaveBenchmark() {
  // Write out the data ranges.
  std::stringstream filename;
  filename << name << ".csv";
  std::ofstream ranges(filename.str(), std::ios::out | std::ios::trunc);
  if (ranges.is_open()) {
    WriteRange(ranges, nn, n);
    WriteRange(ranges, nT, T);
    WriteRange(ranges, nux, ux);
    WriteRange(ranges, nuy, uy);
    WriteRange(ranges, nuz, uz);
    WriteRange(ranges, nBx, Bx);
    WriteRange(ranges, nBy, By);
    WriteRange(ranges, nBz, Bz);
    ranges.close();
  }
  else {
    throw std::ios_base::failure("There was an error writing the ranges.");
  }

  // The format for the binary file is as follows:
  // iterations[0], iterations[1], ..., iterations[size-1]
  // n_errors[0], n_errors[1], ..., n_errors[size-1]
  // ...
  // uz_errors[0], uz_errors[1], ..., uz_errors[size-1]
  unsigned int size = nn*nT*nux*nuy*nuz*nBx*nBy*nBz;
  filename.str(std::string());
  filename << name << ".bmrk";
  std::ofstream data(filename.str(), std::ios::out | std::ios::trunc | std::ios::binary);
  if (data.is_open()) {
    data.write((char *) iterations, size*sizeof(Real));
    data.write((char *) n_errors, size*sizeof(Real));
    data.write((char *) T_errors, size*sizeof(Real));
    data.write((char *) ux_errors, size*sizeof(Real));
    data.write((char *) uy_errors, size*sizeof(Real));
    data.write((char *) uz_errors, size*sizeof(Real));
    data.close();
  }
  else {
    throw std::ios_base::failure("There was an error writing the benchmark.");
  }

  // Write out the auxiliary file. It has the same format as
  // the benchmark file, but it contains the two auxiliary
  // quantities, magnetization and P/rho.
  filename.str(std::string());
  filename << name << ".aux";
  std::ofstream aux(filename.str(), std::ios::out | std::ios::trunc | std::ios::binary);
  if (aux.is_open()) {
    aux.write((char *) magnetization, size*sizeof(Real));
    aux.write((char *) thermal, size*sizeof(Real));
    aux.close();
  }
  else {
    throw std::ios_base::failure("There was an error writing the benchmark.");
  }
}
