//! \file coldeos_compose.cpp
//  \brief Implementation of ColdColdEOSCompOSE

#include <cassert>
#include <cmath>
#include <cstdio>
#include <sstream>
#include <fstream>
#include <iomanip>

#include <hdf5.h>
#include <hdf5_hl.h>

#include "../include/unit_system.hpp"
#include "../include/cold_eos_compose.hpp"

using namespace Primitive;
using namespace std;

#define MYH5CHECK(ierr) \
  if(ierr < 0) { \
    stringstream ss; \
    ss << __FILE__ << ":" << __LINE__ << " error reading EOS table!"; \
    throw runtime_error(ss.str().c_str()); \
  }

ColdEOSCompOSE::ColdEOSCompOSE():
  m_np(0),
  m_table(nullptr),
  m_initialized(false) {
  n_species = NSCALARS;
  eos_units = &Nuclear;
}

ColdEOSCompOSE::~ColdEOSCompOSE() {
  if (m_initialized) {
    delete[] m_table;
  }
}

Real ColdEOSCompOSE::Pressure(Real n) {
  assert (m_initialized);
  return exp(eval_at_n(ECLOGP, n));
}

Real ColdEOSCompOSE::Energy(Real n) {
  assert (m_initialized);
  return exp(eval_at_n(ECLOGE, n));
}

Real ColdEOSCompOSE::dPdn(Real n) {
  assert (m_initialized);
  return eval_at_n(ECDPDN, n);
}

Real ColdEOSCompOSE::SpecificInternalEnergy(Real n) {
  return Energy(n)/(mb*n) - 1;
}

Real ColdEOSCompOSE::Y(Real n, int iy) {
  assert (m_initialized);
  return eval_at_n(ECY+iy, n);
}

Real ColdEOSCompOSE::Enthalpy(Real n) {
  return (Pressure(n) + Energy(n))/n;
}

void ColdEOSCompOSE::ReadColdSliceFromFile(std::string fname, std::string species_names[NSCALARS]) {
  herr_t ierr;
  hid_t file_id;
  hid_t grp_id;
  hsize_t snb;

  // Open input file
  // -------------------------------------------------------------------------
  file_id = H5Fopen(fname.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    MYH5CHECK(file_id);

  // Open the cold_slice group
  grp_id = H5Gopen(file_id, "cold_slice", H5P_DEFAULT);
    MYH5CHECK(grp_id);

  // Get dataset sizes
  // -------------------------------------------------------------------------
  ierr = H5LTget_dataset_info(grp_id, "nb", &snb, NULL, NULL);
    MYH5CHECK(ierr);
  m_np = snb;

  // Allocate memory
  // -------------------------------------------------------------------------
  m_table = new Real[ECNVARS*m_np];
  double * scratch = new double[m_np];

  // Read nb, t, yq
  // -------------------------------------------------------------------------
  ierr = H5LTread_dataset_double(grp_id, "nb", scratch);
    MYH5CHECK(ierr);
  min_n = scratch[0];
  max_n = scratch[m_np-1];
  for (int in = 0; in < m_np; ++in) {
    m_table[index(ECLOGN, in)] = log(scratch[in]);
  }
  m_id_log_nb = 1.0/(m_table[index(ECLOGN, 1)] - m_table[index(ECLOGN, 0)]);

  ierr = H5LTread_dataset_double(grp_id, "t", scratch);
    MYH5CHECK(ierr);
  T = scratch[0];
  // the neutron mass is used as the baryon mass in CompOSE
  ierr = H5LTread_dataset_double(grp_id, "mn", scratch);
    MYH5CHECK(ierr);
  mb = scratch[0];

  // Read other thermodynamics quantities
  // -------------------------------------------------------------------------
  ierr = H5LTread_dataset_double(grp_id, "Q1", scratch);
    MYH5CHECK(ierr);
  for (int in = 0; in < m_np; ++in) {
    m_table[index(ECLOGP, in)] =
      log(scratch[in]) + m_table[index(ECLOGN, in)];
  }

  ierr = H5LTread_dataset_double(grp_id, "Q7", scratch);
    MYH5CHECK(ierr);
  for (int in = 0; in < m_np; ++in) {
    m_table[index(ECLOGE, in)] =
      log(mb*(scratch[in] + 1)) + m_table[index(ECLOGN, in)];
  }

  for (int iy = 0; iy < n_species; ++iy) {
    std::stringstream ss;
    ss << "Y[" << species_names[iy] << "]";
    ierr = H5LTread_dataset_double(grp_id, ss.str().c_str(), scratch);
      MYH5CHECK(ierr);
    for (int in = 0; in < m_np; ++in) {
      m_table[index(ECY+iy, in)] = scratch[in];
    }
  }

  //  fill enthalpy
  for (int in = 0; in < m_np; ++in) {
    m_table[index(ECH, in)] = (exp(m_table[index(ECLOGE, in)]) + exp(m_table[index(ECLOGP, in)]))/m_table[index(ECLOGN, in)];
  }

  //  fill dPdn
  D0_x_2(&m_table[index(ECLOGP, 0)], &m_table[index(ECLOGN, 0)], m_np, &m_table[index(ECDPDN, 0)]);
  for (int in = 1; in < m_np; ++in) {
    m_table[index(ECDPDN, in)] *= exp(m_table[index(ECLOGP, in)] - m_table[index(ECLOGN, in)]);
  }

  // Cleanup
  // -------------------------------------------------------------------------
  delete[] scratch;
  H5Fclose(file_id);

  m_initialized = true;
}

void ColdEOSCompOSE::DumpLoreneEOSFile(std::string fname) {
  // Dump the eos_akmalpr.d file that lorene routines expect
  // Lorene units are n [fm^-3], e [g/cm^3], p [erg/cm^3]
  Real n_conv = eos_units->DensityConversion(Nuclear);
  Real e_conv = eos_units->DensityConversion(CGS) * eos_units->MassConversion(CGS);
  Real p_conv = eos_units->PressureConversion(CGS);

  std::ofstream lorenefile(fname.c_str());
  lorenefile << std::scientific << std::setprecision(15);

  lorenefile << "#\n#\n#\n#\n#\n" << m_np << "\n#\n#\n#\n";

  for (int i = 0; i < m_np; ++i) {
    Real nb = n_conv * exp(m_table[index(ECLOGN, i)]);
    Real e = e_conv * exp(m_table[index(ECLOGE, i)]);
    Real p = p_conv * exp(m_table[index(ECLOGP, i)]);
    lorenefile << i << " " << nb << " " << e << " " << p << std::endl;
  }
}

void ColdEOSCompOSE::weight_idx_ln(Real *w0, Real *w1, int *in, Real log_n) const {
  *in = (log_n - m_table[index(ECLOGN, 0)])*m_id_log_nb;
  *w1 = (log_n - m_table[index(ECLOGN, *in)])*m_id_log_nb;
  *w0 = 1.0 - (*w1);
}

Real ColdEOSCompOSE::eval_at_ln(int iv, Real log_n) const {
  int in;
  Real wn0, wn1;
  weight_idx_ln(&wn0, &wn1, &in, log_n);
  return wn0 * m_table[index(iv, in+0)] + wn1 * m_table[index(iv, in+1)];
}

Real ColdEOSCompOSE::eval_at_n(int iv, Real n) const {
  return eval_at_ln(iv, log(n));
}

Real ColdEOSCompOSE::eval_at_general(int ii, int iv, Real h) const {
  throw std::logic_error("ColdEOSCompOSE::eval_at_general not implemented");
}


//--------------------------------------------------------------------------------------
//! \fn int D0_x_2(double *f, double *x, int n, double *df)
// \brief 1st order centered stencil first derivative, nonuniform grids
int ColdEOSCompOSE::D0_x_2(double *f, double *x, int n, double *df)
{
  int i;
  for(i=1; i<n-1; i++) {
    df[i] = (f[i+1]-f[i-1])/(x[i+1]-x[i-1]);
  }
  i = 0;
  df[i] = (f[i]-f[i+1])/(x[i]-x[i+1]);
  i = n-1;
  df[i] = (f[i]-f[i-1])/(x[i]-x[i-1]);
  return 0;
}
