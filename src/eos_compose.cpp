//! \file eos_compose.cpp
//  \brief Implementation of EOSCompose

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <limits>
#include <sstream>
#include <stdexcept>

#include <hdf5.h>
#include <hdf5_hl.h>

#include <eos_compose.hpp>
#include <numtools_root.hpp>
#include <unit_system.hpp>

using namespace Primitive;
using namespace std;

#define MYH5CHECK(ierr) \
  if(ierr < 0) { \
    stringstream ss; \
    ss << __FILE__ << ":" << __LINE__ << " error reading EOS table!"; \
    throw runtime_error(ss.str().c_str()); \
  }

EOSCompOSE::EOSCompOSE():
  m_id_log_nb(numeric_limits<Real>::quiet_NaN()),
  m_id_log_t(numeric_limits<Real>::quiet_NaN()),
  m_id_yq(numeric_limits<Real>::quiet_NaN()),
  m_nn(0), m_nt(0), m_ny(0),
  m_min_h(numeric_limits<Real>::max()),
  m_log_nb(nullptr),
  m_log_t(nullptr),
  m_yq(nullptr),
  m_table(nullptr),
  m_initialized(false) {
  n_species = 1;
  eos_units = &Nuclear;
}

EOSCompOSE::~EOSCompOSE() {
  if (m_initialized) {
    delete[] m_log_nb;
    delete[] m_log_t;
    delete[] m_yq;
    delete[] m_table;
  }
}

Real EOSCompOSE::TemperatureFromE(Real n, Real e, Real *Y) {
  assert (m_initialized);
  return temperature_from_var(ECLOGE, log(e), n, Y[0]);
}

Real EOSCompOSE::TemperatureFromP(Real n, Real p, Real *Y) {
  assert (m_initialized);
  return temperature_from_var(ECLOGP, log(p), n, Y[0]);
}

Real EOSCompOSE::Energy(Real n, Real T, Real *Y) {
  assert (m_initialized);
  return exp(eval_at_nty(ECLOGE, n, T, Y[0]));
}

Real EOSCompOSE::Pressure(Real n, Real T, Real *Y) {
  assert (m_initialized);
  return exp(eval_at_nty(ECLOGP, n, T, Y[0]));
}

Real EOSCompOSE::Entropy(Real n, Real T, Real *Y) {
  assert (m_initialized);
  return eval_at_nty(ECENT, n, T, Y[0]);
}

Real EOSCompOSE::Enthalpy(Real n, Real T, Real *Y) {
  assert (m_initialized);
  Real const P = Pressure(n, T, Y);
  Real const e = Energy(n, T, Y);
  return (P + e)/n;
}

Real EOSCompOSE::SoundSpeed(Real n, Real T, Real *Y) {
  assert (m_initialized);
  return eval_at_nty(ECCS, n, T, Y[0]);
}

Real EOSCompOSE::SpecificInternalEnergy(Real n, Real T, Real *Y) {
  return Energy(n, T, Y)/(mb*n) - 1;
}

Real EOSCompOSE::BaryonChemicalPotential(Real n, Real T, Real *Y) {
  assert (m_initialized);
  return eval_at_nty(ECMUB, n, T, Y[0]);
}

Real EOSCompOSE::ChargeChemicalPotential(Real n, Real T, Real *Y) {
  assert (m_initialized);
  return eval_at_nty(ECMUQ, n, T, Y[0]);
}

Real EOSCompOSE::ElectronLeptonChemicalPotential(Real n, Real T, Real *Y) {
  assert (m_initialized);
  return eval_at_nty(ECMUL, n, T, Y[0]);
}

int EOSCompOSE::BetaEquilibriumTrapped(Real n, Real e, Real *Yl, Real &T_eq, Real *Y_eq, Real T_guess, Real *Y_guess) {

  const int n_at = 16;
  Real vec_guess[n_at][2] = { 
    {1.00e0, 1.00e0},
    {0.90e0, 1.25e0},
    {0.90e0, 1.10e0},
    {0.90e0, 1.00e0},
    {0.90e0, 0.90e0},
    {0.90e0, 0.75e0},
    {0.75e0, 1.25e0},
    {0.75e0, 1.10e0},
    {0.75e0, 1.00e0},
    {0.75e0, 0.90e0},
    {0.75e0, 0.75e0},
    {0.50e0, 1.25e0},
    {0.50e0, 1.10e0},
    {0.50e0, 1.00e0},
    {0.50e0, 0.90e0},
    {0.50e0, 0.75e0},
  };


  // ierr = 0    Equilibrium found
  // ierr = 1    Equilibrium not found
  int ierr = 1;
  int na = 0; // counter for the number of attempts

  Real x0[2], x1[2]; // T,Ye guess and T,Ye result

  while (ierr!=0 and na<n_at) {
    x0[0] = vec_guess[na][0] * T_guess;
    x0[1] = vec_guess[na][1] * Y_guess[0];

    ierr = trapped_equilibrium_2DNR(n, e, Yl[0], x0, x1);

    na += 1;
  }

  if (ierr==0){ // Success
    T_eq = x1[0];
    Y_eq[0] = x1[1];
  } else { // Failure
    T_eq = T_guess;       // Set results to guesses
    Y_eq[0] = Y_guess[0];
  }

  return ierr;
}

int EOSCompOSE::trapped_equilibrium_2DNR(Real n, Real e, Real Yle, Real x0[2], Real x1[2]) {
  
  const Real eps_lim  = 1.e-14; // tolerance in 2D NR (required for 1e-12 err in T)
  const int n_max     = 100;    // Newton-Raphson max number of iterations
  const int n_cut_max = 8;      // Bisection max number of iterations
  int ierr = 0;

  // initialize the solution
  x1[0] = x0[0];
  x1[1] = x0[1];
  bool KKT = false;

  //compute the initial residuals
  Real y[2] = {0.0};
  func_eq_weak(n,e,Yle,x1,y);

  // compute the error from the residuals
  Real err = error_func_eq_weak(Yle,e,y);

  // initialize the iteration variables
  int n_iter = 0;
  Real J[2][2] = {0.0};
  Real invJ[2][2] = {0.0};
  Real dx1[2] = {0.0};
  Real dxa[2] = {0.0};
  Real norm[2] = {0.0};
  Real x1_tmp[2] = {0.0};

  // loop until a low enough residual is found or until  a too
  // large number of steps has been performed
  while (err>eps_lim && n_iter<=n_max && !KKT) {
    // compute the Jacobian
    ierr = jacobi_eq_weak(n,e,Yle,x1,J);
    if (ierr != 0) {
      return ierr;
    }

    // compute and check the determinant of the Jacobian
    Real det = J[0][0]*J[1][1] - J[0][1]*J[1][0];
    if (det==0.0) {
      ierr = 1;
      return ierr;
    }

    // invert the Jacobian
    inv_jacobi(det,J,invJ);

    // compute the next step
    dx1[0] = - (invJ[0][0]*y[0] + invJ[0][1]*y[1]);
    dx1[1] = - (invJ[1][0]*y[0] + invJ[1][1]*y[1]);

    // check if we are the boundary of the table
    if (x1[0] == min_T) {
      norm[0] = -1.0;
    } else if (x1[0] == max_T) {
      norm[0] = 1.0;
    } else { 
      norm[0] = 0.0;
    }

    if (x1[1] == min_Y[0]) {
      norm[1] = -1.0;
    } else if (x1[1] == max_Y[0]) {
      norm[1] = 1.0;
    } else {
      norm[1] = 0.0;
    }

    // Take the part of the gradient that is active (pointing within the eos domain)
    Real scal = norm[0]*norm[0] + norm[1]*norm[1];
    if (scal <= 0.5) { // this can only happen if norm = (0, 0)
      scal = 1.0;
    }
    dxa[0] = dx1[0] - (dx1[0]*norm[0] + dx1[1]*norm[1])*norm[0]/scal;
    dxa[1] = dx1[1] - (dx1[0]*norm[0] + dx1[1]*norm[1])*norm[1]/scal;

    if ((dxa[0]*dxa[0] + dxa[1]*dxa[1]) < (eps_lim*eps_lim * (dx1[0]*dx1[0] + dx1[1]*dx1[1]))) {
      KKT = true;
      ierr = 2;
      return ierr;
    }

    int n_cut = 0;
    Real fac_cut = 1.0;
    Real err_old = err;

    while (n_cut <= n_cut_max && err >= err_old) {
      // the variation of x1 is divided by an powers of 2 if the
      // error is not decreasing along the gradient direction
      
      x1_tmp[0] = x1[0] + (dx1[0]*fac_cut);
      x1_tmp[1] = x1[1] + (dx1[1]*fac_cut);

      // check if the next step calculation had problems
      if (isnan(x1_tmp[0])) {
        ierr = 1;
        return ierr;
      }

      // tabBoundsFlag = enforceTableBounds(rho, x1_tmp[0], x1_tmp[1]);
      x1_tmp[0] = min(max(x1_tmp[0],min_T),max_T);
      x1_tmp[1] = min(max(x1_tmp[1],min_Y[0]),max_Y[0]);

      // assign the new point
      x1[0] = x1_tmp[0];
      x1[1] = x1_tmp[1];

      // compute the residuals for the new point
      func_eq_weak(n,e,Yle,x1,y);

      // compute the error
      err = error_func_eq_weak(Yle,e,y);

      // update the bisection cut along the gradient
      n_cut += 1;
      fac_cut *= 0.5;
    }

    // update the iteration
    n_iter += 1;
  }
    
  if (n_iter <= n_max) {
    ierr = 0;
  } else {
    ierr = 1;
  }
  
  return ierr;
}

void EOSCompOSE::func_eq_weak(Real n, Real e_eq, Real Yle, Real x[2], Real y[2]) {
  Real T = x[0];

  Real Y[MAX_SPECIES] = {0.0};
  Y[0] = x[1];

  Real mu_l = ElectronLeptonChemicalPotential(n, T, Y);
  Real e = Energy(n, T, Y);
  Real eta = mu_l/T;
  Real eta2 = eta*eta;

  Real t3 = T*T*T;
  Real t4 = t3*T;
  y[0] = Y[0] + pref1*t3*eta*(pi2 + eta2)/n - Yle;
  y[1] = (e+pref2*t4*((cnst5+0.5*eta2*(pi2+0.5*eta2))+cnst6))/e_eq - 1.0;

  return;

}

Real EOSCompOSE::error_func_eq_weak(Real Yle, Real e_eq, Real y[2]) {
  Real err = abs(y[0]/Yle) + abs(y[1]/1.0);
  return err;
}

int EOSCompOSE::jacobi_eq_weak(Real n, Real e_eq, Real Yle, Real x[2], Real J[2][2]) {
  int ierr = 0;

  Real T = x[0];
  Real Y[MAX_SPECIES] = {0.0};
  Y[0] = x[1];

  if (isnan(T)) {
    ierr = 1;
    return ierr;
  }

  Real mu_l = ElectronLeptonChemicalPotential(n, T, Y);
  Real eta = mu_l/T;
  Real eta2 = eta*eta;

  if (isnan(eta)) {
    ierr = 1;
    return ierr;
  }

  Real detadt,detadye,dedt,dedye;
  ierr = eta_e_gradient(n,T,Y,eta,detadt,detadye,dedt,dedye);
  if (ierr != 0){
    return ierr;
  }

  Real T2 = T*T;
  Real T3 = T2*T;
  // Real T4 = T3*T;

  J[0][0] = pref1/n*T2*(3.e0*eta*(pi2+eta2)+T*(pi2+3.e0*eta2)*detadt);
  J[0][1] = 1.e0+pref1/n*T3*(pi2+3.e0*eta2)*detadye;

  J[1][0] = (dedt+pref2*T3*(cnst3+cnst4+2.e0*eta2*(pi2+0.5*eta2)+eta*T*(pi2+eta2)*detadt))/e_eq;
  J[1][1] = (dedye+pref2*T3*eta*(pi2+eta2)*detadye)/e_eq;

  return ierr;
}

int EOSCompOSE::eta_e_gradient(Real n, Real T, Real *Y, Real eta, Real &deta_dT, Real &deta_dYe, Real &de_dT, Real &de_dYe) {
  int ierr=0;

  const Real Ye_delta = 0.005;
  const Real T_delta = 0.01;

  Real Y1[MAX_SPECIES] = {0.0};
  Real Y2[MAX_SPECIES] = {0.0};

  Y1[0] = max(Y[0] - Ye_delta, min_Y[0]);
  Real mu_l1 = ElectronLeptonChemicalPotential(n, T, Y1);
  Real e1 = Energy(n, T, Y1);
  
  Y2[0] = min(Y[0] + Ye_delta, max_Y[0]);
  Real mu_l2 = ElectronLeptonChemicalPotential(n, T, Y2);
  Real e2 = Energy(n, T, Y2);

  Real dmu_l_dYe = (mu_l2-mu_l1)/(Y2[0] - Y1[0]);
  de_dYe         = (e2-e1)/(Y2[0] - Y1[0]);

  Real T1 = max(T - T_delta, min_T);
  mu_l1 = ElectronLeptonChemicalPotential(n, T1, Y);
  e1 = Energy(n, T1, Y);

  Real T2 = min(T + T_delta, max_T);
  mu_l2 = ElectronLeptonChemicalPotential(n, T2, Y);
  e2 = Energy(n, T2, Y);

  Real dmu_l_dT   = (mu_l2 - mu_l1)/(T2 - T1);
  de_dT          = (e2 - e1)/(T2 - T1);
  
  deta_dT  = (dmu_l_dT - eta )/T; // [1/MeV] TODO: Check
  deta_dYe = dmu_l_dYe/T;      // [-]

  if (isnan(deta_dT)||isnan(deta_dYe)||isnan(de_dT)||isnan(de_dYe)) {
    ierr = 1;
  }

  return ierr;
}

void EOSCompOSE::inv_jacobi(Real det, Real J[2][2], Real invJ[2][2]) {
  Real inv_det = 1.0/det;
  invJ[0][0] =  J[1][1]*inv_det;
  invJ[1][1] =  J[0][0]*inv_det;
  invJ[0][1] = -J[0][1]*inv_det;
  invJ[1][0] = -J[1][0]*inv_det;
}


void EOSCompOSE::TrappedNeutrinos(Real n, Real T, Real *Y, Real n_nu[3], Real e_nu[3]) {
  Real mu_le = ElectronLeptonChemicalPotential(n, T, Y);
  Real eta_e = mu_le/T;
  Real eta_e2 = eta_e*eta_e;

  Real eta_m = 0.0;
  Real eta_m2 = 0.0;

  Real eta_t = 0.0;
  Real eta_t2 = 0.0;

  Real T3 = T*T*T;
  Real T4 = T3*T;

  n_nu[0] = pref1 * T3 * (eta_e * (pi2 + eta_e2)); // n_nu_e   - n_anu_e   [fm^-3]
  n_nu[1] = pref1 * T3 * (eta_m * (pi2 + eta_m2)); // n_nu_mu  - n_anu_mu  [fm^-3]
  n_nu[2] = pref1 * T3 * (eta_t * (pi2 + eta_t2)); // n_nu_tau - n_anu_tau [fm^-3]

  e_nu[0] = pref2*T4*(cnst5+0.5*eta_e2*(pi2+0.5*eta_e2)); // e_nu_e   + e_anu_e   [MeV fm^-3]
  e_nu[1] = pref2*T4*(cnst5+0.5*eta_m2*(pi2+0.5*eta_m2)); // e_nu_mu  + e_anu_mu  [MeV fm^-3]
  e_nu[2] = pref2*T4*(cnst5+0.5*eta_t2*(pi2+0.5*eta_t2)); // e_nu_tau + e_anu_tau [MeV fm^-3]

  return;
}

Real EOSCompOSE::MinimumEnthalpy() {
  return m_min_h;
}

Real EOSCompOSE::MinimumPressure(Real n, Real *Y) {
  return Pressure(n, min_T, Y);
}

Real EOSCompOSE::MaximumPressure(Real n, Real *Y) {
  return Pressure(n, max_T, Y);
}

Real EOSCompOSE::MinimumEnergy(Real n, Real *Y) {
  return Energy(n, min_T, Y);
}

Real EOSCompOSE::MaximumEnergy(Real n, Real *Y) {
  return Energy(n, max_T, Y);
}

void EOSCompOSE::ReadTableFromFile(std::string fname) {
  herr_t ierr;
  hid_t file_id;
  hsize_t snb, st, syq;

  // Open input file
  // -------------------------------------------------------------------------
  file_id = H5Fopen(fname.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    MYH5CHECK(file_id);

  // Get dataset sizes
  // -------------------------------------------------------------------------
  ierr = H5LTget_dataset_info(file_id, "nb", &snb, NULL, NULL);
    MYH5CHECK(ierr);
  ierr = H5LTget_dataset_info(file_id, "t", &st, NULL, NULL);
    MYH5CHECK(ierr);
  ierr = H5LTget_dataset_info(file_id, "yq", &syq, NULL, NULL);
    MYH5CHECK(ierr);
  m_nn = snb;
  m_nt = st;
  m_ny = syq;

  // Allocate memory
  // -------------------------------------------------------------------------
  m_log_nb = new Real[m_nn];
  m_log_t = new Real[m_nt];
  m_yq = new Real[m_ny];
  m_table = new Real[ECNVARS*m_nn*m_ny*m_nt];
  double * scratch = new double[m_nn*m_ny*m_nt];

  // Read nb, t, yq
  // -------------------------------------------------------------------------
  ierr = H5LTread_dataset_double(file_id, "nb", scratch);
    MYH5CHECK(ierr);
  min_n = scratch[0];
  max_n = scratch[m_nn-1];
  for (int in = 0; in < m_nn; ++in) {
    m_log_nb[in] = log(scratch[in]);
  }
  m_id_log_nb = 1.0/(m_log_nb[1] - m_log_nb[0]);

  ierr = H5LTread_dataset_double(file_id, "t", scratch);
    MYH5CHECK(ierr);
  min_T = scratch[1];
  max_T = scratch[m_nt-1];
  for (int it = 0; it < m_nt; ++it) {
    m_log_t[it] = log(scratch[it]);
  }
  m_id_log_t = 1.0/(m_log_t[1] - m_log_t[0]);

  ierr = H5LTread_dataset_double(file_id, "yq", scratch);
    MYH5CHECK(ierr);
  min_Y[0] = scratch[0];
  max_Y[0] = scratch[m_ny-1];
  for (int iy = 0; iy < m_ny; ++iy) {
    m_yq[iy] = scratch[iy];
  }
  m_id_yq = 1.0/(m_yq[1] - m_yq[0]);

  // the neutron mass is used as the baryon mass in CompOSE
  ierr = H5LTread_dataset_double(file_id, "mn", scratch);
    MYH5CHECK(ierr);
  mb = scratch[0];

  // Read other thermodynamics quantities
  // -------------------------------------------------------------------------
  ierr = H5LTread_dataset_double(file_id, "Q1", scratch);
    MYH5CHECK(ierr);
  for (int inb = 0; inb < m_nn; ++inb) {
  for (int iyq = 0; iyq < m_ny; ++iyq) {
  for (int it = 0; it < m_nt; ++it) {
    m_table[index(ECLOGP, inb, iyq, it)] =
        log(scratch[index(0, inb, iyq, it)]) + m_log_nb[inb];
  }}}

  ierr = H5LTread_dataset_double(file_id, "Q2", scratch);
    MYH5CHECK(ierr);
  copy(&scratch[0], &scratch[m_nn*m_ny*m_nt], &m_table[index(ECENT, 0, 0, 0)]);

  ierr = H5LTread_dataset_double(file_id, "Q3", scratch);
    MYH5CHECK(ierr);
  for (int in = 0; in < m_nn; ++in) {
  for (int iy = 0; iy < m_ny; ++iy) {
  for (int it = 0; it < m_nt; ++it) {
    m_table[index(ECMUB, in, iy, it)] =
      mb*(scratch[index(0, in, iy, it)] + 1);
  }}}

  ierr = H5LTread_dataset_double(file_id, "Q4", scratch);
    MYH5CHECK(ierr);
  for (int in = 0; in < m_nn; ++in) {
  for (int iy = 0; iy < m_ny; ++iy) {
  for (int it = 0; it < m_nt; ++it) {
    m_table[index(ECMUQ, in, iy, it)] = mb*scratch[index(0, in, iy, it)];
  }}}

  ierr = H5LTread_dataset_double(file_id, "Q5", scratch);
    MYH5CHECK(ierr);
  for (int in = 0; in < m_nn; ++in) {
  for (int iy = 0; iy < m_ny; ++iy) {
  for (int it = 0; it < m_nt; ++it) {
    m_table[index(ECMUL, in, iy, it)] = mb*scratch[index(0, in, iy, it)];
  }}}

  ierr = H5LTread_dataset_double(file_id, "Q7", scratch);
    MYH5CHECK(ierr);
  for (int in = 0; in < m_nn; ++in) {
  for (int iy = 0; iy < m_ny; ++iy) {
  for (int it = 0; it < m_nt; ++it) {
    m_table[index(ECLOGE, in, iy, it)] =
      log(mb*(scratch[index(0, in, iy, it)] + 1)) + m_log_nb[in];
  }}}

  ierr = H5LTread_dataset_double(file_id, "cs2", scratch);
    MYH5CHECK(ierr);
  for (int in = 0; in < m_nn; ++in) {
  for (int iy = 0; iy < m_ny; ++iy) {
  for (int it = 0; it < m_nt; ++it) {
    m_table[index(ECCS, in, iy, it)] = sqrt(scratch[index(0, in, iy, it)]);
  }}}

  // Cleanup
  // -------------------------------------------------------------------------
  delete[] scratch;
  H5Fclose(file_id);

  m_initialized = true;

  // Compute minimum enthalpy
  // -------------------------------------------------------------------------
  for (int in = 0; in < m_nn; ++in) {
    Real const nb = exp(m_log_nb[in]);
    for (int it = 0; it < m_nt; ++it) {
      Real const t = exp(m_log_t[it]);
      for (int iy = 0; iy < m_ny; ++iy) {
        m_min_h = min(m_min_h, Enthalpy(nb, t, &m_yq[iy]));
      }
    }
  }
}

Real EOSCompOSE::temperature_from_var(int iv, Real var, Real n, Real Yq) const {
  int in, iy;
  Real wn0, wn1, wy0, wy1;
  weight_idx_ln(&wn0, &wn1, &in, log(n));
  weight_idx_yq(&wy0, &wy1, &iy, Yq);

  auto f = [=](int it){
    Real var_pt =
      wn0 * (wy0 * m_table[index(iv, in+0, iy+0, it)]  +
             wy1 * m_table[index(iv, in+0, iy+1, it)]) +
      wn1 * (wy0 * m_table[index(iv, in+1, iy+0, it)]  +
             wy1 * m_table[index(iv, in+1, iy+1, it)]);

    return var - var_pt;
  };

  int ilo = 0;
  int ihi = m_nt-1;
  Real flo = f(ilo);
  Real fhi = f(ihi);
  while (flo*fhi>0){
    if (ilo == ihi - 1) {
      break;
    } else {
      ilo += 1;
      flo = f(ilo);
    }
  }
  assert(flo*fhi <= 0);
  while (ihi - ilo > 1) {
    int ip = ilo + (ihi - ilo)/2;
    Real fp = f(ip);
    if (fp*flo <= 0) {
      ihi = ip;
      fhi = fp;
    }
    else {
      ilo = ip;
      flo = fp;
    }
  }
  assert(ihi - ilo == 1);
  Real lthi = m_log_t[ihi];
  Real ltlo = m_log_t[ilo];

  if (flo == 0) {
    return exp(ltlo);
  }
  if (fhi == 0) {
    return exp(lthi);
  }

  Real lt = m_log_t[ilo] - flo*(lthi - ltlo)/(fhi - flo);
  return exp(lt);
}

Real EOSCompOSE::eval_at_nty(int vi, Real n, Real T, Real Yq) const {
  return eval_at_lnty(vi, log(n), log(T), Yq);
}

void EOSCompOSE::weight_idx_ln(Real *w0, Real *w1, int *in, Real log_n) const {
  *in = (log_n - m_log_nb[0])*m_id_log_nb;
  *w1 = (log_n - m_log_nb[*in])*m_id_log_nb;
  *w0 = 1.0 - (*w1);
}

void EOSCompOSE::weight_idx_yq(Real *w0, Real *w1, int *iy, Real yq) const {
  *iy = (yq - m_yq[0])*m_id_yq;
  *w1 = (yq - m_yq[*iy])*m_id_yq;
  *w0 = 1.0 - (*w1);
}

void EOSCompOSE::weight_idx_lt(Real *w0, Real *w1, int *it, Real log_t) const {
  *it = (log_t - m_log_t[0])*m_id_log_t;
  *w1 = (log_t - m_log_t[*it])*m_id_log_t;
  *w0 = 1.0 - (*w1);
}

Real EOSCompOSE::eval_at_lnty(int iv, Real log_n, Real log_t, Real yq) const {
  int in, iy, it;
  Real wn0, wn1, wy0, wy1, wt0, wt1;

  weight_idx_ln(&wn0, &wn1, &in, log_n);
  weight_idx_yq(&wy0, &wy1, &iy, yq);
  weight_idx_lt(&wt0, &wt1, &it, log_t);

  return
    wn0 * (wy0 * (wt0 * m_table[index(iv, in+0, iy+0, it+0)]   +
                  wt1 * m_table[index(iv, in+0, iy+0, it+1)])  +
           wy1 * (wt0 * m_table[index(iv, in+0, iy+1, it+0)]   +
                  wt1 * m_table[index(iv, in+0, iy+1, it+1)])) +
    wn1 * (wy0 * (wt0 * m_table[index(iv, in+1, iy+0, it+0)]   +
                  wt1 * m_table[index(iv, in+1, iy+0, it+1)])  +
           wy1 * (wt0 * m_table[index(iv, in+1, iy+1, it+0)]   +
                  wt1 * m_table[index(iv, in+1, iy+1, it+1)]));
}
