//! \file primitive_utility.cpp
//  \brief Implementation file for primitive_utility.hpp

#include <primitive_utility.hpp>

// MinkowskiMetric {{{
void MinkowskiMetric(Real gd[NSPMETRIC], Real gu[NSPMETRIC]) {
  for (int i = 0; i < NSPMETRIC; i++) {
    gu[i] = 0.0;
    gd[i] = 0.0;
  }
  
  gd[S11] = 1.0;
  gd[S22] = 1.0;
  gd[S33] = 1.0;

  gu[S11] = 1.0;
  gu[S22] = 1.0;
  gu[S33] = 1.0;
}
// }}}

// SchwarzschildMetric {{{
void SchwarzschildMetric(Real gd[NSPMETRIC], Real gu[NSPMETRIC]) {
  for (int i = 0; i < NMETRIC; i++) {
    gu[i] = 0.0;
    gd[i] = 0.0;
  }
  Real R = 1.0;
  Real rs = 1.0;
  Real hp = 1.0 + rs/(4.0*R);
  //Real hm = 1.0 - rs/(4.0*R);
  //Real gt = hm*hm/(hp*hp);
  Real gx = hp*hp*hp*hp;

  //gd[I00] = -gt;
  gd[S11] = gx;
  gd[S22] = gx;
  gd[S33] = gx;

  //gu[I00] = -1.0/gt;
  gu[S11] = 1.0/gx;
  gu[S22] = 1.0/gx;
  gu[S33] = 1.0/gx;
}
// }}}

// ScrewballMinkowskiMetric {{{
void ScrewballMinkowskiMetric(Real gd[NSPMETRIC], Real gu[NSPMETRIC]) {
  for (int i = 0; i < NMETRIC; i++) {
    gd[i] = 0.0;
    gu[i] = 0.0;
  }

  // Screwball coordinates:
  // a = t
  // b = x - y
  // c = x + y
  // d = z + x - 3y
  //gd[I00] = -1.0;
  gd[S11] = 9.0/2.0;
  gd[S12] = -2.0;
  gd[S13] = -2.0;
  gd[S22] = 3.0/2.0;
  gd[S23] = 1.0;
  gd[S33] = 1.0;

  //gu[I00] = 2.0;
  gu[S11] = 2.0;
  gu[S13] = 4.0;
  gu[S22] = 2.0;
  gu[S23] = -2.0;
  gu[S33] = 11.0;
}
// }}}

// ScrewballSchwarzschildMetric {{{
void ScrewballSchwarzschildMetric(Real gd[NSPMETRIC], Real gu[NSPMETRIC]) {
  for (int i = 0; i < NMETRIC; i++) {
    gd[i] = 0.0;
    gu[i] = 0.0;
  }
  Real R = 4.0;
  Real rs = 1.0;
  Real cv = 1.0 - rs/R;
  
  // Eddington-Finkelstein coordinates
  //gd[I00] = -cv;
  //gd[S01] = 2.0;
  gd[S22] = R*R;
  gd[S33] = R*R;

  //gu[I01] = 0.5;
  gu[S11] = cv/4.0;
  gu[S22] = 1.0/(R*R);
  gu[S33] = 1.0/(R*R);
}
// }}}

// ZeroVelocity {{{
void ZeroVelocity(Real prim[NPRIM]) {
  prim[IVX] = 0.0;
  prim[IVY] = 0.0;
  prim[IVZ] = 0.0;
}
// }}}

// StrongVelocity {{{
void StrongVelocity(Real prim[NPRIM]) {
  prim[IVX] = 55.0;
  prim[IVY] = 15.0;
  prim[IVZ] = 30.0;
}
// }}}

// ZeroField {{{
void ZeroField(Real bu[NMAG]) {
  bu[IB1] = 0.0;
  bu[IB2] = 0.0;
  bu[IB3] = 0.0;
}
// }}}

// StrongField {{{
void StrongField(Real bu[NMAG]) {
  bu[IB1] = 40.0;
  bu[IB2] = 60.0;
  bu[IB3] = 70.0;
}
// }}}

// ParticleFractions {{{
void ParticleFractions(Real prim[NPRIM], int s) {
  switch(s) {
    case 3:
      prim[IYF + 2] = 0.25;
    case 2:
      prim[IYF + 1] = 0.25;
    case 1:
      prim[IYF] = 0.25;
      break;
    case 0:
      break;
  }
}
// }}}
