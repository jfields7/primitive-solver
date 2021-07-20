//! \file primitive_utility.hpp
//  \brief Test utilities specific to primitive solvers.

#include <ps_types.hpp>

void MinkowskiMetric(Real gd[NMETRIC], Real gu[NMETRIC]);

void SchwarzschildMetric(Real gd[NMETRIC], Real gu[NMETRIC]);

void ScrewballMinkowskiMetric(Real gd[NMETRIC], Real gu[NMETRIC]);

void ScrewballSchwarzschildMetric(Real gd[NMETRIC], Real gu[NMETRIC]);

void ZeroVelocity(Real prim[NPRIM]);

void StrongVelocity(Real prim[NPRIM]);

void ZeroField(Real bu[NMAG]);

void StrongField(Real bu[NMAG]);

void ParticleFractions(Real prim[NPRIM], int s);
