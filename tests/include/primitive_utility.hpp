#ifndef PRIMITIVE_UTILITY_HPP
#define PRIMITIVE_UTILITY_HPP
//! \file primitive_utility.hpp
//  \brief Test utilities specific to primitive solvers.

#include <ps_types.hpp>

void MinkowskiMetric(Real gd[NSPMETRIC], Real gu[NSPMETRIC]);

void SchwarzschildMetric(Real gd[NSPMETRIC], Real gu[NSPMETRIC]);

void ScrewballMinkowskiMetric(Real gd[NSPMETRIC], Real gu[NSPMETRIC]);

void ScrewballSchwarzschildMetric(Real gd[NSPMETRIC], Real gu[NSPMETRIC]);

void ZeroVelocity(Real prim[NPRIM]);

void StrongVelocity(Real prim[NPRIM]);

void ZeroField(Real bu[NMAG]);

void StrongField(Real bu[NMAG]);

void ParticleFractions(Real prim[NPRIM], int s);

#endif
