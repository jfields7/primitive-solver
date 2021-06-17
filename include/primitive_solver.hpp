#ifndef PRIMITIVE_SOLVER_HPP
#define PRIMITIVE_SOLVER_HPP

/*!
 * \file primitive_solver.hpp
 * \brief contains functions for retrieving primitive variables from conserved variables.
 */

#include <ps_types.hpp>

// Forward declaration
template <typename T> class AthenaArray;

namespace PrimitiveSolver{
  //bool con_to_prim(AthenaArray<Real> &cons, AthenaArray<Real> &prim, AthenaArray<

  // Master function
  Real f(Real x);

  // Upper bound function
  void f_upper(Real &f, Real &df, Real mu, Real h0, Real rsq, Real bsq, Real rbsq);

  // Help functions
  double rbarsq(Real mu, Real rsq, Real rbsq, Real x);
  double drbarsq(Real mu, Real rsq, Real bsq, Real rbsq, Real x);
  double xf(Real mu, Real bsq);
  double dxf(Real x, Real bsq);
};

#endif
