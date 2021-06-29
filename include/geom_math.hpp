#ifndef GEOM_MATH_HPP
#define GEOM_MATH_HPP

//! \file geom_math.hpp
//  \brief Provides some simple tensor algebra functions.
//
//  The provided functions work for three dimensional
//  vectors, one-forms, and metrics.

#include <cmath>
#include <ps_types.hpp>

namespace Primitive {
  //! \brief Calculate the determinant from a 3-metric
  //
  //  \param[in] g3d The spatial metric
  //  \return The determinant of g3d
  inline Real GetDeterminant(const Real g3d[NSPMETRIC]) {
    return g3d[S11]*g3d[S22]*g3d[S33] + 2.0*g3d[S12]*g3d[S13]*g3d[S23] -
          (g3d[S11]*g3d[S23]*g3d[S23] + g3d[S12]*g3d[S12]*g3d[S33] + g3d[S13]*g3d[S13]*g3d[S22]);
  }

  //! \brief Convert a vector to a one-form
  //
  //  \param[out] vd  The output one-form
  //  \param[in]  vu  The input vector
  //  \param[in]  g3d The 3-metric
  inline void LowerVector(Real vd[3], const Real vu[3], const Real g3d[NSPMETRIC]) {
    vd[0] = g3d[S11]*vu[0] + g3d[S12]*vu[1] + g3d[S13]*vu[2];
    vd[1] = g3d[S12]*vu[0] + g3d[S22]*vu[1] + g3d[S23]*vu[2];
    vd[2] = g3d[S13]*vu[0] + g3d[S23]*vu[1] + g3d[S33]*vu[2];
  }

  //! \brief Convert a one-form to a vector
  //
  //  \param[out] vu  The output vector
  //  \param[in]  vd  The input one-form
  //  \param[in]  g3u The inverse 3-metric
  inline void RaiseForm(Real vu[3], const Real vd[3], const Real g3u[NSPMETRIC]) {
    vu[0] = g3u[S11]*vd[0] + g3u[S12]*vd[1] + g3u[S13]*vd[2];
    vu[1] = g3u[S12]*vd[0] + g3u[S22]*vd[1] + g3u[S23]*vd[2];
    vu[2] = g3u[S13]*vd[0] + g3u[S23]*vd[1] + g3u[S33]*vd[2];
  }

  //! \brief Contract a one-form with a vector
  //
  //  \param[in] au  The input vector
  //  \param[in] bd  The input one-form
  //  \return The contraction \f$a^u b_d\f$.
  inline Real Contract(const Real au[3], const Real bd[3]) {
    return au[0]*bd[0] + au[1]*bd[1] + au[2]*bd[2];
  }
}

#endif
