#ifndef NUMTOOLS_ROOT_HPP
#define NUMTOOLS_ROOT_HPP

//! \file numtools_root.hpp
//  \author Jacob Fields
//
//  \brief Declares some functions for root-finding.

#include <ps_types.hpp>
#include <cmath>

namespace NumTools {

class Root {
  private:
  public:
    /// Maximum number of iterations
    unsigned int iterations;
    /// Solver tolerance
    Real tol;
    /// Only used for benchmarking, not thread-safe.
    int last_count;

    Root() : iterations(30), tol(1e-15) {}

    Root(Root const&) = delete;
    void operator=(Root const&) = delete;

    // FalsePosition {{{

    //! \brief Find the root of a functor f using false position.
    //
    // Find the root of a generic functor taking at least one argument. The first
    // argument is assumed to be the quantity of interest. All other arguments are
    // assumed to be constant parameters for the function. The root-finding method
    // is the Illinois variant of false position.
    //
    // \param[in]  f  The functor to find a root for. Its root function must take at
    //                least one argument.
    // \param[in,out]  lb  The lower bound for the root.
    // \param[in,out]  ub  The upper bound for the root.
    // \param[out]  x  The location of the root.
    // \param[in]  args  Additional arguments required by f.

    template<class Functor, class ... Types>
    inline bool FalsePosition(Functor&& f, Real &lb, Real &ub, Real& x, Types ... args) {
      int side = 0;
      Real ftest;
      unsigned int count = 0;
      last_count = 0;
      // Get our initial bracket.
      Real flb = f(lb, args...);
      Real fub = f(ub, args...);
      Real xold;
      x = lb;
      // If one of the bounds is already within tolerance of the root, we can skip all of this.
      if (std::fabs(flb) <= tol) {
        x = lb;
        return true;
      }
      else if (std::fabs(fub) <= tol) {
        x = ub;
        return true;
      }
      if (flb*fub > 0) {
        return false;
      }
      do {
        xold = x;
        // Calculate the new root position.
        x = (fub*lb - flb*ub)/(fub - flb);
        count++;
        // Calculate f at the prospective root.
        ftest = f(x,args...);
        if (std::fabs((x-xold)/x) <= tol) {
          break;
        }
        // Check the sign of f. If f is on the same side as the lower bound, then we adjust
        // the lower bound. Similarly, if f is on the same side as the upper bound, we 
        // adjust the upper bound. If ftest falls on the same side twice, we weight one of
        // the sides to force the new root to fall on the other side. This allows us to
        // whittle down both sides at once and get better average convergence.
        if (ftest*flb >= 0) {
          flb = ftest;
          lb = x;
          if (side == 1) {
            fub /= 2.0;
          }
          side = 1;
        }
        else {
          fub = ftest;
          ub = x;
          if (side == -1) {
            flb /= 2.0;
          }
          side = -1;
        }
      }
      while (count < iterations);
      last_count = count;

      // Return success if we're below the tolerance, otherwise report failure.
      return fabs((x-xold)/x) <= tol;
    }
    
    // }}}

    // Chandrupatla {{{

    //! \brief Find the root of a functor f using Chandrupatla's method
    //
    // Find the root of a generic functor taking at least one argument. The first
    // argument is assumed to be the quantity of interest. All other arguments are
    // assumed to be constant parameters for the function. The root-finding method
    // is Chandrupatla's method, a simpler alternative to Brent's method with
    // comparable performance.
    //
    // \param[in]  f  The functor to find a root for. Its root function must take at
    //                least one argument.
    // \param[in,out]  lb  The lower bound for the root.
    // \param[in,out]  ub  The upper bound for the root.
    // \param[out]  x  The location of the root.
    // \param[in]  args  Additional arguments required by f.

    template<class Functor, class ... Types>
    inline bool Chandrupatla(Functor&& f, Real &lb, Real &ub, Real& x, Types ... args) {
      unsigned int count = 0;
      last_count = 0;
      // Get our initial bracket.
      Real flb = f(lb, args...);
      Real fub = f(ub, args...);
      x = lb;
      // If one of the bounds is already within tolerance of the root, we can skip all of this.
      if (std::fabs(flb) <= tol) {
        x = lb;
        return true;
      }
      else if (std::fabs(fub) <= tol) {
        x = ub;
        return true;
      }
      // Make sure the bracket is valid
      if (flb*fub > 0) {
        return false;
      }
      Real t = 0.5;
      Real x1, x2, x3, f1, f2, f3, ftest;
      Real phi1, xi1;
      x1 = ub;
      x2 = lb;
      f1 = fub;
      f2 = flb;
      do {
        // Estimate the new root position
        x = x1 + t*(x2 - x1);
        count++;
        // Calculate f at the prospective root
        ftest = f(x, args...);
        if (std::fabs((x-x1)/x) <= tol) {
          break;
        }
        // Check the sign of ftest to determine the new bounds
        if (ftest*f1 >= 0) {
          x3 = x1;
          x1 = x;
          f3 = f1;
          f1 = ftest;
        }
        else {
          x3 = x2;
          x2 = x1;
          x1 = x;
          f3 = f2;
          f2 = f1;
          f1 = ftest;
        }
        // Check if we're in the region of validity for quadratic interpolation.
        phi1 = (f1 - f2)/(f3 - f2);
        xi1 = (x1 - x2)/(x3 - x2);
        if (1.0 - std::sqrt(1.0 - xi1) < phi1 && phi1 < std::sqrt(xi1)) {
          // Perform quadratic interpolation
          t = f1/(f3 - f2)*(f3/(f1 - f2) + (x3 - x1)/(x2 - x1)*f2/(f3 - f1));
        }
        else {
          // Perform bisection instead
          t = 0.5;
        }
      }
      while (count < iterations);
      last_count = count;

      // Return success if we're below the tolerance, otherwise report failure.
      return fabs((x-x1)/x) <= tol;
    }

    // }}}

    // NewtonSafe {{{
    /*! \brief Find the root of a function f using a safe Newton solve.
     *
     * A safe Newton solve performs a Newton-Raphson solve, but it also brackets the
     * root using bisection to ensure that the result converges.
     *
     * \param[in]     f     The functor to find a root for. Its root function must take
     *                      at least one argument.
     * \param[in,out] lb    The lower bound for the root of f.
     * \param[in,out] ub    The upper bound for the root of f.
     * \param[out]    x     The root of f.
     * \param[in]     args  Additional arguments required by f.
     */
    template<class Functor, class ... Types>
    inline bool NewtonSafe(Functor&& f, Real &lb, Real &ub, Real& x, Types ... args) {
      Real fx;
      Real dfx;
      Real xold;
      unsigned int count = 0;
      //last_count = 0;
      // We first need to ensure that the bracket is valid.
      Real fub, flb;
      f(flb, dfx, lb, args...);
      f(fub, dfx, ub, args...);
      if (flb*fub > 0) {
        return 0;
      }
      // If one of the roots is already within tolerance, then
      // we don't need to do the solve.
      if (std::fabs(flb) <= tol) {
        x = lb;
        return true;
      }
      else if (std::fabs(fub) <= tol) {
        x = ub;
        return true;
      }
      // Since we already had to evaluate the function at the bounds,
      // we can predict our starting position using false position.
      x = (fub*lb - flb*ub)/(fub - flb);
      do {
        xold = x;
        // Calculate f and df at point x.
        f(fx, dfx, x, args...);
        x = x - fx/dfx;
        // Check that the root is bounded properly.
        if (x > ub || x < lb) {
          // Revert to bisection if the root is not converging.
          x = 0.5*(ub + lb);
          //f(fx, dfx, x, args...);
        }
        // Correct the bounds.
        if (fx*flb > 0) {
          flb = fx;
          lb = xold;
        }
        else if (fx*fub > 0) {
          fub = fx;
          ub = xold;
        }
        count++;
      }
      while (std::fabs((xold-x)/x) > tol && count < iterations);
      //last_count = count;

      // Return success if we're below the tolerance, otherwise report failure.
      return std::fabs((x-xold)/x) <= tol;
    }
    // }}}
};

} // namespace

#endif
