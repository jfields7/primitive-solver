#ifndef NUMTOOLS_ROOT_HPP
#define NUMTOOLS_ROOT_HPP

//! \file numtools_root.hpp
//  \author Jacob Fields
//
//  \brief Declares some functions for root-finding.

#include <iostream>
#include <ps_types.hpp>
#include <cmath>
#include <limits>

namespace NumTools {

class Root {
  private:
  public:
    /// Maximum number of iterations
    unsigned int iterations;

    struct RootResult {
      bool success;
      Real tolerance;
      Real flast;
      Real err;
      unsigned int iterations;
      bool bracketed;
    };

    Root() : iterations(30) {}

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
    inline RootResult FalsePosition(Functor&& f, Real &lb, Real &ub, Real& x, Real tol,
                                    Types ... args) {
      RootResult result{true, tol, 0., 1e10, 0, true};
      int side = 0;
      Real ftest;
      unsigned int count = 0;
      // Get our initial bracket.
      Real flb = f(lb, args...);
      Real fub = f(ub, args...);
      Real xold, dist;
      x = lb;
      // If one of the bounds is already within tolerance of the root, we can skip all of this.
      if (std::fabs(flb)/lb <= tol) {
        x = lb;
        result.flast = flb;
        result.err = std::fabs(flb);
        return result;
      }
      else if (std::fabs(fub)/ub <= tol) {
        x = ub;
        result.flast = fub;
        result.err = std::fabs(fub);
        return result;
      }
      if (flb*fub > 0) {
        result.success = false;
        result.bracketed = false;
        return result;
      }
      Real eps = 0.0;
      //unsigned int iter_exp = std::log2(std::fabs(ub - lb)/(2.*tol));
      unsigned int iter_exp = 10000;
      do {
        xold = x;
        dist = ub - lb;
        // Calculate the new root position.
        if (count >= iter_exp) {
          // If we're converging too slowly, we revert to bisection.
          x = 0.5*(lb + ub);
        } else {
          x = (fub*lb - flb*ub)/(fub - flb);
          // Check that we're not stuck against one of the boundaries
          if ((ub - x)/dist < eps) {
            flb *= 0.5;
            x = (fub*lb - flb*ub)/(fub - flb);
          } else if ((x - lb)/dist < eps) {
            fub *= 0.5;
            x = (fub*lb - flb*ub)/(fub - flb);
          }
        }
        count++;
        // Calculate f at the prospective root.
        ftest = f(x,args...);
        result.err = std::fabs((x-xold)/x);
        if (result.err <= tol) {
          result.iterations = count;
          result.flast = ftest;
          return result;
        }
        // Check the sign of f. If f is on the same side as the lower bound, then we adjust
        // the lower bound. Similarly, if f is on the same side as the upper bound, we 
        // adjust the upper bound. If ftest falls on the same side twice, we weight one of
        // the sides to force the new root to fall on the other side. This allows us to
        // whittle down both sides at once and get better average convergence.
        if (ftest*flb >= 0) {
          if (side == 1) {
            Real m = 1. - ftest/flb;
            fub = (m > 0) ? fub*m : 0.5*fub;
            //fub /= 2.0;
          }
          flb = ftest;
          lb = x;
          side = 1;
        }
        else {
          if (side == -1) {
            Real m = 1. - ftest/fub;
            flb = (m > 0) ? flb*m : 0.5*flb;
            //flb /= 2.0;
          }
          fub = ftest;
          ub = x;
          side = -1;
        }
      }
      while (count < iterations);
      result.iterations = count;
      result.flast = ftest;
      result.err = fabs((x-xold)/x);

      // Return success if we're below the tolerance, otherwise report failure.
      result.success = result.err <= tol;
      return result;
    }
    
    // }}}

    // Hybrid {{{

    //! \brief Find the root of a functor f using a hybrid false position/bisection
    //
    // Find the root of a generic functor taking at least one argument. The first
    // argument is assumed to be the quantity of interest. All other arguments are
    // assumed to be constant parameters for the function. The root-finding method is
    // the ITP method, which is an enhancement of false position designed to guarantee
    // results no worse than bisection.
    //
    // \param[in]  f  The functor to find a root for. Its root function must take at
    //                least one argument.
    // \param[in,out]  lb  The lower bound for the root.
    // \param[in,out]  ub  The upper bound for the root.
    // \param[out]  x  The location of the root.
    // \param[in]  args  Additional arguments required by f.

    template<class Functor, class ... Types>
    inline RootResult ITP(Functor&& f, Real &lb, Real &ub, Real& x, Real tol,
                          Types ... args) {
      RootResult result{true, tol, 0., 1e10, 0, true};
      Real ftest;
      unsigned int count = 0;
      // Get our initial bracket.
      Real flb = f(lb, args...);
      Real fub = f(ub, args...);
      Real xold;
      x = lb;
      // If one of the bounds is already within tolerance of the root, we can skip all of this.
      if (std::fabs(flb) <= tol) {
        x = lb;
        result.flast = flb;
        result.err = std::fabs(flb);
        return result;
      }
      else if (std::fabs(fub) <= tol) {
        x = ub;
        result.flast = fub;
        result.err = std::fabs(fub);
        return result;
      }
      if (flb*fub > 0) {
        result.success = false;
        result.bracketed = false;
        return result;
      }
      // Hyperparameters
      Real k1 = 0.2*std::fabs(ub - lb);
      Real k2 = 2.0;
      Real n0 = 1.0;

      // Preprocessing
      Real nhalf = std::ceil(std::log2((ub - lb)/(2.0*tol)));
      int nmax = (int)(nhalf + n0);
      Real xhalf, xf, r, delta, dx, sigma, xt;
      do {
        // Updating parameters
        xold = x;
        xhalf = 0.5*(ub + lb);
        dx = ub - lb;
        r = tol*std::exp2(nmax - count) - 0.5*dx;
        delta = k1*std::pow(dx,k2);

        // Interpolation
        xf = (fub*lb - flb*ub)/(fub - flb);

        // Truncation
        Real diff = xhalf - xf;
        sigma = (0. < diff) - (diff < 0.);
        if (delta <= std::fabs(xhalf - xf)) {
          xt = xf + sigma*delta;
        } else {
          xt = xhalf;
        }

        // Projection
        if (std::fabs(xt - xhalf) <= r) {
          x = xt;
        } else {
          x = xhalf - sigma*r;
        }
        count++;

        // Terminate if we're within tolerance
        if (std::fabs((x - xold)/x) <= tol) {
          break;
        }

        // Updating interval
        ftest = f(x, args...);
        if (ftest*fub > 0) {
          ub = x;
          fub = ftest;
        } else {
          lb = x;
          flb = ftest;
        }
      }
      while (count < iterations);
      result.iterations = count;
      result.flast = ftest;
      result.err = fabs((x-xold)/x);

      // Return success if we're below the tolerance, otherwise report failure.
      result.success = result.err <= tol;
      return result;
    }
    
    // }}}


    // Brent {{{

    //! \brief Find the root of a functor f using Brent's method.
    //
    // Find the root of a generic functor taking at least one argument. The first
    // argument is assumed to be the quantity of interest. All other arguments are
    // assumed to be constant parameters for the function. The root-finding method
    // is Brent's method, which is a hybrid method that attempts to switch to provide
    // at least linear convergence.
    //
    // \param[in]  f  The functor to find a root for. Its root fucntion must take at
    //                least one argument.
    // \param[in,out]  lb  The lower bound for the root.
    // \param[in,out]  ub  The upper bound for the root.
    // \param[out]  x  The location of the root.
    // \param[in]  args  Additional arguments required by f.

    template<class Functor, class ... Types>
    inline RootResult Brent(Functor&& f, Real &lb, Real &ub, Real& x, Real tol,
                            Types ... args) {
      RootResult result{true, tol, 0., 1e10, 0, true};
      unsigned int count = 0;
      Real flb = f(lb, args...);
      Real fub = f(ub, args...);
      Real ftest;
      Real xold;
      x = lb;
      ftest = flb;
      // If one of the bounds is already within tolerance of the root, we can skip the
      // rest of the root solve.
      if (std::fabs(flb) <= tol) {
        x = lb;
        result.flast = flb;
        result.err = std::fabs(flb);
        return result;
      }
      else if (std::fabs(fub) <= tol) {
        x = ub;
        result.flast = fub;
        result.err = std::fabs(fub);
        return result;
      }
      if (flb*fub > 0) {
        result.success = false;
        result.bracketed = false;
        return result;
      }
      unsigned int iter_exp = std::log2(std::fabs(ub - lb)/(2.*tol));
      Real R, S, T, P, Q;
      Real c = lb, fc = flb;
      Real b = x, fb = flb;
      Real a = ub, fa = fub;
      do {
        xold = x;
        // Estimate the convergence rate
        Real conv = std::fabs((a - c)/(b - a));
        // Calculate the new root position
        if (count >= iter_exp || conv < 2.) {
          // If we're converging too slowly, we revert to bisection.
          x = 0.5*(lb + ub);
        } else {
          // Find the new root using interpolation.
          if (fb == fa) {
            // Use the secant method to start because we have colliding points.
            x = b + fb*(b - c)/(fb - fc);
          }
          else if (fb == fc) {
            // Use the secant method to start because we have colliding points.
            x = b + fb*(b - a)/(fb - fa);
          }
          else {
            // Find the new root with inverse quadratic interpolation.
            R = fb/fc;
            S = fb/fa;
            T = fa/fb;
            P = S*(T*(R - T)*(c - b) - (1. - R)*(b - a));
            Q = (T - 1.)*(R - 1.)*(S - 1.);
            x = b + P/Q;
          }
          // If x is outside the bounds, use bisection instead.
          if (x > ub || x < lb || !std::isfinite(x)) {
            x = 0.5*(lb + ub);
          }
          else {
            std::cout << "Didn't need bisection!\n";
          }
        }
        count++;
        // Calculate f at the prospective root.
        ftest = f(x,args...);
        result.err = std::fabs((x-b)/x);
        if (result.err <= tol) {
          result.iterations = count;
          result.flast = ftest;
          return result;
        }
        // Update the bounds
        if (ftest*flb >= 0.0) {
          flb = ftest;
          lb = x;
        } else {
          fub = ftest;
          ub = x;
        }
        // Cycle the points.
        c = a;
        a = b;
        b = x;
        fc = fa;
        fa = fb;
        fb = ftest;
      } while (count < iterations);
      result.iterations = count;
      result.flast = ftest;
      result.err = std::fabs((x-xold)/x);

      // Return success if we're below the tolerance, otherwise report failure.
      result.success = result.err <= tol;
      return result;
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
    inline bool Chandrupatla(Functor&& f, Real &lb, Real &ub, Real& x, Real tol,
                             Types ... args) {
      unsigned int count = 0;
      //last_count = 0;
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
      //last_count = count;

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
    inline bool NewtonSafe(Functor&& f, Real &lb, Real &ub, Real& x, Real tol,
                           Types ... args) {
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
        // Correct the bounds.
        if (fx*flb > 0) {
          flb = fx;
          lb = xold;
        }
        else if (fx*fub > 0) {
          fub = fx;
          ub = xold;
        }
        x = x - fx/dfx;
        // Check that the root is bounded properly.
        if (x > ub || x < lb) {
          // Revert to bisection if the root is not converging.
          x = 0.5*(ub + lb);
          //f(fx, dfx, x, args...);
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
