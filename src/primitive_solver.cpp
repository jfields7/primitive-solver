#include <primitive_solver.hpp>
#include <cmath>

using namespace PrimitiveSolver;

void PrimitiveSolver::f_upper(Real &f, Real &df, Real mu, 
                              Real h0, Real rsq, Real bsq, Real rbsq){
  Real x = xf(mu, bsq);
  Real sdet = sqrt(h0*h0 + rbarsq(mu, rsq, rbsq, x));
  f = mu*sdet;
  df = sdet + mu*drbarsq(mu, rsq, bsq, rbsq, x)/(2.0*sdet);
}

inline Real PrimitiveSolver::rbarsq(Real mu, Real rsq, Real rbsq, Real x){
  return rsq*x*x + mu*x*(1.0 + x)*rbsq;
}

inline Real PrimitiveSolver::drbarsq(Real mu, Real rsq, Real bsq, Real rbsq, Real x){
  Real dx = dxf(x, bsq);
  return rbsq*x*(1.0 + x) + (mu*rbsq + 2.0*(mu*rbsq + rsq)*x);
}

inline Real PrimitiveSolver::xf(Real mu, Real bsq){
  return 1.0/(1.0 + mu*bsq);
}

inline Real PrimitiveSolver::dxf(Real x, Real bsq){
  return bsq*x*x;
}
