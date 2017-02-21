#ifndef R_IDENTIFIABILITY_INVGAMMA_HEADER
#define R_IDENTIFIABILITY_INVGAMMA_HEADER

#include <R.h>
#include <Rinternals.h>

void
c_dinvgamma(double const *restrict x, R_xlen_t n, double a, double b, double *restrict r);

#endif
