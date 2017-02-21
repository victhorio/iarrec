#include "linvgamma.h"


// C interfaceable function
// Relies on all arguments being valid
void __attribute__((hot))
c_linvgamma(double const *restrict x, R_xlen_t n, double a, double b, double *restrict r)
{
    for (R_xlen_t i = 0; i < n; i++)
        r[i] = a * log(b) - lgamma(a) - (a+1) * log(x[i]) - b/x[i];
}

// R interfaceable wrapper
SEXP
linvgamma(SEXP rx, SEXP ra, SEXP rb)
{
    int np = 0; // Number of PROTECT calls

    // Make sure rx is non-empty
    R_xlen_t n = Rf_xlength(rx);
    if (n <= 0) {
        Rf_error("First argument is empty.");
    }

    // Make sure cx is a REALSXP vector
    SEXP cx;
    if (TYPEOF(rx) == INTSXP) {
        cx = PROTECT(Rf_coerceVector(rx, REALSXP));
        np++;
    } else if (TYPEOF(rx) == REALSXP) {
        cx = rx;
    } else {
        Rf_error("Firt argument is not numeric.");
    }

    // Make sure second and third arguments are REALSXP, atomic and positive
    if (Rf_xlength(ra) != 1 || TYPEOF(ra) != REALSXP || REAL(ra)[0] <= 0)
        Rf_error("Expected second argument to be real valued, atomic and strictly positive.");
    if (Rf_xlength(rb) != 1 || TYPEOF(rb) != REALSXP || REAL(ra)[0] <= 0)
        Rf_error("Expected third argument to be real valued, atomic and strictly positive.");

    SEXP ret = PROTECT(Rf_allocVector(REALSXP, n));
    np++;

    double const *const x = REAL(cx);
    double *const r = REAL(ret);
    double const a = REAL(ra)[0];
    double const b = REAL(rb)[0];

    c_linvgamma(x, n, a, b, r);

    UNPROTECT(np);
    return ret;
}
