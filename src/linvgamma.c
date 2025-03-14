/**
 * linvgamma.c
 * Copyright (c) 2017 Victhor S. Sartório
 * This file and its contents are licensed under the terms of the MIT License
 *
 * Implements the log-density function for the inverse gamma distribution, under the R
 * interfaceable wrapper `linvgamma` and under the C implementation `c_linvgamma`.
 * Implementation of the R function and docstring definition can be found in the R
 * folder in the file `invgamma.R`.
 */

#include "linvgamma.h"


/**
 * c_linvgamma_s
 *
 * @param x The point at which the log-density will be evaluated
 * @param a The scale parameter
 * @param b The shape parameter
 *
 * Note that the function relies that all arguments are valid.
 */
HOTFUNC double
c_linvgamma_s(double x, double a, double b)
{
	return a * log(b) - lgamma(a) - (a+1) * log(x) - b/x;
}


/**
 * c_linvgamma
 *
 * @param x The vector containing the points at which the density will be evaluated
 * @param n The number of elements of the vector x
 * @param a The scale parameter
 * @param b The shape parameter
 * @param r The vector of length (at least) n where the results will be stored
 *
 * Note that the vectors `x` and `r` cannot be aliased.
 * Note that the function relies that all arguments are valid.
 */
HOTFUNC void
c_linvgamma(double const *restrict x, R_xlen_t n, double a, double b, double *restrict r)
{
	for (R_xlen_t i = 0; i < n; i++)
		r[i] = a * log(b) - lgamma(a) - (a+1) * log(x[i]) - b/x[i];
}


/**
 * dinvgamma
 *
 * R interfaceable wrapper for the `c_dinvgamma` implementation.
 * Further documentation on the R wrapper found in `R/invgamma.R`.
 */
SEXP
linvgamma(SEXP rx, SEXP ra, SEXP rb)
{
	int cond_protect = 0; // Number of conditional PROTECT calls

	R_xlen_t n = Rf_xlength(rx);
	if (n <= 0) {
		Rf_error("First argument is empty.");
	}

	// If an INTSXP vector was received, coerce a REALSXP one into cx
	// If an REALSXP vector was received, alias cx to it
	// Otherwise raise an error
	SEXP cx;
	if (TYPEOF(rx) == INTSXP) {
		cx = PROTECT(Rf_coerceVector(rx, REALSXP));
		cond_protect++;
	} else if (TYPEOF(rx) == REALSXP) {
		cx = rx;
	} else {
		Rf_error("Firt argument is not numeric.");
	}

	// Make sure second and third arguments are REALSXP, atomic and positive
	if (Rf_xlength(ra) != 1 || TYPEOF(ra) != REALSXP || REAL(ra)[0] <= 0)
		Rf_error("Expected second argument to be real valued, atomic and positive.");
	if (Rf_xlength(rb) != 1 || TYPEOF(rb) != REALSXP || REAL(rb)[0] <= 0)
		Rf_error("Expected third argument to be real valued, atomic and positive.");

	SEXP ret = PROTECT(Rf_allocVector(REALSXP, n));

	double const *const x = REAL(cx);
	double *const r = REAL(ret);
	double const a = REAL(ra)[0];
	double const b = REAL(rb)[0];

	c_linvgamma(x, n, a, b, r);

	UNPROTECT(1 + cond_protect);
	return ret;
}
