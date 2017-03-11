/**
 * gauss.c
 * Copyright (c) 2017 Victhor S. Sartório
 * This file and its contents are licensed under the terms of the MIT License
 *
 * Implements the log-density function for the gaussian distribution.
 */

#include "gauss.h"


/**
* c_lgauss_s
*
* @param x The point at which the log-density will be evaluated
* @param mu The mean of the distribution
* @param s2 The variance of the distribution
*
* Note that the function relies that all arguments are valid.
*/
HOTFUNC double
c_lgauss_s(double x, double mu, double s2)
{
	return -0.5 * (log(2*PI*s2) + (x-mu)*(x-mu)/s2);
}


/**
* c_lgauss
*
* @param x The vector containing the points at which the log-density will be evaluated
* @param a The scale parameter
* @param b The shape parameter
* @param r The vector of length (at least) n where the results will be stored
*
* Note that the vectors `x` and `r` cannot be aliased.
* Note that the function relies that all arguments are valid.
*/
HOTFUNC void
c_lgauss(double const *restrict x, R_xlen_t n, double mu, double s2, double *restrict r)
{
	for (R_xlen_t i = 0; i < n; i++)
		r[i] = -0.5 * (log(2*PI*s2) + (x[i]-mu)*(x[i]-mu)/s2);
}


/**
* lgauss
*
* R interfaceable wrapper for the `c_lgauss` implementation.
* Further documentation on the R wrapper found in `R/gauss.R`.
*/
SEXP
lgauss(SEXP rx, SEXP rmu, SEXP rs2)
{
	int cond_protect = 0; // Number of conditional PROTECT calls

	R_xlen_t n = Rf_xlength(rx);
	if (n <= 0)
		Rf_error("First argument is empty.");

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

	// Make sure third argument is REALSXP, atomic and positive
	if (Rf_xlength(rs2) != 1 || TYPEOF(rs2) != REALSXP || REAL(rs2)[0] <= 0)
		Rf_error("Expected third argument to be real valued, atomic and positive.");

	SEXP ret = PROTECT(Rf_allocVector(REALSXP, n));

	double const *const x = REAL(cx);
	double *const r = REAL(ret);
	double const mu = REAL(rmu)[0];
	double const s2 = REAL(rs2)[0];

	c_lgauss(x, n, mu, s2, r);

	UNPROTECT(1 + cond_protect);
	return ret;
}
