/**
 * bigauss.c
 * Copyright (c) 2017 Victhor S. Sartório
 * This file and its contents are licensed under the terms of the MIT License
 *
 * Implements the log-density function for the bivariate gaussian distribution.
 */

#include "bigauss.h"

HOTFUNC double
sqr(double a)
{
	return a*a;
}


/**
 * c_lbigauss_s
 *
 * @param x The value of the x variable
 * @param y The value of the y variable
 * @param mx The mean of the x variable
 * @param my The mean of the y variable
 * @param sxx The variance of the distribution of x
 * @param syy The variance of the distribution of y
 * @param sxy The covariance between x and y
 *
 * Note that the function relies that all arguments are valid.
 */
HOTFUNC double
c_lbigauss_s(double x, double y, double mx, double my, double sxx, double syy, double sxy)
{
	const double sx = sqrt(sxx);
	const double sy = sqrt(syy);
	const double r = sxy / (sx * sy);
	const double r2 = sqr(r);
	const double z = sqr(x-mx)/sxx - 2*r*(x-mx)*(y-my)/(sx*sy) + sqr(y-my)/syy;
	return - log(2*PI) - 0.5*log(sxx*syy*(1-r2)) - z / (2 * (1-r2));
}


/**
 * lbigauss
 *
 * R interfaceable wrapper for the `c_lbigauss_s` implementation.
 * Further documentation on the R wrapper found in `R/gauss.R`.
 */
SEXP
lbigauss(SEXP rx, SEXP ry, SEXP rmx, SEXP rmy, SEXP rsxx, SEXP rsyy, SEXP rsxy)
{
	// Syntatic check: Check if types and length are correct
	if (TYPEOF(rx) != REALSXP || XLENGTH(rx) != 1)
		Rf_error("Expected x to be an atomic real.");
	if (TYPEOF(ry) != REALSXP || XLENGTH(ry) != 1)
		Rf_error("Expected y to be an atomic real.");
	if (TYPEOF(rmx) != REALSXP || XLENGTH(rmx) != 1)
		Rf_error("Expected mx to be an atomic real.");
	if (TYPEOF(rmy) != REALSXP || XLENGTH(rmy) != 1)
		Rf_error("Expected my to be an atomic real.");
	if (TYPEOF(rsxx) != REALSXP || XLENGTH(rsxx) != 1)
		Rf_error("Expected sxx to be an atomic real.");
	if (TYPEOF(rsyy) != REALSXP || XLENGTH(rsyy) != 1)
		Rf_error("Expected syy to be an atomic real.");
	if (TYPEOF(rsxy) != REALSXP || XLENGTH(rsxy) != 1)
		Rf_error("Expected sxy to be an atomic real.");

	// Extract values
	const double x = REAL(rx)[0];
	const double y = REAL(ry)[0];
	const double mx = REAL(rmx)[0];
	const double my = REAL(rmy)[0];
	const double sxx = REAL(rsxx)[0];
	const double syy = REAL(rsyy)[0];
	const double sxy = REAL(rsxy)[0];

	// Semantic check: Check if values are valid
	if (sxx <= 0)
		Rf_error("Expected sxx to be positive.");
	if (syy <= 0)
		Rf_error("Expected syy to be positive.");

	SEXP ret = PROTECT(allocVector(REALSXP, 1));
	REAL(ret)[0] = c_lbigauss_s(x, y, mx, my, sxx, syy, sxy);

	UNPROTECT(1);
	return ret;
}