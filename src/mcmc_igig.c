/**
 * mcmc_igig.c
 * Copyright (c) 2017 Victhor S. Sartório
 * This file and its contents are licensed under the terms of the MIT License
 *
 * Implements a MCMC metropolis-hastings routine for inference over the model
 * 		Y ~ N(0, t1 + t2)
 * using IG priors for t1 and t2.
 */

#include "common.h"
#include "gauss.h"
#include "linvgamma.h"
#include "dinvgamma.h"

#include <Rmath.h>

// @TODO: Add support for user-defined mu instead of 0 (make sure to change logposterior)


static double
logposterior(double t1, double t2,       // Points to be evaluated
             const double *restrict Y,   // Vector of observations
             R_xlen_t n,                 // Number of observations
             double a1, double b1,       // Prior of t1
             double a2, double b2)       // Prior of t2
{
	// Calculate value of the prior log-densities
	double prior_t1 = c_linvgamma_s(t1, a1, b1);
	double prior_t2 = c_linvgamma_s(t2, a2, b2);

	// Calculate the log-likelihood of Y given t1 and t2
	long double loglikelihood = 0.0;
	const double s2 = t1+t2;
	for (R_xlen_t i = 0; i < n; i++)
		loglikelihood += c_lgauss_s(Y[i], 0, s2);

	// Return the log-posterior \prorp P(Y|t1,t2) * P(t1) * P(t2)
	return prior_t1 + prior_t2 + loglikelihood;
}


SEXP mcmc_igig(SEXP r_numit,                            // Number of iterations
               SEXP r_t1_start, SEXP r_t2_start,        // Starting points for t1 and t2
               SEXP r_sigma1, SEXP r_sigma2,            // The parameters of the chain
               SEXP r_obs,                              // The observation vector
               SEXP r_a1, SEXP r_b1,                    // The prior for t1
               SEXP r_a2, SEXP r_b2)                    // The prior for t2
{
	// Syntatic check: Check if types and length are correct
	if (TYPEOF(r_numit) != INTSXP || XLENGTH(r_numit) != 1)
		Rf_error("Expected numit to be an atomic integer.");
	if (TYPEOF(r_t1_start) != REALSXP || XLENGTH(r_t1_start) != 1)
		Rf_error("Expected t1_start to be an atomic real.");
	if (TYPEOF(r_t2_start) != REALSXP || XLENGTH(r_t2_start) != 1)
		Rf_error("Expected t2_start to be an atomic real.");
	if (TYPEOF(r_sigma1) != REALSXP || XLENGTH(r_sigma1) != 1)
		Rf_error("Expected r_sigma1 to be an atomic real.");
	if (TYPEOF(r_sigma2) != REALSXP || XLENGTH(r_sigma2) != 1)
		Rf_error("Expected r_sigma2 to be an atomic real.");
	if (TYPEOF(r_obs) != REALSXP || XLENGTH(r_obs) <= 1)
		Rf_error("Expected Y to be a vector with real observations.");
	if (TYPEOF(r_a1) != REALSXP || XLENGTH(r_a1) != 1)
		Rf_error("Expected a1 to be an atomic real.");
	if (TYPEOF(r_b1) != REALSXP || XLENGTH(r_b1) != 1)
		Rf_error("Expected b1 to be an atomic real.");
	if (TYPEOF(r_a2) != REALSXP || XLENGTH(r_a2) != 1)
		Rf_error("Expected a2 to be an atomic real.");
	if (TYPEOF(r_b2) != REALSXP || XLENGTH(r_b2) != 1)
		Rf_error("Expected b2 to be an atomic real.");

	// Semantic check: Check if values are valid when applicable
	if (INTEGER(r_numit)[0] <= 1)
		Rf_error("Expected numit to be higher than 1.");
	if (REAL(r_t1_start)[0] <= 0)
		Rf_error("Expected t1_start to be positive.");
	if (REAL(r_t2_start)[0] <= 0)
		Rf_error("Expected t2_start to be positive.");
	if (REAL(r_sigma1)[0] <= 0)
		Rf_error("Expected sigma1 to be positive.");
	if (REAL(r_sigma2)[0] <= 0)
		Rf_error("Expected sigma2 to be positive.");
	if (REAL(r_a1)[0] <= 0)
		Rf_error("Expected a1 to be positive.");
	if (REAL(r_b1)[0] <= 0)
		Rf_error("Expected b1 to be positive.");
	if (REAL(r_a2)[0] <= 0)
		Rf_error("Expected a2 to be positive.");
	if (REAL(r_b2)[0] <= 0)
		Rf_error("Expected b2 to be positive.");

	// Add alias for numit
	const int numit = INTEGER(r_numit)[0];

	// Build the list to be returned
	SEXP r_t1 = PROTECT(allocVector(REALSXP, numit));
	SEXP r_t2 = PROTECT(allocVector(REALSXP, numit));
	SEXP r_ar = PROTECT(allocVector(REALSXP, 1));
	SEXP ret = PROTECT(allocVector(VECSXP, 3));
	SET_VECTOR_ELT(ret, 0, r_t1);
	SET_VECTOR_ELT(ret, 1, r_t2);
	SET_VECTOR_ELT(ret, 2, r_ar);

	// Add a class to it
	SEXP class_name = PROTECT(allocVector(STRSXP, 1));
	SET_STRING_ELT(class_name, 0, mkChar("MCMCresult2"));
	setAttrib(ret, R_ClassSymbol, class_name);

	// Add names to the list
	SEXP names = PROTECT(allocVector(STRSXP, 3));
	SET_STRING_ELT(names, 0, mkChar("theta_1"));
	SET_STRING_ELT(names, 1, mkChar("theta_2"));
	SET_STRING_ELT(names, 2, mkChar("ar"));
	setAttrib(ret, R_NamesSymbol, names);

	// Set aliases for values
	double *t1 = REAL(r_t1);
	double *t2 = REAL(r_t2);

	// Set starting value
	t1[0] = REAL(r_t1_start)[0];
	t2[0] = REAL(r_t2_start)[0];

	// Set aliases for the parameters
	const double sigma1 = REAL(r_sigma1)[0];
	const double sigma2 = REAL(r_sigma2)[0];

	// Variable to count number of accepted steps
	R_xlen_t accepted = 0;

	// Set alias for the observations and the hyperparameters
	double *Y = REAL(r_obs);
	const R_xlen_t n = XLENGTH(r_obs);
	const double a1 = REAL(r_a1)[0];
	const double b1 = REAL(r_b1)[0];
	const double a2 = REAL(r_a2)[0];
	const double b2 = REAL(r_b2)[0];

	// Iterative process
	for (R_xlen_t i = 1; i < numit; i++) {
		// Get proposal values
		const double p1 = Rf_rnorm(t1[i-1], sigma1);
		const double p2 = Rf_rnorm(t2[i-1], sigma2);

		// If one of them is non-positive, reject straight away
		if (p1 <= 0 || p2 <= 0) {
			t1[i] = t1[i-1];
			t2[i] = t2[i-1];
			continue;
		}

		// Calculate the `alpha` of acceptance
		const double logpost_prp = logposterior(p1, p2, Y, n, a1, b1, a2, b2);
		const double logpost_old = logposterior(t1[i-1], t2[i-1], Y, n, a1, b1, a2, b2);
		const double alpha = exp(logpost_prp - logpost_old);

		// If `alpha` > 1 accept. If `alpha` in [0,1] accept with probability `alpha`
		if (Rf_runif(0, 1) < alpha) {
			t1[i] = p1;
			t2[i] = p2;
			accepted++;
		} else {
			t1[i] = t1[i-1];
			t2[i] = t2[i-1];
		}
	}

	// Calculate acceptance ratio
	REAL(r_ar)[0] = ((double) accepted) / numit;

	UNPROTECT(6);
	return ret;
}
