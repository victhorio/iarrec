/**
 * mcmc_igig.c
 * Copyright (c) 2017 Victhor S. Sartório
 * This file and its contents are licensed under the terms of the MIT License
 *
 * Implements a MCMC metropolis-hastings routine for inference over the model
 *     Y   = \beta_0 + \beta_1 X + \epsilon_r
 *     X^* = X + \epsilon_m
 * using IG priors for \sigma_x and \sigma_m.
 */

// @TODO: Improve that ^ docstring

#include "common.h"
#include "gauss.h"
#include "bigauss.h"
#include "linvgamma.h"
#include "dinvgamma.h"

#include <Rmath.h>


static double
logposterior(double b0, double b1, double mx, double sx, double sr, double sm, // Points
             double mu0, double s0,    // Parameters of b0 prior
             double mu1, double s1,    // Parameters of b1 prior
             double mumu, double smu,  // Parameters of mu prior
             double ax, double bx,     // Parameters of sx prior
             double ar, double br,     // Parameters of sr prior
             double am, double bm,     // Parameters of sm prior
             const double *restrict Y, // Y observations
             const double *restrict X, // X observations
             const R_xlen_t n)         // Number of observations
{
	// Calculate value of the prior log-densities
	double prior_b0 = c_lgauss_s(b0, mu0, s0);
	double prior_b1 = c_lgauss_s(b1, mu1, s1);
	double prior_mx = c_lgauss_s(mx, mumu, smu);
	double prior_sx = c_linvgamma_s(sx, ax, bx);
	double prior_sr = c_linvgamma_s(sr, ar, br);
	double prior_sm = c_linvgamma_s(sm, am, bm);

	// Calculate the log-likelihood of Y given the parameters
	const double ll_muy = b0 + b1 * mx;
	const double ll_mux = mx;
	const double ll_syy = b1*b1*sx+sr;
	const double ll_sxx = sx+sm;
	const double ll_cov = b1*sx;
	long double ll = 0.0;
	for (R_xlen_t i = 0; i < n; i++)
		ll += c_lbigauss_s(Y[i], X[i], ll_muy, ll_mux, ll_syy, ll_sxx, ll_cov);

	// Return the log-posterior \prorp P(Y|theta) * P(theta)
	return prior_b0 + prior_b1 + prior_mx + prior_sx + prior_sr + prior_sm + ll;
}


SEXP
slr_igig(SEXP r_numit,     // Number of iterations
         SEXP r_start,     // Starting points of the chains
         SEXP r_param,     // The parameters of the chain
         SEXP r_Y,         // The observation vector Y
         SEXP r_X,         // The observation vector X
         SEXP r_beta0,     // The beta0 prior parameters
         SEXP r_beta1,     // The beta1 prior parameters
         SEXP r_mux,       // The mux prior parameters
         SEXP r_sx,        // The sigmax prior parameters
         SEXP r_sr,        // The sigmar prior parameters
         SEXP r_sm)        // The sigmam prior parameters
{
	// Syntatic check: Check if types and length of arguments are correct
	if (TYPEOF(r_numit) != INTSXP || XLENGTH(r_numit) != 1)
		Rf_error("Expected numit to be an atomic integer.");
	if (TYPEOF(r_start) != REALSXP || XLENGTH(r_start) != 6)
		Rf_error("Expected start to be a real vector with 6 entries.");
	if (TYPEOF(r_param) != REALSXP || XLENGTH(r_param) != 6)
		Rf_error("Expected param to be a real vector with 6 entries.");
	if (TYPEOF(r_Y) != REALSXP || XLENGTH(r_Y) < 1)
		Rf_error("Expected Y to be a non-empty real vector.");
	if (TYPEOF(r_X) != REALSXP || XLENGTH(r_X) != XLENGTH(r_Y))
		Rf_error("Expected X to be a real vector with conforming length with Y.");
	if (TYPEOF(r_beta0) != REALSXP || XLENGTH(r_beta0) != 2)
		Rf_error("Expected beta0 to be a real valued pair of values.");
	if (TYPEOF(r_beta1) != REALSXP || XLENGTH(r_beta1) != 2)
		Rf_error("Expected beta1 to be a real valued pair of values.");
	if (TYPEOF(r_mux) != REALSXP || XLENGTH(r_mux) != 2)
		Rf_error("Expected mux to be a real valued pair of values.");
	if (TYPEOF(r_sx) != REALSXP || XLENGTH(r_sx) != 2)
		Rf_error("Expected sx to be a real valued pair of values.");
	if (TYPEOF(r_sr) != REALSXP || XLENGTH(r_sr) != 2)
		Rf_error("Expected sr to be a real valued pair of values.");
	if (TYPEOF(r_sm) != REALSXP || XLENGTH(r_sm) != 2)
		Rf_error("Expected sm to be a real valued pair of values.");

	// Extract values
	const int numit = INTEGER(r_numit)[0];
	const double *const start = REAL(r_start);
	const double start_b0 = start[0];
	const double start_b1 = start[1];
	const double start_mx = start[2];
	const double start_sx = start[3];
	const double start_sr = start[4];
	const double start_sm = start[5];
	const double *const jmp = REAL(r_param);
	const double jmp_b0 = jmp[0];
	const double jmp_b1 = jmp[1];
	const double jmp_mx = jmp[2];
	const double jmp_sx = jmp[3];
	const double jmp_sr = jmp[4];
	const double jmp_sm = jmp[5];
	const double *const Y = REAL(r_Y);
	const double *const X = REAL(r_X);
	const double *const prior_b0 = REAL(r_beta0);
	const double mu0 = prior_b0[0];
	const double s0 = prior_b0[1];
	const double *const prior_b1 = REAL(r_beta1);
	const double mu1 = prior_b1[0];
	const double s1 = prior_b1[1];
	const double *const prior_mx = REAL(r_mux);
	const double mumu = prior_mx[0];
	const double smu = prior_mx[1];
	const double *const prior_sx = REAL(r_sx);
	const double ax = prior_sx[0];
	const double bx = prior_sx[1];
	const double *const prior_sr = REAL(r_sr);
	const double ar = prior_sr[0];
	const double br = prior_sr[1];
	const double *const prior_sm = REAL(r_sm);
	const double am = prior_sm[0];
	const double bm = prior_sm[1];

	// Semantic check: Check if values are valid
	if (numit <= 1)
		Rf_error("Expected numit to be higher than 1.");
	if (start_sx <= 0 || start_sr <= 0 || start_sm <= 0)
		Rf_error("Starting values for variances need to be positive.");
	if (jmp_b0 <= 0 || jmp_b1 <= 0 || jmp_mx <= 0 || jmp_sx <= 0 ||
        jmp_sr <= 0 || jmp_sm <= 0)
		Rf_error("Jump parameters need to be positive.");
	if (s0 <= 0 || s1 <= 0 || smu <= 0)
		Rf_error("Variance parameter of normal priors need to be positive.");
	if (ax <= 0 || bx <= 0 || ar <= 0 || br <= 0 || am <= 0 || bm <= 0)
		Rf_error("Parameters of the inverse gamma priors need to be positive.");

	// Add alias for the number of observations
	const R_xlen_t n = XLENGTH(r_Y);

	// Build the list to be returned
	SEXP ret_b0 = PROTECT(allocVector(REALSXP, numit));
	SEXP ret_b1 = PROTECT(allocVector(REALSXP, numit));
	SEXP ret_mx = PROTECT(allocVector(REALSXP, numit));
	SEXP ret_sx = PROTECT(allocVector(REALSXP, numit));
	SEXP ret_sr = PROTECT(allocVector(REALSXP, numit));
	SEXP ret_sm = PROTECT(allocVector(REALSXP, numit));
	SEXP ret_AR = PROTECT(allocVector(REALSXP, 1));
	SEXP ret = PROTECT(allocVector(VECSXP, 7));
	SET_VECTOR_ELT(ret, 0, ret_b0);
	SET_VECTOR_ELT(ret, 1, ret_b1);
	SET_VECTOR_ELT(ret, 2, ret_mx);
	SET_VECTOR_ELT(ret, 3, ret_sx);
	SET_VECTOR_ELT(ret, 4, ret_sr);
	SET_VECTOR_ELT(ret, 5, ret_sm);
	SET_VECTOR_ELT(ret, 6, ret_AR);

	// Add a class to it
	SEXP class_name = PROTECT(allocVector(STRSXP, 1));
	SET_STRING_ELT(class_name, 0, mkChar("MCMCresultSLR"));
	setAttrib(ret, R_ClassSymbol, class_name);

	// Add names to the list
	SEXP names = PROTECT(allocVector(STRSXP, 7));
	SET_STRING_ELT(names, 0, mkChar("beta_0"));
	SET_STRING_ELT(names, 1, mkChar("beta_1"));
	SET_STRING_ELT(names, 2, mkChar("mu_x"));
	SET_STRING_ELT(names, 3, mkChar("sigma2_x"));
	SET_STRING_ELT(names, 4, mkChar("sigma2_r"));
	SET_STRING_ELT(names, 5, mkChar("sigma2_m"));
	SET_STRING_ELT(names, 6, mkChar("ar"));
	setAttrib(ret, R_NamesSymbol, names);

	// Set aliases for the chains
	double *b0 = REAL(ret_b0);
	double *b1 = REAL(ret_b1);
	double *mx = REAL(ret_mx);
	double *sx = REAL(ret_sx);
	double *sr = REAL(ret_sr);
	double *sm = REAL(ret_sm);

	// Set starting values
	b0[0] = start_b0;
	b1[0] = start_b1;
	mx[0] = start_mx;
	sx[0] = start_sx;
	sr[0] = start_sr;
	sm[0] = start_sm;

	// Variable to count number of accepted steps
	R_xlen_t accepted = 0;

	double old_logpost = logposterior(b0[0], b1[0], mx[0], sx[0], sr[0], sm[0],
	                                  mu0, s0, mu1, s1, mumu, smu, ax, bx, ar, br, am, bm,
	                                  Y, X, n);
	double prp_logpost;

	// Iterative process
	for (R_xlen_t i = 1; i < numit; i++) {
		// Get proposal values
		const double p_b0 = Rf_rnorm(b0[i-1], jmp_b0);
		const double p_b1 = Rf_rnorm(b1[i-1], jmp_b1);
		const double p_mx = Rf_rnorm(mx[i-1], jmp_mx);
		const double p_sx = Rf_rnorm(sx[i-1], jmp_sx);
		const double p_sr = Rf_rnorm(sr[i-1], jmp_sr);
		const double p_sm = Rf_rnorm(sm[i-1], jmp_sm);

		// If an impossible value was proposed, reject it
		if (p_sx <= 0 || p_sr <= 0 || p_sm <= 0) {
			b0[i] = b0[i-1];
			b1[i] = b1[i-1];
			mx[i] = mx[i-1];
			sx[i] = sx[i-1];
			sr[i] = sr[i-1];
			sm[i] = sm[i-1];
			continue;
		}

		// Calculate the `alpha` of acceptance
		prp_logpost = logposterior(p_b0, p_b1, p_mx, p_sx, p_sr, p_sm,
		                                        mu0, s0, mu1, s1, mumu, smu,
		                                        ax, bx, ar, br, am, bm,
		                                        Y, X, n);
		const double alpha = exp(prp_logpost - old_logpost);

		// If `alpha` > 1 accept. If `alpha` in [0,1] accept with probability `alpha`
		if (Rf_runif(0, 1) < alpha) {
			b0[i] = p_b0;
			b1[i] = p_b1;
			mx[i] = p_mx;
			sx[i] = p_sx;
			sr[i] = p_sr;
			sm[i] = p_sm;
			old_logpost = prp_logpost;
			accepted++;
		} else {
			b0[i] = b0[i-1];
			b1[i] = b1[i-1];
			mx[i] = mx[i-1];
			sx[i] = sx[i-1];
			sr[i] = sr[i-1];
			sm[i] = sm[i-1];
		}
	}

	// Calculate acceptance ratio
	REAL(ret_AR)[0] = ((double) accepted) / numit;

	UNPROTECT(10);
	return ret;
}
