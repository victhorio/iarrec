### slr.R
### Copyright (c) 2017 Victhor S. Sartório
### This file and its contents are licensed under the terms of the MIT License
###
### This file implements R wrappers for the C functions relating to the file
### `slr_igig.c`.


# @TODO: Improve documentation
# @TODO: Make integer coercion of `numit` inside the C code
#' Simple Linear Regression with Inverse Gamma Priors on the Non-Identifiable Parameters
#'
#' Runs a MCMC to generate a sample for the posterior of the parameters of a simple
#' linear regression model with errors in the covariates.
#'
#' @param numit Number of iterations
#' @param start Vector with starting values for \eqn{\beta_0}, \eqn{\beta_1}, \eqn{\mu_x},
#' \eqn{\sigma_x^2}, \eqn{\sigma_r^2} and \eqn{\sigma_m^2} respectively
#' @param sigma The vector with the random walk parameter for each of the variables
#' @param obs_Y The vector of the Y observations
#' @param obs_X The vector of the X observations
#' @param beta0 Mean and variance of the normal prior for \eqn{\beta_0}
#' @param beta1 Mean and variance of the normal prior for \eqn{\beta_1}
#' @param mux Mean and variance of the normal prior for \eqn{\mu_x}
#' @param sx Parameters of the inverse gamma prior for \eqn{\sigma_x^2}
#' @param sr Parameters of the inverse gamma prior for \eqn{\sigma_r^2}
#' @param sm Parameters of the inverse gamma prior for \eqn{\sigma_m^2}
#' @return An 'MCMCresultSLR' object.
#'
#' @examples
#' x <- rnorm(60, 50, sqrt(20)/2)
#' xs <- x + rnorm(60, 0, sqrt(20)/2)
#' y <- 50 + 50*x + rnorm(60, 0, sqrt(50))
#' O <- slr_igig(75000, rep(1, 6), rep(0.6, 6), y, xs, c(0, 500), c(0, 500), c(0, 500), c(3, 20), c(2, 50), c(3, 20))
#'
#' @useDynLib iarrec
#' @export
slr_igig <- function(numit, start, sigma, obs_Y, obs_X, beta0, beta1, mux, sx, st, sm)
	.Call('slr_igig',
		  as.integer(numit), start, sigma, obs_Y, obs_X, beta0, beta1, mux, sx, st, sm)
