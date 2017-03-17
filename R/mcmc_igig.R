### mcmc_igig.R
### Copyright (c) 2017 Victhor S. Sartório
### This file and its contents are licensed under the terms of the MIT License
###
### This file implements R wrappers for the C functions relating to the file
### `mcmc_igig.R`.


# @TODO: Add documentation
# @TODO: Make integer coercion inside the C code
#' MCMC Routine for Inference with Inverse Gamma priors on \eqn{\sigma_x^2} and
#' \eqn{\sigma_m^2}
#'
#' @useDynLib iarrec
#' @export
mcmc_igig <- function(numit, t1_start, t2_start, sigma1, sigma2, Y, a1, b1, a2, b2)
	.Call('mcmc_igig',
		  as.integer(numit), t1_start, t2_start, sigma1, sigma2, Y, a1, b1, a2, b2)
