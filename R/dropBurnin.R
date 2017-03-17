### dropBurnin.R
### Copyright (c) 2017 Victhor S. Sartório
### This file and its contents are licensed under the terms of the MIT License


# @TODO: Add documentation
#' dropBurnin
#'
#' @export
dropBurnin <- function(x, ...) UseMethod("dropBurnin")


# @TODO: Add documentation
#' dropBurnin.MCMCresult2
#'
#' @export
dropBurnin.MCMCresult2 <- function(x, drop = 0.5) {
	if (!is.numeric(drop) || drop < 0)
		stop("Expected drop argument to be within [0,1] or [1, +Inf]")

	if (drop < 1) {
		n <- length(x$theta_1)
		cutoff <- floor(n * drop) + 1
		x$theta_1 <- x$theta_1[cutoff:n]
		x$theta_2 <- x$theta_2[cutoff:n]
	} else {
		n <- length(x$theta_1)
		x$theta_1 <- x$theta_1[drop:n]
		x$theta_2 <- x$theta_2[drop:n]
	}

	x
}


# @TODO: Add documentation
#' dropBurnin.MCMCresultSLR
#'
#' @export
dropBurnin.MCMCresultSLR <- function(x, drop = 0.5)
{
	if (!is.numeric(drop) || drop < 0)
		stop("Expected drop argument to be within [0,1] or [1, +Inf]")

	if (drop < 1) {
		n <- length(x$beta_0)
		cutoff <- floor(n * drop) + 1
		x$beta_0   <- x$beta_0[cutoff:n]
		x$beta_1   <- x$beta_1[cutoff:n]
		x$mu_x     <- x$mu_x[cutoff:n]
		x$sigma2_x <- x$sigma2_x[cutoff:n]
		x$sigma2_r <- x$sigma2_r[cutoff:n]
		x$sigma2_m <- x$sigma2_m[cutoff:n]
	} else {
		n <- length(x$beta_1)
		x$beta_0   <- x$beta_0[drop:n]
		x$beta_1   <- x$beta_1[drop:n]
		x$mu_x     <- x$mu_x[drop:n]
		x$sigma2_x <- x$sigma2_x[drop:n]
		x$sigma2_r <- x$sigma2_r[drop:n]
		x$sigma2_m <- x$sigma2_m[drop:n]
	}

	x
}
