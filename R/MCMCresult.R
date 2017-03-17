### MCMCresult.R
### Copyright (c) 2017 Victhor S. Sartório
### This file and its contents are licensed under the terms of the MIT License


# @TODO: Add documentation
#' print.MCMCresultSLR
#'
#' @export
print.MCMCresultSLR <- function(x)
	cat(sprintf('MCMC Result containing %d data points and acceptance ratio %f\n',
				length(x[[1]]), x[[7]]))


# @TODO: Add documentation
#' mean.MCMCresultSLR
#'
#' @export
mean.MCMCresultSLR <- function(x)
{
	nms <- names(x)[1:6]
	r <- rep(0, 6)
	for (i in 1:6) {
		r[i] <- mean(x[[i]])
	}
	names(r) <- nms
	r
}


# @TODO: Add documentation
#' print.MCMCresult2
#'
#' @export
print.MCMCresult2 <- function(x)
	cat(sprintf('MCMC Result containing %d data points and acceptance ratio %f\n',
				length(x[[1]]), x[[3]]))


# @TODO: Add documentation
#' mean.MCMCresultSLR
#'
#' @export
mean.MCMCresult2 <- function(x)
{
	nms <- names(x)[1:2]
	r <- rep(0, 2)
	for (i in 1:2) {
		r[i] <- mean(x[[i]])
	}
	names(r) <- nms
	r
}
