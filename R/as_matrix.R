### as_matrix.R
### Copyright (c) 2017 Victhor S. Sartório
### This file and its contents are licensed under the terms of the MIT License


# @TODO: Add documentation
#' as.matrix.MCMCresult2
#'
#' @export
as.matrix.MCMCresult2 <- function(x)
{
	r <- cbind(x$theta_1, x$theta_2)
	colnames(r) <- names(x)[1:2]
	r
}


# @TODO: Add documentation
#' as.matrix.MCMCresultSLR
#'
#' @export
as.matrix.MCMCresultSLR <- function(x)
{
	r <- cbind(x$beta_0, x$beta_1, x$mu_x, x$sigma2_x, x$sigma2_r, x$sigma2_m)
	colnames(r) <- names(x)[1:6]
	r
}
