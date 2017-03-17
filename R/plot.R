### plot.R
### Copyright (c) 2017 Victhor S. Sartório
### This file and its contents are licensed under the terms of the MIT License


#' Plots a MCMC Result from the IARREC package
#'
#' @param x An MCMCresult S3 object
#' @param ... Optional base R graphical parameters
#'
#' @export
plot.MCMCresult2 <- function(x, ...)
{
	qt1 <- quantile(x$theta_1, c(0.25, 0.5, 0.75))
	qt2 <- quantile(x$theta_2, c(0.25, 0.5, 0.75))

	old.pars <- par(no.readonly = TRUE)
	par(...)

	layout(rbind(c(1,2), c(3, 4)))
	plot(x$theta_1, type = 'l', main = 'Chain for Theta 1', ylab = '')
	abline(h = qt1, col = 'orange', lty = 2, lwd = 1.5)
	plot(x$theta_2, type = 'l', main = 'Chain for Theta 2', ylab = '')
	abline(h = qt2, col = 'orange', lty = 2, lwd = 1.5)
	plot(density(x$theta_1), main = 'Density Estimate of Theta 1')
	abline(v = qt1, col = 'orange', lty = 2, lwd = 1.5)
	plot(density(x$theta_2), main = 'Density Estimate of Theta 2')
	abline(v = qt2, col = 'orange', lty = 2, lwd = 1.5)

	par(old.pars)
}


#' Plots a MCMC Result from the IARREC package
#'
#' @param x An MCMCresult S3 object
#' @param pause Indicates if should pause in between plots. Defaults to TRUE.
#' @param ... Optional base R graphical parameters
#'
#' @export
plot.MCMCresultSLR <- function(x, pause = TRUE, ...)
{
	q_beta0 <- quantile(x$beta_0, c(0.25, 0.5, 0.75))
	q_beta1 <- quantile(x$beta_1, c(0.25, 0.5, 0.75))
	q_mux   <- quantile(x$mu_x, c(0.25, 0.5, 0.75))
	q_s2x   <- quantile(x$sigma2_x, c(0.25, 0.5, 0.75))
	q_s2r   <- quantile(x$sigma2_r, c(0.25, 0.5, 0.75))
	q_s2m   <- quantile(x$sigma2_m, c(0.25, 0.5, 0.75))

	old.pars <- par(no.readonly = TRUE)
	par(...)

	layout(rbind(c(1,2), c(3, 4)))
	plot(x$beta_0, type = 'l', main = 'Chain for Beta 0', ylab = '')
	abline(h = q_beta0, col = 'orange', lty = 2, lwd = 1.5)
	plot(x$beta_1, type = 'l', main = 'Chain for Beta 1', ylab = '')
	abline(h = q_beta1, col = 'orange', lty = 2, lwd = 1.5)
	plot(density(x$beta_0), main = 'Density Estimate of Beta 0')
	abline(v = q_beta0, col = 'orange', lty = 2, lwd = 1.5)
	plot(density(x$beta_1), main = 'Density Estimate of Beta 1')
	abline(v = q_beta1, col = 'orange', lty = 2, lwd = 1.5)

	if (pause)
		readline("Press <enter> to continue...")

	layout(rbind(c(1,2), c(3, 4)))
	plot(x$mu_x, type = 'l', main = 'Chain for Mu x', ylab = '')
	abline(h = q_mux, col = 'orange', lty = 2, lwd = 1.5)
	plot(x$sigma2_x, type = 'l', main = 'Chain for Sigma2 X', ylab = '')
	abline(h = q_s2x, col = 'orange', lty = 2, lwd = 1.5)
	plot(density(x$mu_x), main = 'Density Estimate of Mu x')
	abline(v = q_mux, col = 'orange', lty = 2, lwd = 1.5)
	plot(density(x$sigma2_x), main = 'Density Estimate of Sigma2 x')
	abline(v = q_s2x, col = 'orange', lty = 2, lwd = 1.5)

	if (pause)
		readline("Press <enter> to continue...")

	layout(rbind(c(1,2), c(3, 4)))
	plot(x$sigma2_r, type = 'l', main = 'Chain for Sigma2 r', ylab = '')
	abline(h = q_s2r, col = 'orange', lty = 2, lwd = 1.5)
	plot(x$sigma2_m, type = 'l', main = 'Chain for Sigma2 m', ylab = '')
	abline(h = q_s2m, col = 'orange', lty = 2, lwd = 1.5)
	plot(density(x$sigma2_r), main = 'Density Estimate of Sigma2 r')
	abline(v = q_s2r, col = 'orange', lty = 2, lwd = 1.5)
	plot(density(x$sigma2_m), main = 'Density Estimate of Sigma2 m')
	abline(v = q_s2m, col = 'orange', lty = 2, lwd = 1.5)

	par(old.pars)
}
