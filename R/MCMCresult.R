# @TODO: Add documentation
#' print.MCMCresult
#'
#' @export
print.MCMCresult <- function(x)
	cat(sprintf('MCMC Result containing %d data points and acceptance ratio %f\n',
				length(x[[1]]), x[[3]]))


#' Plots a MCMC Result from the IARREC package
#'
#' @param x An MCMCresult S3 object
#' @param drop Indicates how much of the samples to ignore.
#'     If it's a value between [0,1) then it ignores the first \eqn{drop * n} samples from
#'     the result. If it's a value within [1, +Inf) then the first \eqn{drop} values are
#'     ignored. Ignores nothing by default.
#' @param ... Optional base R graphical parameters
#'
#' @export
plot.MCMCresult <- function(x, drop = 0, ...)
{
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
