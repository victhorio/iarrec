### lgauss.R
### Copyright (c) 2017 Victhor S. Sartório
### This file and its contents are licensed under the terms of the MIT License
###
### This file implements R wrappers for the C functions relating with the Gaussian
### distribution.


#' Log-Density of the Gaussian Distribution
#'
#' Vectorized log-density function for the Gaussian distribution
#'
#' @param x A numeric vector of the points at which the function is to be evaluated
#' @param mu The mean of the distribution
#' @param sigma2 The variance of the distribution
#' @return The log-density of a Gaussian variable at the points x
#'
#' @examples
#' lgauss(seq(0.1, 100, 0.2), 3, 5)
#'
#' @useDynLib identifiability
#' @export
lgauss <- function(x, mu, sigma2)
	.Call('lgauss', x, mu, sigma2)
