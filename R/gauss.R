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
#' @useDynLib iarrec
#' @export
lgauss <- function(x, mu, sigma2)
	.Call('lgauss', x, mu, sigma2)


#' Log-Density of the Bivariate Gaussian Distribution
#'
#' Atomic log-density function for the Gaussian distribution
#'
#' @param x The point x of the function
#' @param y The point y of the function
#' @param mx The mean of the X distribution
#' @param my The mean of the Y distribution
#' @param sxx The variance of the X distribution
#' @param syy The variance of the Y distribution
#' @param sxy The covariance between X and Y
#' @return The log-density of a bivariate Gaussian variable (X,Y) at the point (x,y)
#'
#' @examples
#' lbigauss(0.3, 0.8, 6, 3, 8, 12, 4)
#'
#' @useDynLib iarrec
#' @export
lbigauss <- function(x, y, mx, my, sxx, syy, sxy)
	.Call('lbigauss', x, y, mx, my, sxx, syy, sxy)
