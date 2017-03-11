### invgamma.R
### Copyright (c) 2017 Victhor S. Sartório
### This file and its contents are licensed under the terms of the MIT License
###
### This file implements R wrappers for the C functions relating with the inverse gamma
### distribution.


#' Log-Density of the Inverse Gamma Distribution
#'
#' Vectorized log-density function for the inverse gamma distribution
#'
#' @param x A numeric vector of the points at which the function is to be evaluated
#' @param alpha The alpha parameter of the inverse gamma distribution
#' @param beta The beta parameter of the inverse gamma distribution
#' @return The log-density of an inverse-gamma variable at the points x
#'
#' @examples
#' linvgamma(seq(0.1, 100, 0.2), 3, 5)
#'
#' @useDynLib iarrec
#' @export
linvgamma <- function(x, alpha, beta)
	.Call('linvgamma', x, alpha, beta)


#' Density of the Inverse Gamma Distribution
#'
#' Vectorized density function for the inverse gamma distribution
#'
#' @param x A numeric vector of the points at which the function is to be evaluated
#' @param alpha The alpha parameter of the inverse gamma distribution
#' @param beta The beta parameter of the inverse gamma distribution
#' @return The density of an inverse-gamma variable at the points x
#'
#' @examples
#' dinvgamma(seq(0.1, 100, 0.2), 3, 5)
#'
#' @useDynLib iarrec
#' @export
dinvgamma <- function(x, alpha, beta)
	.Call('dinvgamma', x, alpha, beta)
