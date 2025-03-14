% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gauss.R
\name{lgauss}
\alias{lgauss}
\title{Log-Density of the Gaussian Distribution}
\usage{
lgauss(x, mu, sigma2)
}
\arguments{
\item{x}{A numeric vector of the points at which the function is to be evaluated}

\item{mu}{The mean of the distribution}

\item{sigma2}{The variance of the distribution}
}
\value{
The log-density of a Gaussian variable at the points x
}
\description{
Vectorized log-density function for the Gaussian distribution
}
\examples{
lgauss(seq(0.1, 100, 0.2), 3, 5)

}
