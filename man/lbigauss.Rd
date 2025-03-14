% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gauss.R
\name{lbigauss}
\alias{lbigauss}
\title{Log-Density of the Bivariate Gaussian Distribution}
\usage{
lbigauss(x, y, mx, my, sxx, syy, sxy)
}
\arguments{
\item{x}{The point x of the function}

\item{y}{The point y of the function}

\item{mx}{The mean of the X distribution}

\item{my}{The mean of the Y distribution}

\item{sxx}{The variance of the X distribution}

\item{syy}{The variance of the Y distribution}

\item{sxy}{The covariance between X and Y}
}
\value{
The log-density of a bivariate Gaussian variable (X,Y) at the point (x,y)
}
\description{
Atomic log-density function for the Gaussian distribution
}
\examples{
lbigauss(0.3, 0.8, 6, 3, 8, 12, 4)

}
