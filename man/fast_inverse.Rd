% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FALightning.R
\name{fast_inverse}
\alias{fast_inverse}
\title{Fast Computation of Inverse Covariance Matrix}
\usage{
fast_inverse(lambda, phi)
}
\arguments{
\item{lambda}{Factor loading matrix (dimension P by P')}

\item{phi}{Vector of noise variance}
}
\value{
The inverse of diag(phi) + lambda*lambda^T
}
\description{
Fast Computation of Inverse Covariance Matrix
}
\note{
For Internal Use. Deprecated as of Jan 5, 2023
}
