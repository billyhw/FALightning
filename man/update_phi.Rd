% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FALightning.R
\name{update_phi}
\alias{update_phi}
\title{Phi Update}
\usage{
update_phi(dxx, n, xez, lambda)
}
\arguments{
\item{dxx}{Column sums-of-squares of the sample matrix}

\item{n}{Sample size}

\item{xez}{Covariance between observed and latent factors}

\item{lambda}{Factor loading matrix (dimension P by P')}
}
\value{
The updated phi vector
}
\description{
Phi Update
}
\note{
For Internal Use
}
