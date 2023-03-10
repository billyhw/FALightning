% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FALightning.R
\name{m_step}
\alias{m_step}
\title{EM Algorithm for Factor Analsysis: M-step}
\usage{
m_step(dxx, e_obj, n)
}
\arguments{
\item{dxx}{Column sums-of-squares of the sample matrix}

\item{e_obj}{A list of expected values from the E-Step}

\item{n}{Sample size}
}
\value{
A list containing:
\describe{
\item{lambda}{The updated lambda matrix}
\item{phi}{The updated phi vector}
}
}
\description{
EM Algorithm for Factor Analsysis: M-step
}
\note{
For Internal Use
}
