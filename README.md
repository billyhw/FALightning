# FALightning

<img src="/images/zenitsu.jpg" width="400" height="300">

# Overview

This R package contains a fast implementation of factor analysis based on the EM-algorithm. The implementation is based on (Ghahramani & Hinton 1997. EM Algorithm for Mixtures of Factor Analyzers. University of Toronto Technical Report).

The speed of this implemenation relies on (1) R-optimized implementation procedure and (2) the use of the matrix inversion formula (a.k.a. the Woodbury Matrix Identity):

$(\Lambda \Lambda^T + \Phi)^{-1} = \Phi^{-1} - \Phi^{-1}\Lambda(I + \Lambda^T \Phi^{-1} \Lambda)^{-1}\Lambda^T \Phi^{-1}

This formula allows the inverse of a large P-by-P covariance matrix, with P > 1000 often, to be evaluated by inverting instead the much smaller P'-by-P' matrix in the right hand side of the formula. Here P' is the number of factors, usually much less than P and rarely over 100 in practice.

# Installation

```
#install.packages("devtools")
devtools::install_github("billyhw/FALightning")
```
