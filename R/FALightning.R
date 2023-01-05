#' Fast Computation of Inverse Covariance Matrix
#'
#' @param lambda Factor loading matrix (dimension P by P')
#' @param phi Vector of noise variance
#' @return The inverse of diag(phi) + lambda*lambda^T
#' @note For Internal Use. Deprecated as of Jan 5, 2023
fast_inverse = function(lambda, phi) {
  phi_inv = 1/phi
  phi_lambda = lambda * phi_inv
  inverse = solve(diag(ncol(lambda)) + crossprod(lambda, phi_lambda))
  inverse = -tcrossprod(phi_lambda %*% inverse, phi_lambda)
  diag(inverse) = diag(inverse) + phi_inv
  inverse
}

#' Coefficients for Evaluating the Expected Latent Factors
#'
#' @param lambda Factor loading matrix (dimension P by P')
#' @param inverse Inverse covariance matrix
#' @return A set of coefficients for evaluating the expected latent factors
#' @note For Internal Use. Deprecated as of Jan 5, 2023
get_beta = function(lambda, inverse) crossprod(lambda, inverse)

#' Coefficients for Evaluating the Expected Latent Factors: Fast Calculation
#'
#' @description
#' This function merged two older functions get_beta() and fast_inverse().
#' This allows speed ups as the entire inverse covariance is never necessary
#' in the EM fitting procedure. Speed up is thus achieved by only
#' performing the calculations necessary for the coefficients.
#'
#' @param lambda Factor loading matrix (dimension P by P')
#' @param phi Vector of noise variance
#' @return A set of coefficients for evaluating the expected latent factors
#' @note For Internal Use
fast_get_beta = function(lambda, phi) {
  phi_lambda = lambda / phi
  inverse = solve(diag(ncol(lambda)) + crossprod(lambda, phi_lambda))
  beta = crossprod(lambda, phi_lambda)
  beta = tcrossprod(beta %*% inverse, phi_lambda)
  t(phi_lambda) - beta
}

#' Expected Latent Scores
#'
#' @param x Sample matrix (dimension N by P)
#' @param beta Coefficients from get_beta()
#' @return Expected latent scores conditioned on observed data
#' @note For Internal Use
expected_scores = function(x, beta) tcrossprod(x, beta)

#' Expected Covariance of the Latent Factors
#'
#' @param ez Expected latent scores from expected_scores()
#' @param beta Coefficients from get_beta()
#' @param lambda Factor loading matrix (dimension P by P')
#' @return The covariance of the expected latent scores
#' @note For Internal Use
expected_cov = function(ez, beta, lambda) {
  beta_lambda = beta %*% lambda
  zz = crossprod(ez, ez)
  diag(nrow(ez), nrow(zz)) - nrow(ez) * beta_lambda + zz
}

#' Covariance Between Observed Samples and Latent Factors
#'
#' @param x Sample matrix (dimension N by P)
#' @param ez Expected latent scores from expected_scores()
#' @return Covariance between observed and latent scores
#' @note For Internal Use
get_xEz = function(x, ez) crossprod(x, ez)

#' Lambda Update
#'
#' @param xez Covariance between observed and latent factors
#' @param ezz Expected covariance of the latent factors from expected_cov()
#' @return The updated lambda matrix
#' @note For Internal Use
update_lambda = function(xez, ezz) xez %*% solve(ezz)

#' Phi Update
#'
#' @param dxx Column sums-of-squares of the sample matrix
#' @param n Sample size
#' @param xez Covariance between observed and latent factors
#' @param lambda Factor loading matrix (dimension P by P')
#' @return The updated phi vector
#' @note For Internal Use
update_phi= function(dxx, n, xez, lambda) (dxx - rowSums(lambda*xez))/n

#' EM Algorithm for Factor Analsysis: E-step
#'
#' @param x Sample matrix (dimension N by P')
#' @param lambda Factor loading matrix (dimension P by P')
#' @param phi Vector of noise variance
#' @return A list containing:
#' \describe{
#' \item{ez}{Expected latent scores from expected_scores()}
#' \item{ezz}{Expected covariance of the latent factors from expected_cov()}
#' \item{xez}{Covariance between observed and latent factors}
#' }
#' @note For Internal Use
e_step = function(x, lambda, phi) {
  beta = fast_get_beta(lambda, phi)
  ez = expected_scores(x, beta)
  ezz = expected_cov(ez, beta, lambda)
  xez = get_xEz(x, ez)
  return(ls = list(ez = ez, ezz = ezz, xez = xez))
}

#' EM Algorithm for Factor Analsysis: M-step
#'
#' @param dxx Column sums-of-squares of the sample matrix
#' @param e_obj A list of expected values from the E-Step
#' @param n Sample size
#' @param lambda The factor loadings
#' @return A list containing:
#' \describe{
#' \item{lambda}{The updated lambda matrix}
#' \item{phi}{The updated phi vector}
#' }
#' @note For Internal Use
m_step = function(dxx, e_obj, n, lambda) {
  lambda = update_lambda(e_obj$xez, e_obj$ezz)
  phi = update_phi(dxx, n, e_obj$xez, lambda)
  return(ls = list(lambda = lambda, phi = phi))
}

#' EM Algorithm for Factor Analsysis
#'
#' @param x Sample matrix (dimension N by P)
#' @param n_factors An integer: the number of factors
#' @param n_iter Number of EM iterations
#' @param tol Tolerance for convergence
#' @param rotation A function for factor rotation, e.g. varimax or oblimin from the GPArotation package
#' @param breath_of_lightning Whether to use irlba for fast approximate svd initialization
#' @param verbose Whether to display EM updates
#' @param ... Other parameters passed to the "rotation" function
#' @return A list containing:
#' \describe{
#'   \item{loadings}{If rotation = NULL, the factor loadings estimate. Otherwise, the object returned by the "rotation" function.}
#'   \item{phi}{The phi vector estimate}
#'   \item{scores}{The expected factor scores}
#'   \item{crit}{The log-likelihood criterion per iteration}
#'   \item{converge_status}{Convergence status: 0 = converged, 1 = tolerance not reached, 2 = likelihood has decreased}
#' }
#'
#' @examples
#' set.seed(8)
#' z = matrix(rnorm(2000), 1000, 2)
#' lambda_orig = matrix(rnorm(20, sd = 1), nrow = 10, ncol = 2)
#' x = z %*% t(lambda_orig) + matrix(rnorm(10000, sd = sqrt(0.1)), 1000, 10)
#' fit = factor_analyzer(x, 2, rotation = NULL)
#' plot(fit$crit, type = "l")
#' fit$loadings
#' fit_varimax = factor_analyzer(x, 2)
#' fit_varimax$loadings
#' fit_promax = factor_analyzer(x, 2, rotation = promax)
#' fit_promax$loadings
#' @export
#' @import irlba
factor_analyzer = function(x, n_factors, n_iter = 200, tol = 1e-4, rotation = varimax,
                           breath_of_lightning = F, verbose = F, ...) {

  if (breath_of_lightning) svd_fit = irlba::irlba(x, nu = 0, nv = n_factors)
  else svd_fit = svd(x, nu = 0, nv = n_factors)
  lambda = as.matrix(svd_fit$v[,1:n_factors])
  lambda = t(t(lambda)*(svd_fit$d[1:n_factors]/sqrt(nrow(x))))
  dxx = colSums(x^2)
  phi = pmax(dxx/nrow(x) - rowSums(lambda^2), 1e-2)
  crit = rep(NA, n_iter)

  for (i in 1:n_iter) {
    e_obj = e_step(x, lambda, phi)
    crit[i] = loglik(x, dxx, lambda, phi)
    if (i > 1) if ((crit[i]-crit[i-1]) < tol) break
    m_obj = m_step(dxx, e_obj, nrow(x), lambda)
    lambda = m_obj$lambda
    phi = m_obj$phi
    if (verbose) message("iter = ", i, ", crit = ", crit[i])
  }

  if (is.null(rotation)) loadings = lambda
  else loadings = rotation(lambda, ...)

  crit = crit[!is.na(crit)]
  converge_status = 0
  if ((crit[i]-crit[i-1]) > tol) converge_status = 1
  if (any(diff(crit) < 0)) coverge_status = 2

  return(ls = list(loadings = loadings, phi = phi, scores = e_obj$ez,
                   crit = crit, converge_status = converge_status))

}

#' Log-Likelihood
#'
#' @param x Sample matrix (dimension N by P)
#' @param dxx Column sums-of-squares of the sample matrix
#' @param lambda Factor loading matrix (dimension P by P')
#' @param phi Vector of noise variance
#' @return The log-likelihood criterion
#' @note For Internal Use
loglik = function(x, dxx, lambda, phi) {
  phi_lambda = lambda / phi
  log_det = sum(log(phi)) + log(det(diag(1, ncol(lambda)) + crossprod(lambda, phi_lambda)))
  trace_1 = crossprod(dxx, 1/phi)
  inverse = solve(diag(ncol(lambda)) + crossprod(lambda, phi_lambda))
  xlp = x %*% phi_lambda
  trace_2 = xlp %*% inverse
  trace_2 = sum(trace_2 * xlp)
  -nrow(x)/2 * (log_det + (trace_1 - trace_2)/nrow(x))
}
