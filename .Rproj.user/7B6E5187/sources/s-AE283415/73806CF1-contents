#' Fast Computation of Inverse Covariance Matrix
#'
#' @param lambda Factor loading matrix (dimension P by P')
#' @param phi Vector of noise variance
#' @return The inverse of diag(phi) + lambda*lambda^T
#' @note For Internal Use
fast_inverse = function(lambda, phi) {
  phi_inv = 1/phi
  phi_lambda = lambda * phi_inv
  inverse = solve(diag(ncol(lambda)) + crossprod(lambda, phi_lambda))
  diag(phi_inv) - tcrossprod(phi_lambda %*% inverse, phi_lambda)
}

#' Coefficients for Evaluating the Expected Latent Factors
#'
#' @param lambda Factor loading matrix (dimension P by P')
#' @param inverse Inverse covariance matrix
#' @return A set of coefficients for evaluating the expected latent factors
#' @note For Internal Use
get_beta = function(lambda, inverse) crossprod(lambda, inverse)

#' Expected Latent Scores
#'
#' @param x Sample matrix (dimension N by P')
#' @param beta Coefficients from get_beta()
#' @return Expected latent scores conditioned on observed data
#' @note For Internal Use
expected_scores = function(x, beta) x %*% t(beta)

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
#' @param x Sample matrix (dimension N by P')
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
#' @param x Sample matrix (dimension N by P')
#' @param xez Covariance between observed and latent factors
#' @param lambda Factor loading matrix (dimension P by P')
#' @return The updated phi vector
#' @note For Internal Use
update_phi= function(x, xez, lambda) diag(crossprod(x) - tcrossprod(lambda, xez))/nrow(x)

#' EM Algorithm for Factor Analsysis: E-step
#'
#' @param x Sample matrix (dimension N by P')
#' @param beta Coefficients from get_beta()
#' @return A list containing:
#'
#' ez: Expected latent scores from expected_scores()
#'
#' ezz: Expected covariance of the latent factors from expected_cov()
#'
#' xez: Covariance between observed and latent factors
#'
#' @note For Internal Use
e_step = function(x, lambda, phi) {
  inverse = fast_inverse(lambda, phi)
  beta = get_beta(lambda, inverse)
  ez = expected_scores(x, beta)
  ezz = expected_cov(ez, beta, lambda)
  xez = get_xEz(x, ez)
  return(ls = list(ez = ez, ezz = ezz, xez = xez))
}

#' EM Algorithm for Factor Analsysis: M-step
#'
#' @param x Sample matrix (dimension N by P')
#' @param e_obj A list of expected values from the E-Step
#' @return A list containing:
#'
#' lambda: The updated lambda matrix
#'
#' phi: The updated phi vector
#'
#' @note For Internal Use
m_step = function(x, e_obj) {
  lambda = update_lambda(e_obj$xez, e_obj$ezz)
  phi = update_phi(x, e_obj$xez, lambda)
  return(ls = list(lambda = lambda, phi = phi))
}

#' EM Algorithm for Factor Analsysis
#'
#' @param x Sample matrix (dimension N by P')
#' @param n_factors An integer: the number of factors
#' @param n_iter Number of EM iterations
#' @param rotation A function for factor rotation, e.g. varimax or oblimin from the GPArotation package
#' @param ... Other parameters passed to the "rotation" function
#' @return A list containing:
#' \describe{
#'   \item{loadings}{If rotation = NULL, the factor loadings estimate. Otherwise, the object returned by the "rotation" function.}
#'   \item{phi}{The phi vector estimate}
#'   \item{ez}{The expected factor scores}
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
factor_analyzer = function(x, n_factors, n_iter = 200, rotation = varimax, ...) {

  lambda = svd(x)$v[,1:n_factors]
  phi = diag(cov(x))
  crit = rep(NA, n_iter)

  for (i in 1:n_iter) {
    e_obj = e_step(x, lambda, phi)
    m_obj = m_step(x, e_obj)
    lambda = m_obj$lambda
    phi = m_obj$phi
    crit[i] = mean((tcrossprod(lambda) + diag(phi) - cov(x))^2)
  }

  e_obj = e_step(x, lambda, phi)

  if (is.null(rotation)) loadings = lambda
  else loadings = rotation(lambda, ...)

  return(ls = list(loadings = loadings, phi = phi, scores = e_obj$ez, crit = crit))

}
