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
#' @return A list containing:
#' \describe{
#' \item{lambda}{The updated lambda matrix}
#' \item{phi}{The updated phi vector}
#' }
#' @note For Internal Use
m_step = function(dxx, e_obj, n) {
  lambda = update_lambda(e_obj$xez, e_obj$ezz)
  phi = update_phi(dxx, n, e_obj$xez, lambda)
  return(ls = list(lambda = lambda, phi = phi))
}

#' EM Algorithm for Factor Analsysis
#'
#' @description
#' Performs Factor Analysis given a multivariate sample. For very large data sets
#' (e.g. both sample size and dimension > 2000), the function can achieve
#' substantial speed-ups by model initialization with the "
#' "augmented implicitly restarted Lanczos bidiagonalization algorithm" (IRLBA).
#' User can turn this on by setting "breath_of_lightning = T".
#'
#' @param x Sample matrix (dimension N by P)
#' @param n_factors An integer: the number of factors
#' @param n_iter Number of EM iterations
#' @param tol Tolerance for convergence
#' @param ini_method Initialization by PCA via SVD ("pca") or by random factor loading ("random")
#' @param breath_of_lightning Whether to use irlba for fast approximate svd initialization
#' @param verbose Whether to display EM updates
#' @return A list containing:
#' \describe{
#'   \item{loadings}{The factor loadings estimate.}
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
#' fit = factor_analyzer(x, 2)
#' plot(fit$crit, type = "l")
#' fit$loadings
#' varimax(fit$loadings)
#' promax(fit$loadings)
#' @export
#' @import irlba
factor_analyzer = function(x, n_factors, n_iter = 200, tol = 1e-6,
                           ini_method = c("pca", "random"), breath_of_lightning = F, verbose = F) {

  ini_method = match.arg(ini_method)
  if (ini_method == "pca") {
    if (breath_of_lightning) svd_fit = irlba::irlba(x, nu = 0, nv = n_factors)
    else svd_fit = svd(x, nu = 0, nv = n_factors)
    lambda = as.matrix(svd_fit$v[,1:n_factors])
    lambda = t(t(lambda)*(svd_fit$d[1:n_factors]/sqrt(nrow(x))))
  }
  else lambda = matrix(rnorm(ncol(x) * n_factors, sd = 0.1), ncol(x), n_factors)
  dxx = colSums(x^2)
  phi = pmax(dxx/nrow(x) - rowSums(lambda^2), 1e-2)
  crit = rep(NA, n_iter)

  for (i in 1:n_iter) {
    e_obj = e_step(x, lambda, phi)
    crit[i] = loglik(x, dxx, lambda, phi)
    if (i > 1) if ((crit[i]-crit[i-1]) < tol) break
    m_obj = m_step(dxx, e_obj, nrow(x))
    lambda = m_obj$lambda
    phi = m_obj$phi
    if (verbose) message("iter = ", i, ", crit = ", crit[i])
  }

  # if (is.null(rotation)) loadings = lambda
  # else loadings = rotation(lambda, ...)

  loadings = lambda

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
  -nrow(x)*ncol(x)*log(2*pi)/2 - nrow(x)/2 * (log_det + (trace_1 - trace_2)/nrow(x))
}

#' Akaike Information Criterion (AIC) for Factor Analysis
#'
#' @param fit A fitted object from factor_analyzer()
#' @return The AIC
#' @note The AIC formula used is from Akaike (1987) Pyschometrika, 52(3):317-322.
#' It is not exactly the same as the typical -2(loglik) + 2(num_parameter).
#' @examples
#' set.seed(8)
#' z = matrix(rnorm(2000), 1000, 2)
#' lambda_orig = matrix(rnorm(20, sd = 1), nrow = 10, ncol = 2)
#' x = z %*% t(lambda_orig) + matrix(rnorm(10000, sd = sqrt(0.1)), 1000, 10)
#' aic_vec = rep(0, 5)
#' for (i in seq_along(aic_vec)) aic_vec[i] = aic_fa(factor_analyzer(x, i, rotation = NULL))
#' which.min(aic_vec)
#' @export
aic_fa = function(fit) {
  p = nrow(fit$loadings)
  k = ncol(fit$loadings)
  -2*fit$crit[length(fit$crit)] + 2*(p*(k + 1)) - k*(k-1)
}


#' Likeilood Ratio Test for Goodness of Fit
#'
#' @param cov_x The sample covariance matrix
#' @param n The sample size
#' @param fit A fitted object from factor_analyzer()
#' @return A list contain the chi-square statistics (chi_sq), the degree-of-freedom (df), and the p-value.
#' @examples
#' set.seed(8)
#' z = matrix(rnorm(3000), 1000, 3)
#' lambda_orig = matrix(rnorm(30, sd = 1), nrow = 10, ncol = 3)
#' x = z %*% t(lambda_orig) + matrix(rnorm(10000, sd = sqrt(0.1)), 1000, 10)
#' s = cov(x)/nrow(x)*(nrow(x) - 1)
#' p_vec = rep(0, 5)
#' for (i in seq_along(p_vec)) p_vec[i] = lrt_fa(s, nrow(x), factor_analyzer(x, i, rotation = NULL))$p_val
#' p_vec
#' @export
lrt_fa = function(cov_x, n, fit) {
  p = ncol(cov_x)
  lambda = fit$loadings
  phi = fit$phi
  m = ncol(lambda)
  phi_lambda = lambda / phi
  log_det_fa = sum(log(phi)) + log(det(diag(1, ncol(lambda)) + crossprod(lambda, phi_lambda)))
  log_det_s = log(det(cov_x))
  # chi_sq = (n - 1 - (2*p + 5)/6 - 2*m/3) * (log_det_fa - log_det_s)
  inverse_fa = fast_inverse(lambda, phi)
  chi_sq = n * (log_det_fa - log_det_s + sum(cov_x*inverse_fa) - p)
  df = ((p-m)^2 - p - m)/2
  return(ls = list(chi_sq = chi_sq, df = df, p_val = pchisq(chi_sq, df = df, lower.tail = F)))
}

#' Factor Scores
#'
#' @param fit A fitted object from factor_analyzer()
#' @param newdata A new data set for scores computation
#' @return The projected scores of newdata
#' @examples
#' set.seed(8)
#' z = matrix(rnorm(3000), 1000, 3)
#' z_2 = matrix(rnorm(3000), 1000, 3)
#' lambda_orig = matrix(rnorm(30, sd = 1), nrow = 10, ncol = 3)
#' x = z %*% t(lambda_orig) + matrix(rnorm(10000, sd = sqrt(0.1)), 1000, 10)
#' newdata = z_2 %*% t(lambda_orig) + matrix(rnorm(10000, sd = sqrt(0.1)), 1000, 10)
#' fit = factor_analyzer(x, 2, rotation = NULL)
#' scores = get_scores(fit, newdata)
#' @export
get_scores = function(fit, newdata) {
  beta = fast_get_beta(fit$loadings, fit$phi)
  expected_scores(newdata, beta)
}



#' EM Algorithm for Factor Analsysis with Covariance Input
#'
#' @description
#' Performs Factor Analysis given a covariance/correlation matrix.
#'
#' @param cov_x Sample covariance/correlation matrix
#' @param n_factors An integer: the number of factors
#' @param covar Whether input cov_x is a covariance (T) or correlation (F) matrix
#' @param n_sample Number of samples to draw from Multivariate Gaussian
#' @param ini_method The EM initialization method ("pca" recommended)
#' @param breath_of_lighting Whether to use irlba for initialization (recommended)
#' @param ... Other parameters passed to the "factor_analyzer" function
#' @return A fitted factor analyzer object
#'
#' @note The function runs EM on a simulated multivariate Gaussian sample
#' to obtain approximate maximum likelihood estimates. n_sample, the number of
#' samples to simulate, should be large (200000 based on our benchmark experiments).
#' As EM will be performed on this large data set, initializing EM by
#' pca + breath_of_lightning is recommended.
#'
#' @examples
#' \dontrun{
#' set.seed(8)
#' z = matrix(rnorm(2000), 1000, 2)
#' lambda_orig = matrix(rnorm(20, sd = 1), nrow = 10, ncol = 2)
#' x = z %*% t(lambda_orig) + matrix(rnorm(10000, sd = sqrt(0.1)), 1000, 10)
#'
#' fit = factor_analyzer(x, 2, rotation = NULL)
#' varimax(fit$loadings)
#' cov_x = cov(x)*(nrow(x)-1)/nrow(x)
#' fit_cov = factor_analyzer_cov(cov_x, 2, covar = T, n_sample = 200000)
#' varimax(fit_cov$loadings)
#'
#' x = scale(x)
#' fit = factor_analyzer(x, 2, rotation = NULL)
#' varimax(fit$loadings)
#' cor_x = cor(x)
#' fit_cor = factor_analyzer_cov(cor_x, 2, covar = F, n_sample = 200000)
#' varimax(fit_cor$loadings)
#' }
#' @export
#' @import mvtnorm
factor_analyzer_cov = function(cov_x, n_factors, covar = NULL, n_sample = 10000,
                           ini_method = "pca", breath_of_lightning = T, ...) {

  if (is.null(covar)) stop("covar must be specified (Is cov_x a covariance (T) or correlation matrix (F)?)")
  x_samp = mvtnorm::rmvnorm(n_sample, sigma = cov_x)
  if (!covar) x_samp = scale(x_samp)
  factor_analyzer(x_samp, n_factors = n_factors,
                  ini_method = ini_method, breath_of_lightning = breath_of_lightning, ...)

}
