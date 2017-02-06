#'
#' Generate synthetic data from a dedicated factor model
#'
#' This function simulates data from a dedicated factor model. The parameters of
#' the model are either passed by the user or simulated by the function.
#'
#' @param N
#'        Number of observations in data set.
#' @param dedic
#'        Vector of indicators. The number of manifest variables is equal to the
#'        length of this vector, and the number of factors is equal to the
#'        number of unique nonzero elements. Each integer element indicates on
#'        which latent factor the corresponding variable loads uniquely.
#' @param alpha
#'        Vector of factor loadings, should be of same length as \code{dedic}.
#'        If missing, values are simulated (see details below).
#' @param sigma
#'        Idiosyncratic variances, should be of same length as \code{dedic}.
#'        If missing, values are simulated (see details below).
#' @param R
#'        Covariance matrix of the latent factors.
#'        If missing, values are simulated (see details below).
#' @param R.corr
#'        If TRUE, covariance matrix \code{R} is rescaled to be a correlation
#'        matrix.
#' @param max.corr
#'        Maximum correlation allowed between the latent factors.
#' @param R.max.trial
#'        Maximum number of trials allowed to sample from the truncated
#'        distribution of the covariance matrix of the latent factors
#'        (accept/reject sampling scheme, to make sure \code{max.corr} is not
#'        exceeded).
#'
#' @details The function simulates data from the following dedicated factor
#' model, for \eqn{i = 1, ..., N}:
#'   \deqn{Y_i = \alpha \theta_i + \epsilon_i}
#'   \deqn{\theta_i \sim \mathcal{N}(0, R)}{\theta_i ~ N(0, R)}
#'   \deqn{\epsilon_i \sim \mathcal{N}(0, \Sigma)}{\epsilon_i ~ N(0, \Sigma)}
#' where the \eqn{K}-vector \eqn{\theta_i} contains the latent factors, and
#' \eqn{\alpha} is the \eqn{(M \times K)}{(M*K)}-matrix of factor loadings. Each
#' row \eqn{m} of \eqn{\alpha} contains only zeros, besides its element
#' indicated by the \eqn{m}th element of \code{dedic} that is equal to the
#' \eqn{m}th element of \code{alpha} (denoted \eqn{\alpha_m^\Delta} below).
#' The \eqn{M}-vector \eqn{\epsilon_i} is the vector of error terms, with
#' \eqn{\Sigma = diag(}\code{sigma}\eqn{)}. \eqn{M} is equal to the length of
#' the vector \code{dedic}, and \eqn{K} is equal to the maximum value of this
#' vector.
#'
#' Only \code{N} and \code{dedic} are required, all the other parameters can be
#' missing, completely or partially. Missing values (\code{NA}) are
#' independently sampled from the following distributions, for each manifest
#' variable \eqn{m = 1, ..., M}:
#'
#' Factor loadings:
#' \deqn{\alpha_m^\Delta = (-1)^{\phi_m}\sqrt{a_m}}{
#'       \alpha_m^\Delta = (-1)^\phi_m\sqrt(a_m)}
#' \deqn{\phi_m \sim \mathcal{B}er(0.5)}{\phi_m ~ Ber(0.5)}
#' \deqn{a_m \sim \mathcal{U}nif (0.04, 0.64)}{a_m ~ Unif (0.04, 0.64)}
#'
#' Idiosyncratic variances:
#' \deqn{\sigma^2_m \sim \mathcal{U}nif (0.2, 0.8)}{
#'       \sigma^2_m ~ Unif (0.2, 0.8)}
#'
#' For the variables that do not load on any factors (i.e., for which the
#' corresponding elements of \code{dedic} are equal to 0), it is specified that
#' \eqn{\alpha_m^\Delta = 0} and \eqn{\sigma^2_m = 1}.
#'
#' Covariance matrix of the latent factors:
#' \deqn{\Omega \sim \mathcal{I}nv-\mathcal{W}ishart(K+5, I_K)}{
#'       \Omega ~ Inv-Wishart(K+5, I_K)}
#' which is rescaled to be a correlation matrix if \code{R.corr = TRUE}:
#' \deqn{R = \Lambda^{-1/2} \Omega \Lambda^{-1/2}}{
#'       R = \Lambda^-1/2 \Omega \Lambda^-1/2}
#' \deqn{\Lambda = diag(\Omega)}
#'
#' Note that the distribution of the covariance matrix is truncated such that
#' all the off-diagonal elements of the implied correlation matrix \eqn{R} are
#' below \code{max.corr} in absolute value. The truncation is also applied if
#' the covariance matrix is used instead of the correlation matrix (i.e., if
#' \code{R.corr = FALSE}).
#'
#' The distributions and the corresponding default values used to simulate the
#' model parameters are specified as in the Monte Carlo study of CFSHP, see
#' section 4.1 (p.43).
#'
#' @return The function returns a data frame with \code{N} observations
#' simulated from the corresponding dedicated factor model.
#' The parameters used to generate the data are saved as attributes:
#' \code{dedic}, \code{alpha}, \code{sigma} and \code{R}.
#'
#' @author Rémi Piatek \email{remi.piatek@@econ.ku.dk}
#'
#' @references
#' G. Conti, S. Frühwirth-Schnatter, J.J. Heckman,
#' R. Piatek (2014): ``Bayesian Exploratory Factor Analysis'',
#' \emph{Journal of Econometrics}, 183(1), pages 31-57,
#' \url{http://dx.doi.org/10.1016/j.jeconom.2014.06.008}.
#'
#' @examples
#' # generate 1000 observations from model with 4 factors and 20 variables
#' # (5 variables loading on each factor)
#' dat <- simul.dedic.facmod(N = 1000, dedic = rep(1:4, each = 5))
#'
#' # generate data set with 5000 observations from the following model:
#' dedic <- rep(1:3, each = 4)        # 3 factors and 12 manifest variables
#' alpha <- rep(c(1, NA, NA, NA), 3)  # set first loading to 1 for each factor,
#'                                    #   sample remaining loadings from default
#' sigma <- rep(0.5, 12)              # idiosyncratic variances all set to 0.5
#' R <- toeplitz(c(1, .6, .3))        # Toeplitz matrix
#' dat <- simul.dedic.facmod(N = 5000, dedic, alpha, sigma, R)
#'
#' @export simul.dedic.facmod
#' @import checkmate
#' @importFrom stats rWishart rgamma rnorm runif cov2cor

simul.dedic.facmod <- function(N, dedic, alpha, sigma, R, R.corr = TRUE,
                               max.corr = 0.85, R.max.trial = 1000)
{

  # sanity checks for N and dedic, relabel appropriately
  assertInt(N, lower = 1)
  assertIntegerish(dedic, lower = 0, any.missing = FALSE, min.len = 3)
  dedic <- relabel.dedic(dedic)
  nmeas <- length(dedic)
  nfac <- max(dedic)

  # factor loading matrix
  if (missing(alpha)) {
    alpha <- rep(NA, nmeas)
  }
  assertNumeric(alpha, finite = TRUE, len = nmeas)
  nna <- sum(is.na(alpha))
  if (nna > 0) {
    alpha[is.na(alpha)] <- sample(c(-1, 1), nna, replace = TRUE) *
                             sqrt(runif(nna, min = 0.04, max = 0.64))
  }
  alpha.mat <- matrix(0, nmeas, nfac)
  for (i in 1:nmeas) {
    if (dedic[i] > 0)
      alpha.mat[i, dedic[i]] <- alpha[i]
  }

  # error terms
  if (missing(sigma)) sigma <- rep(NA, nmeas)
  assertNumeric(sigma, lower = 0.001, finite = TRUE, len = nmeas)
  nna <- sum(is.na(sigma))
  if (nna > 0) sigma[is.na(sigma)] <- runif(nna, min = 0.2, max = 0.8)
  sigma[dedic == 0] <- 1
  eps <- matrix(rnorm(N * nmeas), N, nmeas) %*% diag(sqrt(sigma))

  # latent factors
  assertLogical(R.corr, len = 1, any.missing = FALSE)
  assertNumber(max.corr, lower = 0, upper = 1)
  assertInt(R.max.trial, lower = 1, upper = 10^6)
  if (missing(R)) {
    for (i in 1:R.max.trial) {
      Omega <- rWishart(n = 1, df = nfac + 5, Sigma = diag(nfac))
      R <- cov2cor(Omega[, , 1])
      rho <- max(abs(R[lower.tri(R)]))
      if (rho <= max.corr)
        break
    }
    if (rho > max.corr)
      stop("R.max.trial exceeded, could not sample R")
    if (!R.corr)
      R <- Omega
  } else {
    if (!is.pos.semidefinite.matrix(R))
      stop("R is not positive semidefinite.")
    if (R.corr & !all(diag(R) == 1))
      stop("R is not a correlation matrix")
  }
  theta <- matrix(rnorm(N * nfac), N, nfac) %*% chol(R)

  # manifest variables
  Y <- theta %*% t(alpha.mat) + eps
  colnames(Y) <- paste0("Y", 1:nmeas)
  rownames(Y) <- 1:N

  # label parameters
  names(dedic) <- paste0("Y", 1:nmeas)
  names(alpha) <- paste0("alpha:", 1:nmeas)
  names(sigma) <- paste0("sigma:", 1:nmeas)
  rownames(R)  <- colnames(R) <- paste0("R:", 1:nfac)

  # return simulated data
  output <- as.data.frame(Y)
  attr(output, "dedic") <- dedic
  attr(output, "alpha") <- alpha
  attr(output, "sigma") <- sigma
  attr(output, "R")     <- R
  return(output)

}
