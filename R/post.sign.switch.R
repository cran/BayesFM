#'
#' Perform sign switchting on posterior MCMC sample
#'
#' This function performs a sign switch on the MCMC draws to restore the
#' consistency of the signs of the factors loadings and of the correlations of
#' the latent factors \emph{a posteriori}.
#'
#' @param mcmc
#'        Object of class '\code{befa}'.
#' @param benchmark
#'        Vector of integers of length equal to the maximum number of latent
#'        factors. Each element indicates which factor loading is used as a
#'        benchmark for the sign switch. If \code{NULL}, the factor loadings
#'        with the highest posterior probabilities of being different from zero
#'        in each column of the factor loading matrix are used as benchmarks.
#' @param benchmark.threshold
#'        Minimum posterior probability for a factor loading to be considered
#'        as a benchmark.
#'
#' @details The signs of the factor loadings, as well as of the corresponding
#' correlations of the latent factors, are switched for each MCMC iteration such
#' that the factor loadings defined as \code{benchmark}s are positive. The sign
#' switch can only be performed if \code{\link{post.column.switch}} has been run
#' before. See section 4.3 (p.42) of CFSHP for more details.
#'
#' If a latent factor has no benchmarks, or if its benchmark is equal to zero at
#' some MCMC iteration, then no sign switch is performed on the corresponding
#' loadings and correlations for this particular factor or MCMC iteration.
#'
#' Note that in complicated models where the sampler visits several models with
#' different numbers of latent factors, it may not be relevant to use the
#' default value of \code{benchmark}, as the posterior probabilities that the
#' factor loadings are different from zero would be computed across models.
#' Instead, the user might consider finding the highest posterior probability
#' model first, and use its top elements in each column of the factor loading
#' matrix as benchmarks to perform the sign switch.
#'
#' @return This function returns the same '\code{befa}' object, where the signs
#' of the factor loadings and of the factor correlations have been switched
#' appropriately to restore the identification of the factor model with respect
#' to sign switching.
#'
#' @author Rémi Piatek \email{remi.piatek@@econ.ku.dk}
#'
#' @references
#' G. Conti, S. Frühwirth-Schnatter, J.J. Heckman,
#' R. Piatek (2014): ``Bayesian Exploratory Factor Analysis'',
#' \emph{Journal of Econometrics}, 183(1), pages 31-57,
#' \url{http://dx.doi.org/10.1016/j.jeconom.2014.06.008}.
#'
#' @seealso \code{\link{post.column.switch}} for column switching of the factor
#' loading matrix and of the correlation matrix of the latent factors to restore
#' identification \emph{a posteriori}.
#'
#' @examples
#' set.seed(6)
#' Y <- simul.dedic.facmod(N = 200, dedic = rep(1:3, each = 5))
#' mcmc <- befa(Y, Kmax = 5, iter = 1000)
#' mcmc <- post.column.switch(mcmc)
#'
#' # factor loadings corresponding to variables 1, 6, 11, 12 and 13 are
#' # used as benchmarks:
#' mcmc1 <- post.sign.switch(mcmc, benchmark = c(1, 6, 11, 12, 13))
#'
#' # factor loadings with the highest posterior probability of being different
#' # from zero in each column are used as benchmark:
#' mcmc2 <- post.sign.switch(mcmc)
#'
#' @export post.sign.switch
#' @import checkmate

post.sign.switch <- function(mcmc, benchmark = NULL, benchmark.threshold = 0.5)
{

  assertClass(mcmc, classes = "befa")
  if (attr(mcmc, "post.sign.switch")) {
    warning("sign switching already performed, nothing done.")
    return(mcmc)
  }
  if (!attr(mcmc, "post.column.switch")) {
    stop("post.column.switch should be run first.")
  }
  assertNumber(benchmark.threshold, lower = 0, upper = 1)

  Kmax   <- attr(mcmc, "Kmax")
  nvar   <- ncol(mcmc$dedic)
  iter   <- nrow(mcmc$dedic)
  R.npar <- Kmax * (Kmax - 1)/2

  # factor loadings used as benchmarks
  if (is.null(benchmark)) {
    alpha.post.prob <- matrix(0, Kmax, nvar)
    for (k in 1:Kmax) {
      alpha.post.prob[k, ] <- colMeans(mcmc$dedic == k)
    }
    benchmark <- max.col(alpha.post.prob, ties.method = "first")
    for (k in 1:Kmax) {
      if (alpha.post.prob[k, benchmark[k]] < benchmark.threshold) {
        benchmark[k] <- NA  # not a benchmark if post prob < threshold
      }
    }
  } else {
    assertIntegerish(benchmark, lower = 0, upper = nvar, any.missing = FALSE,
                     len = Kmax)
  }

  # index matrix used to reconstruct matrix from lower triangular elements
  R.mat <- diag(Kmax) * (R.npar + 1)
  R.mat[lower.tri(R.mat)] <- 1:R.npar
  R.mat[upper.tri(R.mat)] <- t(R.mat)[upper.tri(R.mat)]

  for (i in 1:iter) {
    # switch signs of factor loadings
    dedic <- mcmc$dedic[i, ]
    switch <- sign(mcmc$alpha[i, benchmark])
    switch[is.na(switch)] <- 1  # no sign switch for factors with no benchmark
    switch.meas <- rep(0, nvar)
    for (j in 1:nvar) {
      if (dedic[j] == 0) next
      switch.meas[j] <- switch[dedic[j]]
    }
    mcmc$alpha[i, ] <- mcmc$alpha[i, ] * switch.meas

    # switch signs of rows and columns of correlation matrix
    R <- c(mcmc$R[i, ], 1)
    R <- matrix(R[R.mat], nrow = Kmax)
    R <- diag(switch) %*% R %*% diag(switch)
    mcmc$R[i, ] <- R[lower.tri(R)]
  }

  attr(mcmc, "post.sign.switch") <- TRUE
  return(mcmc)

}
