#'
#' Simulate prior distribution of factor correlation matrix
#'
#' This function produces a sample of correlation matrices drawn from their
#' prior distribution induced in the identified version of the factor model,
#' given the prior distribution specified on the corresponding covariance
#' matrices of the factors in the expanded model.
#'
#' @inheritParams befa
#' @param Kmax Maximum number of latent factors.
#' @param nrep Number of Monte Carlo replications.
#'
#' @details Covariance matrices are sampled from the prior distribution in the
#' expanded model, and transformed to produce the corresponding correlation
#' matrices. See section 2.3.5 of CFSHP (p.36-37), as well as the details of
#' the function \code{\link{befa}}.
#'
#' To compare several prior specifications, different values of the parameters
#' \code{nu0} and \code{S0} can be specified. The function then simulates for
#' each pair of these parameters. \code{nu0} and \code{S0} should therefore be
#' scalars or vectors of same length.
#'
#' @return A list of length equal to the number of pairs of parameters
#' \code{nu0} and \code{S0}, where each element of the list is an array of
#' dimension (\code{Kmax}, \code{Kmax}, \code{nrep}) that contains the
#' correlation matrices of the latent factors drawn from the prior.
#'
#' @author Rémi Piatek \email{remi.piatek@@gmail.com}
#'
#' @references G. Conti, S. Frühwirth-Schnatter, J.J. Heckman, R. Piatek (2014):
#' ``Bayesian Exploratory Factor Analysis'', \emph{Journal of Econometrics},
#' 183(1), pages 31-57, \doi{10.1016/j.jeconom.2014.06.008}.
#'
#' @examples
#' # partial reproduction of figure 1 in CFSHP (p.38)
#' # note: use larger number of replications nrep to increase smoothness
#' Kmax <- 10
#' Rsim <- simul.R.prior(Kmax, nu0 = Kmax + c(1, 2, 5), S0 = .5, nrep = 1000)
#' summary(Rsim)
#' plot(Rsim)
#'
#' @export simul.R.prior
#' @import checkmate
#' @importFrom stats rWishart rgamma

simul.R.prior <- function(Kmax, nu0 = Kmax + 1, S0 = 1, HW.prior = TRUE,
                          nrep = 10^5, verbose = TRUE) {

  # sanity checks
  checkArgs <- makeAssertCollection()
  assertInt(Kmax, lower = 1, add = checkArgs)
  assertNumeric(nu0, lower = Kmax, finite = TRUE, any.missing = FALSE,
                min.len = 1, add = checkArgs)
  assertNumeric(S0, lower = 0.000001, finite = TRUE, any.missing = FALSE,
                min.len = 1, add = checkArgs)
  if (length(nu0) == 1) nu0 <- rep(nu0, length(S0))
  if (length(S0) == 1)  S0  <- rep(S0,  length(nu0))
  if (length(nu0) != length(S0))
    checkArgs$push('nu0 and S0 must be of same length')
  assertFlag(HW.prior, add = checkArgs)
  assertCount(nrep, positive = TRUE, add = checkArgs)
  assertFlag(verbose, add = checkArgs)
  reportAssertions(checkArgs)

  Rmat <- list()
  par  <- unique(cbind(nu0, S0))
  for (i in 1:nrow(par))
  {
    nu0i <- par[i, 1]
    S0i  <- par[i, 2]
    cat('simulating prior for nu0 =', nu0i, 'and S0 =', S0i, '\n')
    omg <- rWishart(nrep, df = nu0i, Sigma = diag(Kmax))
    if (HW.prior) {
      nus <- nu0i - Kmax + 1
      Sk  <- replicate(Kmax, rgamma(nrep, shape = .5, rate = 1/(2 * nus * S0i)))
      Sk  <- sqrt(1/Sk)
      omg <- lapply(seq(nrep), function(x)
               diag(Sk[x,]) %*% solve(omg[,,x]) %*% diag(Sk[x,]))
    } else {
      S0d <- sqrt(1/S0i) * diag(Kmax)
      omg <- lapply(seq(nrep), function(x) S0d %*% solve(omg[,,x]) %*% S0d)
    }
    omg <- lapply(omg, cov2cor)
    lab <- paste0('nu0 = ', nu0i, ', S0 = ', S0i)
    Rmat[[lab]] <- simplify2array(omg)
  }

  attr(Rmat, 'call') <- match.call()
  attr(Rmat, 'Kmax') <- Kmax
  attr(Rmat, 'nrep') <- nrep
  class(Rmat) <- 'simul.R.prior'
  return(Rmat)

}
