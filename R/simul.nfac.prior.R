#'
#' Simulate prior distribution of number of latent factors
#'
#' This function produces a sample from the prior distribution of the number of
#' latent factors. It depends on the prior parameters used for the distribution
#' of the indicators, on the size of the model (number of manifest variables
#' and maximum number of latent factors), and on the identification restriction
#' (minimum number of manifest variables dedicated to each factor).
#'
#' @inheritParams befa
#' @param nvar Number of manifest variables.
#' @param Kmax Maximum number of latent factors.
#' @param kappa Concentration parameter of the Dirichlet prior distribution on
#'        the indicators.
#' @param nrep  Number of Monte Carlo replications.
#'
#' @details This function simulates the prior distribution of the number of
#' latent factors for models that fulfill the identification restriction
#' restriction that at least \code{Nid} manifest variables (or no variables) are
#' loading on each latent factor. Several (scalar) parameters \code{kappa} can
#' be passed to the function to simulate the prior for different prior parameter
#' values and compare the results.
#'
#' An accept/reject sampling scheme is used: a vector of probabilities is drawn
#' from a Dirichlet distribution with concentration parameter \code{kappa}, and
#' the \code{nvar} manifest variables are randomly allocated to the \code{Kmax}
#' latent factors. If each latent factor has at least \code{Nid} dedicated
#' variables or no variables at all, the identification requirement is fulfilled
#' and the draw is accepted. The number of factors loaded by at least \code{Nid}
#' manifest variables is returned as a draw from the prior distribution.
#'
#' Note that this function does not use the two-level prior distribution
#' implemented in CFSHP, where manifest variables can be discarded from the
#' model according to a given probability. Therefore, this function only help
#' understand the prior distribution conditional on all the manifest variables
#' being included into the model.
#'
#' @return A list of length equal to the number of parameters specified in
#' \code{kappa} is returned, where each element of the list contains:
#' \itemize{
#'   \item \code{nfac}: Vector of integers of length equal to the number of
#'   accepted draws.
#'   \item \code{acc}: Acceptance rate of the accept/reject sampling scheme.
#' }
#'
#' @author Rémi Piatek \email{remi.piatek@@econ.ku.dk}
#'
#' @references G. Conti, S. Frühwirth-Schnatter, J.J. Heckman, R. Piatek (2014):
#' ``Bayesian Exploratory Factor Analysis'', \emph{Journal of Econometrics},
#' 183(1), pages 31-57, \url{http://dx.doi.org/10.1016/j.jeconom.2014.06.008}.
#'
#' @examples
#' # replicate first row of table 2 in CFSHP (p.44)
#' # note: use larger number of replications nrep to improve accuracy
#' prior.nfac <- simul.nfac.prior(nvar = 15, Kmax = 5, kappa = c(.3, .7, 1))
#' summary(prior.nfac)
#' plot(prior.nfac)
#'
#' @export simul.nfac.prior
#' @import checkmate
#' @useDynLib BayesFM, .registration = TRUE, .fixes = "F_"

simul.nfac.prior <- function(nvar, Kmax, Nid = 3, kappa = 1/Kmax, nrep = 10^6)
{

  # sanity checks
  checkArgs <- makeAssertCollection()
  assertInt(nvar, lower = 1, add = checkArgs)
  assertInt(Kmax,  lower = 1, upper = nvar, add = checkArgs)
  assertInt(Nid,   lower = 1, upper = floor(nvar/Kmax), add = checkArgs)
  assertNumeric(kappa, lower = 10^-7, finite = TRUE, any.missing = FALSE,
                min.len = 1, unique = TRUE, add = checkArgs)
  assertCount(nrep, positive = TRUE, add = checkArgs)
  reportAssertions(checkArgs)

  out <- list()
  for (kap in kappa) {

    seed <- round(runif(1) * 10^9)
    sim  <- .Fortran(F_simnfacprior,
                     as.integer (nvar),
                     as.integer (Kmax),
                     as.integer (Nid),
                     as.double  (rep(kap, Kmax)),
                     as.integer (nrep),
                     as.integer (seed),
                     nfac    = integer(nrep),
                     restrid = logical(nrep),
                     PACKAGE = 'BayesFM')

    lab <- paste('kappa =', kap)
    out[[lab]] <- list(nfac = sim$nfac[sim$restrid],
                       acc  = sum(sim$restrid) / nrep)

  }

  attr(out, 'call')  <- match.call()
  attr(out, 'nvar')  <- nvar
  attr(out, 'Kmax')  <- Kmax
  attr(out, 'Nid')   <- Nid
  attr(out, 'kappa') <- kappa
  attr(out, 'nrep')  <- nrep
  class(out) <- 'simul.nfac.prior'
  return(out)

}
