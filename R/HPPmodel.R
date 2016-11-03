#'
#' Highest posterior probability model
#'
#' This function finds the highest posterior probability (HPP) model, which
#' corresponds to the indicator matrix most often visited by the sampler across
#' MCMC iterations.
#'
#' @param mcmc
#'        Object of class '\code{befa}'.
#'
#' @details The HPP model can only be found if identification with respect to
#' column switching has been restored \emph{a posteriori}. An error message is
#' returned if this is not the case.
#'
#'  A warning is issued if the HPP model found is not unique.
#'
#' @return The function returns an object of class '\code{hppm.befa}'
#' containing:
#' \itemize{
#'   \item \code{draws}: Posterior draws of model parameters corresponding to
#'   the HPP model.
#'   \item \code{dedic}: Vector of indicators of the the HPP model.
#'   \item \code{prob}: Posterior probability of the HPP model.
#' }
#'
#' @author Rémi Piatek \email{remi.piatek@@econ.ku.dk}
#'
#' @seealso \code{\link{post.column.switch}}
#'
#' @examples
#' set.seed(6)
#' Y <- simul.dedic.facmod(N = 200, dedic = rep(1:3, each = 5))
#' mcmc <- befa(Y, Kmax = 5, iter = 1000)
#' mcmc <- post.column.switch(mcmc)
#' hppm <- HPPmodel(mcmc)
#' hppm$prob
#' hppm$dedic
#'
#' @export HPPmodel
#' @import checkmate
#' @importFrom plyr count

HPPmodel <- function(mcmc)
{

  assertClass(mcmc, classes = "befa")
  if (!attr(mcmc, "post.column.switch")) {
    stop("post.column.switch should be run first.")
  }

  dedic <- mcmc$dedic
  nmeas <- ncol(dedic)
  niter <- nrow(dedic)

  tab <- plyr::count(as.data.frame(dedic))
  freq <- tab[, nmeas + 1]
  hppm <- which.max(freq)
  hppm.freq <- freq[hppm]/niter
  hppm.dedic <- tab[hppm, 1:nmeas]

  if (sum(freq == hppm.freq) > 1) {
    warning("HPP model not unique, first one found returned.")
  }

  hd <- matrix(hppm.dedic, niter, nmeas, byrow = TRUE)
  id <- apply(dedic == hd, 1, all)

  output <- list(draws = mcmc$draws[id, ], dedic = hppm.dedic, prob = hppm.freq)
  attr(output, "title") <- "HPP model for BEFA posterior sample"
  attr(output, "post.sign.switch") <- attr(mcmc, "post.sign.switch")
  class(output) <- "hppm.befa"
  return(output)

}
