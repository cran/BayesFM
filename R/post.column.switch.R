#'
#' Perform column switchting on posterior MCMC sample
#'
#' This function reorders the columns of the factor loading matrix for each MCMC
#' draw, as well as the rows and columns of the correlation matrix of the
#' factors, to restore the identification of the model \emph{a posteriori} with
#' respect to the column switching problem.
#'
#' @param mcmc
#'        Object of class '\code{befa}'.
#'
#' @details The reordering of the columns of the factor loading matrix is based
#' on the top elements of the columns (i.e., the first row containing a nonzero
#' factor loading in each nonzero column of \eqn{\alpha}, starting from the top
#' of the matrix). At each MCMC iteration, the nonzero columns of \eqn{\alpha}
#' are reordered such that the top elements appear in increasing order.
#' The rows and columns of the correlation matrix \eqn{R} of the factors are
#' switched accordingly. See section 4.3 of CFSHP (p.42) for more details.
#'
#' @return Same '\code{befa}' object as the one passed to the function, where
#' the indicators in the matrix \code{dedic}, as well as the rows and columns of
#' the correlation matrix of the factors saved in \code{draws}, have been
#' switched appropriately to restore the identification of the factor model with
#' respect to column switching.
#'
#' @author Rémi Piatek \email{remi.piatek@@econ.ku.dk}
#'
#' @references
#' G. Conti, S. Frühwirth-Schnatter, J.J. Heckman,
#' R. Piatek (2014): ``Bayesian Exploratory Factor Analysis'',
#' \emph{Journal of Econometrics}, 183(1), pages 31-57,
#' \url{http://dx.doi.org/10.1016/j.jeconom.2014.06.008}.
#'
#' @seealso \code{\link{post.sign.switch}} to restore identification a
#' posteriori with respect to the sign switching problem.
#'
#' @examples
#' set.seed(6)
#' Y <- simul.dedic.facmod(N = 200, dedic = rep(1:3, each = 5))
#' mcmc <- befa(Y, Kmax = 5, iter = 1000)
#' mcmc <- post.column.switch(mcmc)
#'
#' @export post.column.switch
#' @import checkmate

post.column.switch <- function(mcmc)
{

  assertClass(mcmc, classes = "befa")
  if (attr(mcmc, "post.column.switch")) {
    warning("column switching already performed, nothing done.")
    return(mcmc)
  }

  # issue warning if M-H acceptance rate too low
  if (mean(mcmc$MHacc) < 0.2) {
    warning(paste("M-H acceptance rate of sampler is low (< 0.20).",
                  "Check convergence and mixing!"))
  }

  Kmax   <- attr(mcmc, "Kmax")
  nvar   <- ncol(mcmc$dedic)
  iter   <- nrow(mcmc$dedic)
  R.npar <- Kmax * (Kmax - 1)/2

  # index matrix used to reconstruct matrix from lower triangular elements
  R.mat <- diag(Kmax) * (R.npar + 1)
  R.mat[lower.tri(R.mat)] <- 1:R.npar
  R.mat[upper.tri(R.mat)] <- t(R.mat)[upper.tri(R.mat)]

  v <- 1:Kmax
  for (i in 1:iter) {
    d <- mcmc$dedic[i, ]

    # relabel indicators
    mcmc$dedic[i, ] <- relabel.dedic(d)

    # reorder rows and columns of correlation matrix
    u <- unique(d[d != 0])
    r <- c(u, v[!v %in% u])
    R <- c(mcmc$R[i, ], 1)
    R <- matrix(R[R.mat], nrow = Kmax)
    R <- R[r, r]
    mcmc$R[i, ] <- R[lower.tri(R)]
  }

  attr(mcmc, "post.column.switch") <- TRUE
  return(mcmc)

}
