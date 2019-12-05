#'
#' Summarize 'befa' object
#'
#' Generic function summarizing the posterior results of a 'befa' object.
#' Optional arguments can be specified to customize the summary.
#'
#' @param object
#'        Object of class 'befa'.
#' @param ...
#'        The following extra arguments can be specified:
#' \itemize{
#'   \item \code{what}: How to summarize the posterior distribution?
#'         \itemize{
#'           \item \code{what = 'maxp'} (default): Only factor loadings with
#'                 highest posterior probability of being different from zero or
#'                 discarded from the model (if \code{dedic = 0}) are
#'                 summarized.
#'           \item \code{what = 'all'}: All factor loadings with corresponding
#'                 posterior probability to be allocated to a given factor (or
#'                 to be discarded from the model) larger than \code{min.prob}
#'                 are summarized.
#'           \item \code{what = 'hppm'}: Highest posterior probability models
#'                 with probability larger than \code{min.prob} are summarized.
#'         }
#'   \item \code{byfac}: Sort factor loadings by factors if \code{TRUE},
#'         otherwise by manifest variables if \code{FALSE} (default).
#'   \item \code{hpd.prob}: Probability used to compute the highest posterior
#'         density intervals of the posterior distribution of the model
#'         parameters (default: 0.95).
#'   \item \code{min.prob}: If \code{what = 'all'}, only factor loadings with
#'         posterior probability of being dedicated to a given factor (or
#'         discarded from the model) larger than this value are displayed.
#'         If \code{what = 'hppm'}, only highest posterior probability models
#'         with probability larger than this value are displayed. (default:
#'         0.20)
#' }
#'
#' @details
#' This function summarizes the posterior distribution of the parameters.
#' The algorithm may visit different configurations of the indicator matrix
#' \eqn{\Delta} during sampling, where the manifest variables are allocated to
#' different latent factors. When the posterior distribution of the factor
#' loadings is summarized separately for each manifest variable
#' (\code{what = 'maxp'} or \code{what = 'all'}), the function provides the
#' latent factor each manifest variable is allocated to (\code{dedic}), and the
#' corresponding posterior probability (\code{prob}). If \code{dedic = 0}, then
#' \code{prob} corresponds to the posterior probability that the manifest
#' variable is discarded. Discarded variables are listed last if
#' \code{byfac = TRUE}. Low probability cases can be discarded by setting
#' \code{min.prob} appropriately (default is 0.20).
#'
#' Idiosyncratic variances, factor correlation matrix and regression
#' coefficients (if any) are summarized across all MCMC iterations if
#' \code{what = 'all'} or \code{what = 'maxp'}, and within each HPP model if
#' \code{what = 'hppm'}.
#'
#' \strong{Highest posterior probability model.}
#' The HPP model is the model with a given allocation of the measurements to the
#' latent factors (i.e., a given indicator matrix \eqn{\Delta}) that is visited
#' most often by the algorithm.
#'
#' When specifying \code{what = 'hppm'}, the function sorts the models according
#' to the posterior frequencies of their indicator matrices in decreasing order.
#' Therefore, the first model returned (labeled 'm1') corresponds to the HPP
#' model.
#' Low probability models can be discarded by setting \code{min.prob}
#' appropriately(default is 0.20, implying that only models with a posterior
#' probability larger than 0.20 are displayed).
#'
#' HPP models can only be found if identification with respect to column
#' switching has been restored \emph{a posteriori}. An error message is returned
#' if this is not the case.
#'
#' @return If called directly, the summary is formatted and displayed on the
#' standard output. Otherwise if saved in an object, a list of the following
#' elements is returned:
#' \itemize{
#'  \item \code{MHacc}: Metropolis-Hastings acceptance rate.
#'  \item \code{alpha}: Data frame (or list of data frames if
#'        \code{what = 'hppm'}) containing posterior summary statistics for the
#'        factor loadings.
#'  \item \code{sigma}: Data frame (or list of matrices if \code{what = 'hppm'})
#'        containing posterior summary statistics for the idiosyncratic
#'        variances.
#'  \item \code{R}: Data frame (or list of data frames if \code{what = 'hppm'})
#'        containing posterior summary statistics for the factor correlations.
#'  \item \code{beta}: Data frame (or list of data frames if
#'        \code{what = 'hppm'}) containing posterior summary statistics for the
#'        regression coefficients (if any).
#'  \item \code{nfac} (only if \code{what = 'maxp'} or \code{what = 'all'}):
#'        Table of posterior frequencies of numbers of factors.
#'  \item \code{hppm} (only if \code{what = 'hppm'}): List of the following
#'        elements summarizing the different HPP models, sorted in decreasing
#'        order of their posterior probabilities:
#'        \itemize{
#'          \item \code{prob}: Vector of posterior probabilities.
#'          \item \code{nfac}: Vector of numbers of factors.
#'          \item \code{dedic}: Data frame of factor indicators.
#'        }
#' }
#' Data frames of posterior summary statistics include the means (\code{mean}),
#' standard deviations (\code{sd}) and highest posterior density intervals
#' (\code{hpd.lo} and \code{hpd.up}, for the probability specified in
#' \code{hpd.prob}) of the corresponding parameters.
#'
#' For the factor loadings, the matrix may also include a column labeled
#' '\code{dedic}' indicating to which factors the corresponding manifest
#' variables are dedicated (a zero value means that the manifest variable does
#' not load on any factor), as well as a column labeled '\code{prob}' showing
#' the corresponding posterior probabilities that the manifest variables load on
#' these factors.
#'
#' Summary results are returned as lists of data frames for HPP models, where
#' the elements of the list are labeled as '\code{m1}, '\code{m2}', etc.
#'
#' @author Rémi Piatek \email{remi.piatek@@econ.ku.dk}
#'
#' @seealso \code{\link{plot.befa}} to plot posterior results.
#'
#' @examples
#' set.seed(6)
#'
#' # generate fake data with 15 manifest variables and 3 factors
#' Y <- simul.dedic.facmod(N = 100, dedic = rep(1:3, each = 5))
#'
#' # run MCMC sampler and post process output
#' # notice: 1000 MCMC iterations for illustration purposes only,
#' #  increase this number to obtain reliable posterior results!
#' mcmc <- befa(Y, Kmax = 5, iter = 1000)
#' mcmc <- post.column.switch(mcmc)
#' mcmc <- post.sign.switch(mcmc)
#'
#' # summarize posterior results
#' summary(mcmc)
#'
#' # summarize highest posterior probability (HPP) model
#' hppm.sum <- summary(mcmc, what = 'hppm')
#'
#' # print summary with 6-digit precision
#' print(hppm.sum, digits = 6)
#'
#' # extract posterior means of the factor loadings in HPP model
#' alpha.mean <- hppm.sum$alpha$m1$mean
#' print(alpha.mean)
#'
#' @examples \dontshow{
#' summary(mcmc, what = 'maxp', byfac = TRUE)
#' summary(mcmc, what = 'all')
#' summary(mcmc, what = 'all', byfac = TRUE)
#' summary(mcmc, what = 'all', min.prob = 0)
#' summary(mcmc, what = 'all', min.prob = 0, byfac = TRUE)
#' summary(mcmc, what = 'hppm', byfac = TRUE)
#' summary(mcmc, what = 'hppm', min.prob = 0)
#' summary(mcmc, what = 'hppm', min.prob = 0, byfac = TRUE)
#' }
#'
#' @export
#' @import checkmate
#' @importFrom stats sd
#' @importFrom coda as.mcmc HPDinterval
#' @importFrom plyr count

summary.befa <- function(object, ...)
{

  if (class(object) != 'befa')
    stop('object passed to print.befa should be of class befa')

  # extra arguments
  args <- list(...)
  min.prob <- ifelse (is.null(args$min.prob), .20, args$min.prob)
  hpd.prob <- ifelse (is.null(args$hpd.prob), .95, args$hpd.prob)
  byfac    <- ifelse (is.null(args$byfac), FALSE, args$byfac)
  assertNumber(min.prob, lower = 0, upper = 1)
  assertNumber(hpd.prob, lower = 0, upper = 1)
  assertFlag(byfac)
  what <- match.arg(args$what, choices = c('maxp', 'all', 'hppm'))

  # container for summary object
  output <- list()

  ### column-switching and/or sign-switching problems
  if (!attr(object, "post.column.switch")) {
    warning("MCMC output not processed by function 'post.column.switch'. ",
            "Posterior results for factor loadings and factor correlations ",
            "may not be interpretable! Make sure column-switching problem has ",
            "been fixed before summarizing.\n",
            immediate. = TRUE)
    if (what == 'hppm')
      stop("Highest posterior probability model (hppm) can only be searched ",
           " for if columns have be reordered a posteriori. Use function ",
           "'post.column.switch' first.\n")
  }

  if (!attr(object, "post.sign.switch")) {
    warning("MCMC output not processed by function 'post.sign.switch'. ",
            "Posterior results for factor loadings and factor correlations ",
            "may not be interpretable! Make sure sign-switching problem has ",
            "been fixed before summarizing.\n",
            immediate. = TRUE)
  }

  ### Metropolis-Hastings acceptance rate
  output$MHacc <- mean(object$MHacc)
  if (output$MHacc < .2)
    warning('Metropolis-Hastings acceptance rate is low (< 0.20). ',
            'Check convergence and mixing before summarizing.',
            immediate. = TRUE)

  ### function summarizing draws from posterior
  iter <- attr(object, 'iter')
  summarize.mcmc <- function(z, prob = FALSE) {
    if (length(z) == 1) {
      sd  <- 0
      hpd <- c(NA, NA)
    } else {
      sd  <- sd(z)
      hpd <- coda::HPDinterval(coda::as.mcmc(z), prob = hpd.prob)
    }
    res <- c(mean = mean(z), sd = sd, hpd.lo = hpd[1], hpd.up = hpd[2])
    if (prob) res <- c(prob = length(z)/iter, res)
    return(res)
  }

  ### function sorting manifest variables according to the factors they load on
  ### (variables loading on no factors are listed last)
  sort.byfac <- function(a) {
    a[a[, 'dedic'] == 0, 'dedic'] <- NA             # make 0 -> NA
    a <- a[order(a[, 'dedic'], na.last = TRUE), ]   # sort, putting NAs last
    a[is.na(a[, 'dedic']), 'dedic'] <- 0            # make NA -> 0
    return(a)
  }

  ### highest posterior probability models
  if (what == 'hppm') {

    alpha <- object$alpha
    dedic <- object$dedic

    dedic.tab <- plyr::count(as.data.frame(dedic))
    dedic.tab <- dedic.tab[order(dedic.tab$freq, decreasing = TRUE),]
    if (dedic.tab$freq[1] < iter*min.prob)
      stop('Probability of HPP model lower than min.prob. ',
           'Try to decrease the value of min.prob')
    dedic.tab  <- dedic.tab[dedic.tab$freq >= iter*min.prob,]

    hppm.freq  <- dedic.tab$freq/iter
    hppm.dedic <- dedic.tab[-length(dedic.tab)]
    hppm.dedic <- as.data.frame(t(hppm.dedic))
    colnames(hppm.dedic) <- paste0('m', 1:ncol(hppm.dedic))
    hppm.nfac  <- apply(hppm.dedic, 2, count.unique.nonzero)

    output$hppm <- list(prob  = hppm.freq,
                        nfac  = hppm.nfac,
                        dedic = hppm.dedic)

    output$alpha <- list()
    output$sigma <- list()
    output$R     <- list()
    if (!is.null(object$beta)) output$beta <- list()

    for (i in seq_along(hppm.freq)) {

      if (hppm.freq[i] < min.prob) break

      hppm.id <- apply(t(dedic) == hppm.dedic[,i], 2, all)

      a <- apply(alpha[hppm.id, , drop = FALSE], 2, summarize.mcmc)
      a <- cbind(dedic = hppm.dedic[,i], t(a))
      if (byfac) a <- sort.byfac(a)
      rownames(a) <- colnames(alpha)

      lab <- paste0('m', i)
      output$alpha[[lab]] <- a
      output$sigma[[lab]] <- t(apply(object$sigma[hppm.id, , drop = FALSE], 2, 
                                     summarize.mcmc))
      output$R[[lab]]     <- t(apply(object$R[hppm.id, , drop = FALSE], 2, 
                                     summarize.mcmc))
      if (!is.null(object$beta))
        output$beta[[lab]] <- t(apply(object$beta[hppm.id, , drop = FALSE], 2, 
                                      summarize.mcmc))

    }

    if (length(output$alpha) == 0)
      stop("No hpp model found! Try to decrease value of 'min.prob'.")

    # factor loadings: insert NA values when not dedicated to any factor
    for (i in 1:length(output$alpha))
      for (j in 1:nrow(output$alpha[[i]])) {
        cat(i, j, "\n")
        if (output$alpha[[i]][j, 'dedic'] == 0)
          output$alpha[[i]][j, c('mean', 'sd', 'hpd.lo', 'hpd.up')] <- NA
      }

  ### general case (not HPP models)
  } else {

    alpha <- mat2list(object$alpha)
    dedic <- mat2list(object$dedic)

    # posterior number of latent factors
    iter <- nrow(object$dedic)
    nfac <- table(object$nfac)/iter
    output$nfac <- nfac[nfac >= .05]   # discard if prob < 0.05

    # summarize posterior results for loadings
    if (what == 'all') {

      bysum <- function(x, d) by(x, d, summarize.mcmc, prob = TRUE)
      a <- mapply(bysum, alpha, dedic, SIMPLIFY = FALSE)
      a <- lapply(a, function(x) do.call(rbind, x))
      a <- lapply(a, function(x) cbind(dedic = as.integer(rownames(x)), x))
      for(i in 1:length(a)) rownames(a[[i]]) <- rep(names(a)[i], nrow(a[[i]]))
      a <- do.call(rbind, a)
      a <- a[a[, 'prob'] >= min.prob,]   # drop rows for which prob < min.prob

    } else if (what == 'maxp') {

      get.max <- function(x) as.integer(names(which.max(table(x))))
      dhp     <- sapply(dedic, get.max)
      hpsum   <- function(x, d, dhp) summarize.mcmc(x[d == dhp], prob = TRUE)
      a <- t(mapply(hpsum, alpha, dedic, dhp))
      a <- cbind(dedic = dhp, a)
      rownames(a) <- names(alpha)

    }
    if (byfac) a <- sort.byfac(a)
    output$alpha <- a

    # summarize posterior results for remaining parameters
    output$sigma <- t(apply(object$sigma, 2, summarize.mcmc))
    output$R     <- t(apply(object$R, 2, summarize.mcmc))
    if (!is.null(object$beta))
      output$beta <- t(apply(object$beta, 2, summarize.mcmc))

    # factor loadings: insert NA values when not dedicated to any factor
    for (i in 1:nrow(output$alpha))
      if (output$alpha[i, 'dedic'] == 0)
        output$alpha[i, c('mean', 'sd', 'hpd.lo', 'hpd.up')] <- NA

  }

  # convert parameter summaries to data frames,
  # rename duplicate row names if necessary (esp. for factor loadings)
  mdf <- function(x) as.data.frame(x, row.names = make.unique(rownames(x)))
  if (what == 'hppm') {
    output$alpha <- lapply(output$alpha, mdf)
    output$sigma <- lapply(output$sigma, mdf)
    output$R     <- lapply(output$R,     mdf)
    if (!is.null(output$beta)) output$beta <- lapply(output$beta, mdf)
  } else {
    output$alpha <- mdf(output$alpha)
    output$sigma <- mdf(output$sigma)
    output$R     <- mdf(output$R)
    if (!is.null(output$beta)) output$beta <- mdf(output$beta)
  }

  ### save attributes and return object
  attr(output, 'Kmax')     <- attr(object, 'Kmax')
  attr(output, 'Nid')      <- attr(object, 'Nid')
  attr(output, 'iter')     <- attr(object, 'iter')
  attr(output, 'burnin')   <- attr(object, 'burnin')
  attr(output, 'what')     <- what
  attr(output, 'min.prob') <- min.prob
  attr(output, 'hpd.prob') <- hpd.prob
  attr(output, 'byfac')    <- byfac

  class(output) <- "summary.befa"
  return(output)

}
