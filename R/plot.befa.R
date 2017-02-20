#'
#' Plot object of class 'befa'
#'
#' This function makes different plots that are useful to assess the posterior
#' results: a trace plot of the number of latent factors (also showing
#' Metropolis-Hastings acceptance across MCMC replications), a histogram
#' of the posterior probabilities of the number of factors, heatmaps for the
#' inficator probabilities, the factor loading matrix, and the correlation
#' matrix of the latent factors.
#'
#' @inheritParams summary.befa
#' @param x Object of class 'befa'.
#'
#' @details This function makes graphs based on the summary results returned by
#' \code{\link{summary.befa}}. It therefore accepts the same optional arguments
#' as this function.
#'
#' @author Rémi Piatek \email{remi.piatek@@econ.ku.dk}
#'
#' @seealso \code{\link{summary.befa}} to summarize posterior results.
#'
#' @examples
#' set.seed(6)
#'
#' # generate fake data with 15 manifest variables and 3 factors
#' Y <- simul.dedic.facmod(N = 200, dedic = rep(1:3, each = 5))
#'
#' # run MCMC sampler and post process output
#' # notice: 1000 MCMC iterations for illustration purposes only,
#' #  increase this number to obtain reliable posterior results!
#' mcmc <- befa(Y, Kmax = 5, iter = 1000)
#' mcmc <- post.column.switch(mcmc)
#' mcmc <- post.sign.switch(mcmc)
#'
#' # plot results for highest posterior probability model
#' plot(mcmc, what = 'hppm')
#'
#' @export
#' @importFrom ggplot2 ggplot aes_string theme labs element_text
#' @importFrom ggplot2 geom_line geom_tile geom_text geom_rug geom_bar
#' @importFrom ggplot2 scale_fill_gradient2 guide_colorbar guides unit

plot.befa <- function(x, ...)
{

  args <- list(...)
  show.val <- ifelse(is.null(args$show.val), TRUE, args$show.val)
  what <- ifelse(is.null(args$what), 'maxp', args$what)
  assertFlag(show.val)
  if (!what %in% c('maxp', 'hppm'))
    stop('plot.befa() only implemented for what = "maxp" and what = "hppm"')

  Kmax <- attr(x, 'Kmax')

  ##############################################################################
  ### trace of number of factors

  dat <- data.frame(nfac  = factor(x$nfac, levels = 0:Kmax),
                    iter  = as.numeric(names(x$nfac)),
                    MHacc = as.numeric(x$MHacc))

  p.nfac <- ggplot(dat, aes_string(x = 'iter', y = 'nfac')) +
              geom_line(colour = 'steelblue')

  p.nfac <- p.nfac + labs(title = paste0('trace plot of number of factors\n',
                             '(accepted Metropolis-Hastings draws at bottom)'),
                         x = 'MCMC iterations',
                         y = 'number of factors')

  # add Metropolis-Hastings acceptance
  p.nfac <- p.nfac + geom_rug(aes_string(y = 'MHacc'), sides = 'b',
                              colour = 'darkcyan')


  # posterior probabilities of number of factors
  nft <- table(factor(x$nfac, levels = 0:Kmax))
  dat <- data.frame(nfac = as.factor(0:Kmax),
                    freq = as.numeric(nft / length(x$nfac)))
  p.hnfac <- ggplot(dat, aes_string(x = 'nfac')) +
               geom_bar(aes_string(weight = 'freq'), fill = 'steelblue') +
               labs(title = 'posterior probabilities of number of factors',
                    x = 'number of factors',
                    y = 'frequencies')


  ##############################################################################
  # summarize and plot

  x <- summary(x, ...)
  if (what == 'hppm') {
    alpha <- x$alpha$m1
    dedic <- x$alpha$m1$dedic
    R     <- x$R$m1
  } else {
    alpha <- x$alpha
    dedic <- x$alpha$dedic
    R     <- x$R
  }
  nvar  <- length(dedic)

  ### matrix of indicator probabilities

  if (what != 'hppm') {     # skip for HPP model

    pind <- matrix(NA, nvar, Kmax)
    rownames(pind) <- sapply(strsplit(rownames(alpha), ':'), '[', 2)
    colnames(pind) <- paste0('f', 1:Kmax)
    for (i in 1:nvar)
      pind[i, dedic[i]] <- alpha$prob[i]

    # which factors are loaded by at least one measurement?
    acti <- apply(!is.na(pind), 2, any)

    # heatmap for active factors only
    p.indic <- make.heatmap(pind[, acti],
                            title = 'indicator probabilities of being nonzero',
                            xlab = 'latent factors (active factors only)',
                            ylab = 'manifest variables',
                            show.val)

  }

  ### matrix of factor loadings

  # construct matrix from factor loadings and indicators
  # (remove the 'alpha:' part from variable names to simplify plot)
  alpha.post <- matrix(NA, nvar, Kmax)
  rownames(alpha.post) <- sapply(strsplit(rownames(alpha), ':'), '[', 2)
  colnames(alpha.post) <- paste0('f', 1:Kmax)
  for (i in 1:nvar)
    alpha.post[i, dedic[i]] <- alpha$mean[i]

  # which factors are loaded by at least one measurement?
  acti <- !apply(is.na(alpha.post), 2, all)

  # heatmap for active factors only
  p.alpha <- make.heatmap(alpha.post[, acti],
                          title = 'factor loading matrix',
                          xlab = 'latent factors (active factors only)',
                          ylab = 'manifest variables',
                          show.val)

  ### correlation matrix of the factors

  # construct matrix from lower diagonal elements
  Rmat <- .5 * diag(Kmax)
  Rmat[lower.tri(Rmat)] <- R$mean
  Rmat <- Rmat + t(Rmat)
  rownames(Rmat) <- colnames(Rmat) <- paste0('f', 1:Kmax)

  # heatmap for active factors only
  p.Rmat <- make.heatmap(Rmat[acti, acti],
                         title = 'correlation matrix of the factors',
                         xlab = 'latent factors (active factors only)',
                         ylab = '',
                         show.val)

  ##############################################################################

  print(p.nfac)
  invisible(readline(prompt = "Press <Enter> to show next graph..."))
  print(p.hnfac)
  invisible(readline(prompt = "Press <Enter> to show next graph..."))
  if (what != 'hppm') {
    print(p.indic)
    invisible(readline(prompt = "Press <Enter> to show next graph..."))
  }
  print(p.alpha)
  invisible(readline(prompt = "Press <Enter> to show next graph..."))
  print(p.Rmat)

}


################################################################################

make.heatmap <- function(x, title, xlab, ylab, show.val) {

  # prepare data
  xcol <- colnames(x)
  xrow <- rownames(x)
  dat <- data.frame(xvar = factor(rep(xcol, each = nrow(x)), levels = xcol),
                    yvar = factor(rep(xrow, ncol(x)), levels = rev(xrow)),
                    val  = c(round(x, digits = 2)))

  # make heatmap
  p <- ggplot(dat, aes_string(x = 'xvar', y = 'yvar')) +
    geom_tile(aes_string(fill = 'val')) +
    scale_fill_gradient2(low  = 'steelblue',
                         mid  = 'white',
                         high = 'red',
                         na.value = NA)

  # add title and labels
  p <- p + labs(title = title, x = xlab, y = ylab, fill = 'posterior\nmean')

  # adjust font sizes
  # p <- p + theme(axis.title  = element_text(size=12),
  #                axis.text   = element_text(size=12),
  #                plot.title  = element_text(size=16),
  #                legend.text = element_text(size=10))

  # legend
  p <- p + guides(fill = guide_colorbar(barheight = unit(6, 'cm')))

  # add values to heatmap?
  if (show.val)
    p <- p + geom_text(aes_string(label = 'val'),
                       colour = 'white',
                       na.rm = TRUE)

  return(p)

}
