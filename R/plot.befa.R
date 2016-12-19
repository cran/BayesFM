#' @export
#' @importFrom graphics axis barplot layout par plot points

plot.befa <- function(x, ...)
{

  cat('Trace plot and frequency of number of active factors are displayed in ',
      'blue. \nMetropolis-Hastings proposals accepted across MCMC iterations ',
      'are shown in green (prob = ', round(mean(x$MHacc), digits = 3), ').\n\n',
      sep = '')

  nfac  <- x$nfac
  MHacc <- which(x$MHacc)

  k0 <- min(nfac)
  k1 <- max(nfac)

  # pos.MHacc: location of bars on vertical axis for accepted M-H proposals
  if (k0 == k1) {
    pos.MHacc <- k0 - .5
    k1 <- k0 + .5
  } else {
    pos.MHacc <- k0 - .1*(k1 - k0)
  }

  layout(matrix(c(1, 2), nrow = 1))
  par(cex = 1)

  plot(nfac,
       main = 'Trace',
       xlab = 'MCMC iterations',
       ylab = 'factors',
       type = 'l',
       lwd  = 2,
       col  = 'royalblue',
       ylim = c(pos.MHacc, k1),
       yaxt = 'n')
  points(MHacc, y = rep(pos.MHacc, length(MHacc)),
         pch = '|',
         cex = .5,
         col = 'darkcyan')
  axis (2, at = k0:k1)

  freq <- table(nfac)/length(nfac)
  barplot(freq,
          main   = 'Frequency',
          xlab   = 'factors',
          ylab   = 'freq',
          ylim   = c(0, .1*ceiling(10*max(freq))),
          col    = 'royalblue',
          border = NA)

}
