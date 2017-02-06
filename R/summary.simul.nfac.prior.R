#' @export

summary.simul.nfac.prior <- function(object, ...)
{
  Kmax <- attr(object, 'Kmax')

  comp.freq <- function(x) table(factor(x, levels = 1:Kmax))/length(x)

  acc  <- lapply(object, '[[', 'acc')
  nfac <- lapply(object, '[[', 'nfac')
  nfac <- lapply(nfac, comp.freq)

  out <- list(nfac = nfac, acc = acc)
  class(out) <- 'summary.simul.nfac.prior'
  return(out)
}
