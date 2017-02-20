#' @export

print.simul.nfac.prior <- function(x, ...)
{
  cat('sample from prior distribution of number of latent factors\n')
  cat('\ncall:\n  ')
  print(attr(x, 'call'))
  cat('\n')
  cat('nvar =',  attr(x, 'nvar'), '\n')
  cat('Kmax  =', attr(x, 'Kmax'), '\n')
  cat('Nid   =', attr(x, 'Nid'), '\n')
  cat('kappa =', attr(x, 'kappa'), '\n')
  cat('nrep  =', attr(x, 'nrep'), '\n\n')
  cat('use summary() and plot() functions to summarize and plot',
      'prior distribution of number of latent factors', sep = '\n')
}
