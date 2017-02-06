#' @export

print.simul.R.prior <- function(x, ...)
{
  cat('sample from prior distribution of factor correlation matrix\n')
  cat('\ncall:\n  ')
  print(attr(x, 'call'))
  cat('\nKmax =', attr(x, 'Kmax'), '\n')
  cat('\nparameter values:\n  ')
  cat(names(x), sep = '\n  ')
  cat('\nnrep =', attr(x, 'nrep'), '\n\n')
  cat('use summary() and plot() functions to summarize and plot prior density',
      'of maximum correlation and/or of minimum eigenvalue', sep = '\n')
}
