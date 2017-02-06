#' @export

print.summary.simul.R.prior <- function(x, ...)
{

  cat('prior distribution of correlation matrix of latent factors\n\n')

  cat('maximum correlation (in absolute value):\n\n')
  print(x$maxcor$stat,  print.gap = 3)
  cat('\n')
  print(x$maxcor$quant, print.gap = 3)
  cat('\n')

  cat('minimum eigenvalue of correlation matrix:\n\n')
  print(x$mineig$stat,  print.gap = 3)
  cat('\n')
  print(x$mineig$quant, print.gap = 3)

}
