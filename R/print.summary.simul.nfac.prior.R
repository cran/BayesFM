#' @export

print.summary.simul.nfac.prior <- function(x, ...)
{
  cat('prior probabilities of numbers of factors:\n')
  out <- cbind(do.call(rbind, x$nfac), acc = unlist(x$acc))
  out <- round(out, digits = 3)
  print(out, print.gap = 3)
}
