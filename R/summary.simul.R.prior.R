#' @export
#' @importFrom stats median sd quantile

summary.simul.R.prior <- function(object, ...)
{

  sum.prior <- function(x, FUN) {
    val <- lapply(x, function(z) apply(z, 3, FUN))
    stat <- lapply(val, function(z) return(c(mean   = mean(z),
                                            median = median(z),
                                            sd     = sd(z),
                                            min    = min(z),
                                            max    = max(z))))
    quant <- lapply(val, quantile, probs = c(.05, .1, .25, .75, .9, .95))
    res <- list(stat  = do.call(rbind, stat),
                quant = do.call(rbind, quant))
    res <- lapply(res, round, digits = 3)
    return(res)
  }

  out <- list(maxcor = sum.prior(object, get.maxcor),
              mineig = sum.prior(object, get.mineig))
  class(out) <- 'summary.simul.R.prior'
  return(out)

}
