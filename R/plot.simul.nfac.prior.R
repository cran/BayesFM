#' @export
#' @importFrom ggplot2 ggplot aes_string geom_bar labs theme element_blank

plot.simul.nfac.prior <- function(x, ...)
{
  nfac <- summary(x)$nfac
  dat  <- data.frame(freq = unlist(nfac),
                     nfac = c(as.integer(sapply(nfac, names))),
                     lab  = rep(names(nfac), sapply(nfac, length)))

  ggplot(dat, aes_string(x = 'nfac', y = 'freq')) +
    geom_bar(aes_string(fill = 'lab'), position = 'dodge', stat = 'identity') +
    labs(x = 'number of factors', y = 'prior probability') +
    theme(legend.position  = 'bottom',
          legend.direction = 'horizontal',
          legend.title     = element_blank())
}
