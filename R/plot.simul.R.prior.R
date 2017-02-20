#' @export
#' @importFrom ggplot2 ggplot aes_string geom_density theme labs element_blank
#' @importFrom ggplot2 ggplot_gtable ggplot_build
#' @importFrom gridExtra grid.arrange arrangeGrob

plot.simul.R.prior <- function(x, ...)
{

  if (!is.list(x)) x <- as.list(x)
  nsim <- length(x)

  # function making individual plot
  plot.dens <- function(Rsim, FUN, xlab) {

    val <- lapply(Rsim, function(z) apply(z, 3, FUN))
    val <- do.call(cbind, val)
    dat <- data.frame(val = c(val),
                      lab = rep(colnames(val), each = nrow(val)))

    # make plot
    p <- ggplot(dat, aes_string(x = 'val', color = 'lab')) +
           geom_density(kernel = 'gaussian') +
           labs(x = xlab, y = 'prior density') +
           theme(legend.position  = 'bottom',
                 legend.direction = 'horizontal',
                 legend.title     = element_blank())
    return(p)

  }

  # make plots
  plots <- list()
  plots$maxcor <- plot.dens(x, FUN = get.maxcor, xlab = 'max(|R|)')
  plots$mineig <- plot.dens(x, FUN = get.mineig, xlab = 'min(eigen(R))')

  # get common legend
  legend <- ggplot_gtable(ggplot_build(plots[[1]]))
  legid  <- which(sapply(legend$grobs, function(p) p$name) == 'guide-box')
  legend <- legend$grobs[[legid]]

  # make plot
  grid.arrange(arrangeGrob(plots[[1]] + theme(legend.position = 'none'),
                           plots[[2]] + theme(legend.position = 'none'),
                           nrow = 1),
               legend, nrow = 2, heights = c(10, 1),
               top = 'correlation matrix of the latent factors')

}
