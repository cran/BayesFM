#' @export

print.summary.befa <- function(x, ...)
{

  cat('BEFA - Bayesian Exploratory Factor Analysis',
      'Summary of posterior results\n', sep = '\n')

  args <- list(...)
  digits <- args$digits
  if (is.null(digits)) digits <- 3

  cat('Maximum number of factors (Kmax) =', attr(x, 'Kmax'), '\n')
  cat('Identification restriction (Nid) =', attr(x, 'Nid'), '\n\n')

  cat('MCMC iterations =', attr(x, 'iter'), '\n')
  cat('burn-in period  =', attr(x, 'burnin'), '\n\n')

  cat('Metropolis-Hastings acceptance rate =', round(x$MHacc, digits = 3))
  cat('\n\n')

  what     <- attr(x, 'what')
  hpd.prob <- attr(x, 'hpd.prob')

  if (what == 'hppm') {

    cat('-----------------------------------------------------------------',
        'Indicator matrix\n',
        'Dedicated structures visited by the sampler, with corresponding ',
        'posterior probabilities and number of factors (K)\n', sep = '\n')
    cat(paste0('prob(', colnames(x$hppm$dedic), ') = ',
               sprintf('%.3f', x$hppm$prob),
               ',  K = ', x$hppm$nfac),
        sep = '\n')
    cat('\n')
    print(x$hppm$dedic, print.gap = 3)
    cat('\n')
  } else {

    cat('Posterior frequency of number of latent factors:\n')
    freq <- sprintf('%3.2f%%', 100*x$nfac)
    nf   <- sprintf('%2i', as.integer(names(x$nfac)))
    cat(paste('  K =', nf, '   ', freq), sep='\n')
    cat('\n')

  }

  # print MCMC summary results for model parameters
  cat('-----------------------------------------------------------------',
      'Model parameters\n',
      sep = '\n')

  print.mcmc <- function(z, title) {
    res <- round(z, digits = digits)
    collab <- c('mean', 'sd',
                sprintf(paste0('[%', digits, '.0f%%'), 100*hpd.prob),
                sprintf(paste0('%', -(digits+1), 's]'), 'hpd'))
    if (ncol(z) == 5) collab <- c('dedic', collab)
    if (ncol(z) == 6) collab <- c('dedic', 'prob', collab)
    colnames(res) <- collab
    if (!missing(title)) cat(title, ':\n\n', sep='')
    print(as.matrix(res), print.gap = 3, na.print = ' ')
    cat('\n')
  }

  print.mcmc.all <- function(alpha, sigma, R, beta) {
    print.mcmc(alpha, title = 'Factor loadings')
    print.mcmc(sigma, title = 'Idiosyncratic variances')
    print.mcmc(R,     title = 'Factor correlation matrix')
    if (!is.null(beta)) print.mcmc(beta, title = 'Regression coefficients')
  }

  if (what == 'hppm') {

    for (i in seq_along(x$hppm$prob)) {
      cat('-----------------------------------------------------------------\n')
      cat('Highest posterior probability model ', i, '\n\n',
          'prob = ', round(x$hppm$prob[i], digits = 3), '\n',
          '   K = ', x$hppm$nfac[i], '\n\n',
          sep = '')
      print.mcmc.all(x$alpha[[i]], x$sigma[[i]], x$R[[i]], x$beta[[i]])
    }

  } else {

    print.mcmc.all(x$alpha, x$sigma, x$R, x$beta)

  }

}
