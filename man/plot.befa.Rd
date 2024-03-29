% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.befa.R
\name{plot.befa}
\alias{plot.befa}
\title{Plot object of class 'befa'}
\usage{
\method{plot}{befa}(x, ...)
}
\arguments{
\item{x}{Object of class 'befa'.}

\item{...}{The following extra arguments can be specified:
\itemize{
  \item \code{what}: How to summarize the posterior distribution?
        \itemize{
          \item \code{what = 'maxp'} (default): Only factor loadings with
                highest posterior probability of being different from zero or
                discarded from the model (if \code{dedic = 0}) are
                summarized.
          \item \code{what = 'all'}: All factor loadings with corresponding
                posterior probability to be allocated to a given factor (or
                to be discarded from the model) larger than \code{min.prob}
                are summarized.
          \item \code{what = 'hppm'}: Highest posterior probability models
                with probability larger than \code{min.prob} are summarized.
        }
  \item \code{byfac}: Sort factor loadings by factors if \code{TRUE},
        otherwise by manifest variables if \code{FALSE} (default).
  \item \code{hpd.prob}: Probability used to compute the highest posterior
        density intervals of the posterior distribution of the model
        parameters (default: 0.95).
  \item \code{min.prob}: If \code{what = 'all'}, only factor loadings with
        posterior probability of being dedicated to a given factor (or
        discarded from the model) larger than this value are displayed.
        If \code{what = 'hppm'}, only highest posterior probability models
        with probability larger than this value are displayed. (default:
        0.20)
}}
}
\value{
No return value, called for side effects (plots the posterior results 
returned by \code{\link{summary.befa}}).
}
\description{
This function makes different plots that are useful to assess the posterior
results: a trace plot of the number of latent factors (also showing
Metropolis-Hastings acceptance across MCMC replications), a histogram
of the posterior probabilities of the number of factors, heatmaps for the
inficator probabilities, the factor loading matrix, and the correlation
matrix of the latent factors.
}
\details{
This function makes graphs based on the summary results returned by
\code{\link{summary.befa}}. It therefore accepts the same optional arguments
as this function.
}
\examples{
set.seed(6)

# generate fake data with 15 manifest variables and 3 factors
Y <- simul.dedic.facmod(N = 100, dedic = rep(1:3, each = 5))

# run MCMC sampler and post process output
# notice: 1000 MCMC iterations for illustration purposes only,
#  increase this number to obtain reliable posterior results!
mcmc <- befa(Y, Kmax = 5, iter = 1000)
mcmc <- post.column.switch(mcmc)
mcmc <- post.sign.switch(mcmc)

# plot results for highest posterior probability model
plot(mcmc, what = 'hppm')

}
\seealso{
\code{\link{summary.befa}} to summarize posterior results.
}
\author{
Rémi Piatek \email{remi.piatek@gmail.com}
}
