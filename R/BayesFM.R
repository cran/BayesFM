#' BayesFM: Package for Bayesian Factor Modeling
#'
#' The long-term goal of this package is to provide a collection of procedures
#' to perform Bayesian inference on a variety of factor models.
#'
#' @details Currently, this package includes: Bayesian Exploratory Factor
#' Analysis (\code{befa}), as developed in Conti et al. (2014), an approach to
#' dedicated factor analysis with stochastic search on the structure of the
#' factor loading matrix. The number of latent factors, as well as the
#' allocation of the observed variables to the factors, are not fixed a priori
#' but determined during MCMC sampling. More approaches will be included in
#' future releases of this package.
#'
#' @note You are very welcome to send me any comments or suggestions for 
#' improvements, and to share with me any problems you may encounter with the 
#' use of this package.
#'
#' @author Rémi Piatek \email{remi.piatek@@gmail.com}
#'
#' @references G. Conti, S. Frühwirth-Schnatter, J.J. Heckman, R. Piatek (2014):
#' ``Bayesian Exploratory Factor Analysis'', \emph{Journal of Econometrics},
#' 183(1), pages 31-57, \doi{10.1016/j.jeconom.2014.06.008}.
#'
#' @docType package
#' @name BayesFM
NULL

.onAttach <- function(libname, pkgname) {
  if (interactive() || getOption("verbose")) {
    msg <- sprintf(paste(
      "###",
      "### Package %s (%s) loaded",
      "###",
      "### Please report any bugs, and send suggestions or feedback",
      "### to %s",
      "###", sep = "\n"),
      pkgname,
      utils::packageDescription(pkgname)$Version,
      utils::maintainer(pkgname))
    packageStartupMessage(msg)
  }
}

.onUnload <- function(libpath) {
  library.dynam.unload("BayesFM", libpath)
}
