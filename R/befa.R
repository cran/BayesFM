#'
#' Bayesian Exploratory Factor Analysis
#'
#' This function implements the Bayesian Exploratory Factor Analysis
#' (\code{befa}) approach developed in Conti et al. (CFSHP, 2014). It runs a
#' MCMC sampler for a factor model with dedicated factors, where each manifest
#' variable is allowed to load on at most one latent factor. The allocation of
#' the manifest variables to the latent factors is not fixed \emph{a priori} but
#' determined stochastically during sampling. The minimum number of variables
#' dedicated to each factor can be controlled by the user to achieve the desired
#' level of identification. The manifest variables can be continuous or
#' dichotomous, and control variables can be introduced as covariates.
#'
#' @param model
#'        This argument specifies the manifest variables and the covariates used
#'        in the model (if any). It can be specified in two different ways:
#'        \itemize{
#'        \item A numeric matrix or a data frame containing the manifest
#'          variables. This corresponds to a model without covariates,
#'          and the argument \code{data} is not required.
#'        \item A list of model formulas. Each element of the list specifies
#'          a manifest variable and its corresponding control variables (e.g.,
#'          '\code{Y1 ~ X1 + X2}' to use \code{X1} and \code{X2} as control
#'          variables for \code{Y1}).
#'          If a formula has no left-hand side variable, covariates on the
#'          right-hand side are included in all equations (e.g., '\code{~ X3}'
#'          means that \code{X3} is used as a control variable for all the
#'          manifest variables). Argument \code{data} can be passed to the
#'          function in that case, otherwise parent data frame is used.
#'        }
#'        Binary manifest variables should be specified as logical vectors in
#'        the data frame to be treated as dichotomous. \code{NA} values are
#'        accepted in manifest variables only.
#' @param data
#'        Data frame. If missing, parent data frame if used.
#' @param burnin
#'        Burn-in period of the MCMC sampler.
#' @param iter
#'        Number of MCMC iterations saved for posterior inference (after
#'        burn-in).
#' @param Nid
#'        Minimum number of manifest variables dedicated to each latent factor
#'        for identification.
#' @param Kmax
#'        Maximum number of latent factors. If missing, the maximum number of
#'        factors that satisfies the identification condition determined by
#'        \code{Nid} and the Ledermann bound is specified (see CFSHP,
#'        section 2.2).
#' @param A0
#'        Scaling parameters of the variance of the Normal prior on the nonzero
#'        factor loadings. Either a scalar or a numeric vector of length equal
#'        to the number of manifest variables.
#' @param B0
#'        Variances of the Normal prior on the regression coefficients. Either a
#'        scalar or a numeric vector of length equal to the number of
#'        manifest variables.
#' @param c0
#'        Shape parameters of the Inverse-Gamma prior on the idiosyncratic
#'        variances. Either a scalar or a numeric vector of length equal to the
#'        number of manifest variables.
#' @param C0
#'        Scale parameters of the Inverse-Gamma prior on the idiosyncratic
#'        variances. Either a scalar or a numeric vector of length equal to the
#'        number of manifest variables.
#' @param HW.prior
#'        If \code{TRUE}, implement Huang-Wand (2013) prior on the covariance
#'        matrix of the factors in the expanded model, otherwise use an
#'        Inverse-Wishart prior if \code{FALSE}, see CFSHP section 2.3.5.
#' @param nu0
#'        Degrees of freedom of the Inverse-Wishart prior on the covariance
#'        matrix of the latent factors in the expanded model.
#' @param S0
#'        Scale parameters of the Inverse-Wishart prior on the covariance matrix
#'        of latent factors in the expanded model:
#'        \itemize{
#'          \item if \code{HW.prior = TRUE}, scale parameter of the Gamma
#'          hyperprior distribution on the individual scales of the
#'          Inverse-Wishart prior.
#'          \item if \code{HW.prior = FALSE}, diagonal elements of the scale
#'          matrix of the Inverse-Wishart prior on the covariance matrix of
#'          the latent factors in the expanded model.
#'        }
#'        Either a scalar or a numeric vector of length equal to \code{Kmax}.
#' @param kappa0
#'        First shape parameter of the Beta prior distribution on the
#'        probability \eqn{\tau_0} that a manifest variable does not load on any
#'        factor.
#' @param xi0
#'        Second shape parameter of the Beta prior distribution on the
#'        probability \eqn{\tau_0} that a manifest variable does not load on any
#'        factor.
#' @param kappa
#'        Concentration parameters of the Dirichlet prior distribution on the
#'        indicators. Either a scalar or a numeric vector of length equal to
#'        \code{Kmax}.
#' @param indp.tau0
#'        If \code{TRUE}, specify the alternative prior specification with
#'        independent parameters \eqn{\tau_{0m}}{\tau_0m} across manifest
#'        variables \eqn{m = 1, ..., M}, otherwise use a common parameter
#'        \eqn{\tau_0} if \code{FALSE}.
#' @param rnd.step
#'        If \code{TRUE}, select randomly the number of intermediate steps in
#'        non-identified models at each MCMC iteration, otherwise use a fixed
#'        number of steps if \code{FALSE}.
#' @param n.step
#'        Controls the number of intermediate steps in non-identified models:
#'        \itemize{
#'          \item if \code{rnd.step = TRUE}, average number of steps. The number
#'          of steps is sampled at each MCMC iteration from
#'          1+Poisson(\code{n.step}-1).
#'          \item if \code{rnd.step = FALSE}, fixed number of steps.
#'        }
#' @param search.delay
#'        Number of MCMC iterations run with fixed indicator matrix (specified
#'        in \code{dedic.start}) at beginning of MCMC sampling.
#' @param R.delay
#'        Number of MCMC iterations run with fixed correlation matrix (specified
#'        in \code{dedic.start}) at beginning of MCMC sampling.
#' @param alpha.start
#'        Starting values for the factor loadings. Numeric vector of length
#'        equal to the number of manifest variables. If missing, sampled from a
#'        Normal distribution with zero mean and variance \code{A0}.
#' @param dedic.start
#'        Starting values for the indicators. Vector of integers of length equal
#'        to the number of manifest variables. Each element takes a value among
#'        0, 1, ..., \code{Kmax}. If missing, random allocation of the manifest
#'        variables to the maximum number of factors \code{Kmax}, with a minimum
#'        of \code{Nid} manifest variables loading on each factor.
#' @param sigma.start
#'        Starting values for the idiosyncratic variances. Numeric vector of
#'        length equal to the number of manifest variables. Sampled from prior
#'        if missing.
#' @param beta.start
#'        Starting values for the regression coefficients. Numeric vector of
#'        length equal to the total number of regression coefficients,
#'        concatenated for all the manifest variables. Sampled from prior if
#'        missing.
#' @param R.start
#'        Starting values for the correlation matrix of the latent factors.
#'        Numeric matrix with \code{Kmax} rows and columns, and unit diagonal
#'        elements. If missing, identity matrix is used.
#' @param verbose
#'        If \code{TRUE}, display information on the progress of the function.
#'
#' @details \strong{Model specification.} The model is specified as follows, for
#' each observation \eqn{i = 1, ..., N}:
#'   \deqn{Y^*_i = \beta X_i + \alpha \theta_i + \epsilon_i}{
#'         Y*_i = \beta X_i + \alpha \theta_i + \epsilon_i}
#'   \deqn{\theta_i \sim \mathcal{N}(0, R)}{\theta_i ~ N(0, R)}
#'   \deqn{\epsilon_i \sim \mathcal{N}(0, \Sigma)}{\epsilon_i ~ N(0, \Sigma)}
#'   \deqn{\Sigma = diag(\sigma^2_1, ..., \sigma^2_M)}
#' where \eqn{Y^*_i}{Y*_i} is the \eqn{M}-vector containing the latent
#' variables underlying the corresponding \eqn{M} manifest variables
#' \eqn{Y_i}, which can be continuous such that
#' \eqn{Y_{im} = Y^*_{im}}{Y_im = Y*_im}, or binary with
#' \eqn{Y_{im} = 1[Y^*_{im} > 0]}{Y_im = 1[Y*_im > 0]}, for \eqn{m = 1, ..., M}.
#' The \eqn{K}-vector \eqn{\theta_i} contains the latent factors, and
#' \eqn{\alpha} is the \eqn{(M \times K)}{(M*K)}-matrix of factor loadings.
#' The \eqn{M}-vector \eqn{\epsilon_i} is the vector of error terms.
#' Covariates can be included in the \eqn{Q}-vector \eqn{X_i} and are related to
#' the manifest variables through the \eqn{(M \times Q)}{(M*Q)}-matrix of
#' regression coefficients \eqn{\beta}. Intercept terms are automatically
#' included, but can be omitted in some or all equations using the usual syntax
#' for R formulae (e.g., 'Y1 ~ X1 - 1' specifies that that Y1 is regressed on X1
#' and no intercept is included in the corresponding equation).
#'
#' The number of latent factors \eqn{K} is specified as \code{Kmax}. However,
#' during MCMC sampling the stochastic search process on the matrix \eqn{\alpha}
#' may produce zero columns, thereby reducing the number of active factors.
#'
#' The covariance matrix \eqn{R} of the latent factors is assumed to be a
#' correlation matrix for identification.
#'
#' Each row of the factor loading matrix \eqn{\alpha} contains at most one
#' nonzero element (dedicated factor model). The allocation of the manifest
#' variables to the latent factors is indicated by the binary matrix
#' \eqn{\Delta} with same dimensions as \eqn{\alpha}, such that each row
#' \eqn{\Delta_m} indicates which factor loading is different from zero, e.g.:
#' \deqn{\Delta_m = (0, .., 0, 1, 0, ..., 0) \equiv e_k}{
#'       \Delta_m = (0, .., 0, 1, 0, ..., 0) = e_k}
#' indicates that variable \eqn{m} loads on the \eqn{k}th factor, where
#' \eqn{e_k} is a \eqn{K}-vector that contains only zeros, besides its \eqn{k}th
#' element that equals 1.
#'
#' \strong{Identification.} The function verifies that the maximum number of
#' latent factors \code{Kmax} does not exceed the Ledermann bound. It also
#' checks that \code{Kmax} is consistent with the identification restriction
#' specified with \code{Nid} (enough variables should be available to load on
#' the factors when \code{Kmax} is reached). The default value for \code{Kmax}
#' is the minimum between the Ledermann bound and the maximum number of factors
#' that can be loaded by \code{Nid} variables. The user is free to select the
#' level of identification, see CFSHP section 2.2 (non-identified models are
#' allowed with \code{Nid = 1}).
#'
#' Note that identification is achieved only with respect to the scale of the
#' latent factors. Non-identifiability problems may affect the posterior sample
#' because of column switching and sign switching of the factor loadings.
#' These issues can be addressed \emph{a posteriori} with the functions
#' \code{\link{post.column.switch}} and \code{\link{post.sign.switch}}.
#'
#' \strong{Prior specification.}
#' The indicators are assumed to have the following probabilities,
#' for \eqn{k = 1, ..., K}:
#' \deqn{Prob(\Delta_m = e_k \mid \tau_k) = \tau_k}{
#'       Prob(\Delta_m = e_k | \tau_k) = \tau_k}
#' \deqn{\tau = (\tau_0, \tau_1, ..., \tau_K)}
#' If \code{indp.tau0 = FALSE}, the probabilities are specified as:
#' \deqn{\tau = [\tau_0, (1-\tau_0)\tau^*_1, ..., (1-\tau_0)\tau^*_K]}{
#'       \tau = [\tau_0, (1-\tau_0)\tau*_1, ..., (1-\tau_0)\tau*_K]}
#' \deqn{\tau_0 \sim \mathcal{B}eta(\kappa_0, \xi_0)}{
#'       \tau_0 ~ Beta(\kappa_0, \xi_0)}
#' \deqn{\tau^* = (\tau^*_1, ..., \tau^*_K) \sim \mathcal{D}ir(\kappa)}{
#'       \tau* = (\tau*_1, ..., \tau*_K) ~ Dir(\kappa)}
#' with \eqn{\kappa_0} = \code{kappa0}, \eqn{\xi_0} = \code{xi0} and
#' \eqn{\kappa} = \code{kappa}.
#' Alternatively, if \code{indp.tau0 = TRUE}, the probabilities are specified
#' as:
#' \deqn{\tau_m = [\tau_{0m}, (1-\tau_{0m})\tau^*_1, ...,
#'         (1-\tau_{0m})\tau^*_K]}{
#'       \tau_m = [\tau_0m, (1-\tau_0m)\tau*_1, ..., (1-\tau_0m)\tau*_K]}
#' \deqn{\tau_{0m} \sim \mathcal{B}eta(\kappa_0, \xi_0)}{
#'       \tau_0m ~ Beta(\kappa_0, \xi_0)}
#' for each manifest variable \eqn{m = 1, ..., M}.
#'
#' A normal-inverse-Gamma prior distribution is assumed on the nonzero factor
#' loadings and on the idiosyncratic variances:
#' \deqn{\sigma^2_m \sim \mathcal{I}nv-\mathcal{G}amma(c_{0m}, C_{0m})}{
#'       \sigma^2_m ~ Inv-Gamma(c0_m, C0_m)}
#' \deqn{\alpha_m^\Delta \mid \sigma^2_m \sim \mathcal{N}(0, A_{0m}\sigma^2_m)}{
#'       \alpha_m^\Delta | \sigma^2_m ~ N(0, A0_m * \sigma^2_m)}
#' where \eqn{\alpha_m^\Delta} denotes the only nonzero loading in the \eqn{m}th
#' row of \eqn{\alpha}.
#'
#' For the regression coefficients, a multivariate Normal prior distribution is
#' assumed on each row \eqn{m = 1, ..., M} of \eqn{\beta}:
#' \deqn{\beta_m \sim \mathcal{N}(0, B_0 I_Q)}{\beta_m ~ N(0, B_0 I_Q)}
#' The covariates can be different across manifest variables, implying zero
#' restrictions on the matrix \eqn{\beta}. To specify covariates, use a list
#' of formulas as \code{model} (see example below). Intercept terms can be
#' introduced using
#'
#' To sample the correlation matrix \eqn{R} of the latent factors, marginal data
#' augmentation is implemented (van Dyk and Meng, 2001), see CFSHP section 2.2.
#' Using the transformation \eqn{\Omega = \Lambda^{1/2} R \Lambda^{1/2}}, the
#' parameters \eqn{\Lambda = diag(\lambda_1, ..., \lambda_K)} are used as
#' \emph{working parameters}. These parameters correspond to the variances of
#' the latent factors in an expanded version of the model where the factors do
#' not have unit variances. Two prior distributions can be specified on the
#' covariance matrix \eqn{\Omega} in the expanded model:
#' \itemize{
#'   \item If \code{HW.prior = FALSE}, inverse-Wishart distribution:
#'   \deqn{\Omega \sim \mathcal{I}nv-\mathcal{W}ishart(\nu_0, diag(S_0))}{
#'         \Omega ~ Inv-Wishart(\nu_0, diag(S0))}
#'   with \eqn{\nu_0} = \code{nu0} and \eqn{S_0} = \code{S0}.
#'   \item If \code{HW.prior = TRUE}, Huang-Wand (2013) prior:
#'   \deqn{\Omega \sim \mathcal{I}nv-\mathcal{W}ishart(\nu_0, W), \qquad
#'         W = diag(w_1, ..., w_K)}{
#'         \Omega ~ Inv-Wishart(nu0, W), W = diag(w_1, ..., w_K)}
#'   \deqn{w_k \sim \mathcal{G}amma\left(\frac{1}{2},
#'         \frac{1}{2\nu^*S_{0k}}\right)}{w_k ~ Gamma(1/2, 1/(2\nu*S0_k))}
#'   with \eqn{\nu^*}{\nu*} = \code{nu0} - \code{Kmax} + 1, and the shape and
#'   rate parameters are specified such that the mean of the gamma distribution
#'   is equal to \eqn{\nu^* S_{0k}}{\nu* S0_k}, for each \eqn{k = 1, ..., K}.
#' }
#'
#' \strong{Missing values.} Missing values (\code{NA}) are allowed in the
#' manifest variables. They are drawn from their corresponding conditional
#' distributions during MCMC sampling. Control variables with missing values
#' can be passed to the function. However, all the observations with at least
#' one missing value in the covariates are discarded from the sample (a warning
#' message is issued in that case).
#'
#' @return The function returns an object of class '\code{befa}' containing the
#' MCMC draws of the model parameters saved in the following matrices (each
#' matrix has '\code{iter}' rows):
#' \itemize{
#'   \item \code{alpha}: Factor loadings.
#'   \item \code{sigma}: Idiosyncratic variances.
#'   \item \code{R}: Correlation matrix of the latent factors (off-diagonal
#'         elements only).
#'   \item \code{beta}: regression coefficients (if any).
#'   \item \code{dedic}: indicators (integers indicating on which factors the
#'         manifest variable load).
#'   }
#' The returned object also contains:
#' \itemize{
#'   \item \code{nfac}: Vector of number of 'active' factors across MCMC
#'         iterations (i.e., factors loaded by at least \code{Nid} manifest
#'         variables).
#'   \item \code{MHacc}: Logical vector indicating accepted proposals of
#'         Metropolis-Hastings algorithm.
#' }
#' The parameters \code{Kmax} and \code{Nid} are saved as object attributes, as
#' well as the function call and the number of mcmc iterations (\code{burnin}
#' and \code{iter}), and two logical variables indicating if the returned object
#' has been post processed to address the column switching problem
#' (\code{post.column.switch}) and the sign switching problem
#' (\code{post.sign.switch}).
#'
#' @author Rémi Piatek \email{remi.piatek@@econ.ku.dk}
#'
#' @references G. Conti, S. Frühwirth-Schnatter, J.J. Heckman, R. Piatek (2014):
#' ``Bayesian Exploratory Factor Analysis'', \emph{Journal of Econometrics},
#' 183(1), pages 31-57, \url{http://dx.doi.org/10.1016/j.jeconom.2014.06.008}.
#'
#' @references A. Huang, M.P. Wand (2013):
#' ``Simple Marginally Noninformative Prior Distributions for Covariance
#' Matrices'', \emph{Bayesian Analysis}, 8(2), pages 439-452,
#' \url{http://dx.doi.org/10.1214/13-BA815}.
#'
#' @references D.A. van Dyk, X.-L. Meng (2001):
#' ``The Art of Data Augmentation'',
#' \emph{Journal of Computational and Graphical Statistics}, 10(1), pages 1-50,
#' \url{http://dx.doi.org/10.1198/10618600152418584}.
#'
#' @seealso \code{\link{post.column.switch}} and \code{\link{post.sign.switch}}
#' for column switching and sign switching of the factor loading matrix and of
#' the correlation matrix of the latent factors to restore identification
#' \emph{a posteriori}.
#'
#' @seealso \code{\link{summary.befa}} and \code{\link{plot.befa}} to summarize
#' and plot the posterior results.
#'
#' @seealso \code{\link{simul.R.prior}} and \code{\link{simul.nfac.prior}} to
#' simulate the prior distribution of the correlation matrix of the factors and
#' the prior distribution of the indicator matrix, respectively. This is useful
#' to perform prior sensitivity analysis and to understand the role of the
#' corresponding parameters in the factor search.
#'
#' @examples
#' #### model without covariates
#'
#' set.seed(6)
#'
#' # generate fake data with 15 manifest variables and 3 factors
#' N <- 100    # number of observations
#' Y <- simul.dedic.facmod(N, dedic = rep(1:3, each = 5))
#'
#' # run MCMC sampler
#' # notice: 1000 MCMC iterations for illustration purposes only,
#' #  increase this number to obtain reliable posterior results!
#' mcmc <- befa(Y, Kmax = 5, iter = 1000)
#'
#' # post process MCMC draws to restore identification
#' mcmc <- post.column.switch(mcmc)
#' mcmc <- post.sign.switch(mcmc)
#' \donttest{
#' summary(mcmc)  # summarize posterior results
#' plot(mcmc)     # plot posterior results
#'
#' # summarize highest posterior probability (HPP) model
#' summary(mcmc, what = 'hppm')
#'
#' #### model with covariates
#'
#' # generate covariates and regression coefficients
#' Xcov <- cbind(1, matrix(rnorm(4*N), ncol = 4))
#' colnames(Xcov) <- c('(Intercept)', paste0('X', 1:4))
#' beta <- rbind(rnorm(15), rnorm(15), diag(3) %x% t(rnorm(5)))
#'
#' # add covariates to previous model
#' Y <- Y + Xcov %*% beta
#'
#' # specify model
#' model <- c('~ X1',                        # X1 covariate in all equations
#'            paste0('Y', 1:5,   ' ~ X2'),   # X2 covariate for Y1-Y5 only
#'            paste0('Y', 6:10,  ' ~ X3'),   # X3 covariate for Y6-Y10 only
#'            paste0('Y', 11:15, ' ~ X4'))   # X4 covariate for Y11-Y15 only
#' model <- lapply(model, as.formula)        # make list of formulas
#'
#' # run MCMC sampler, post process and summarize
#' mcmc <- befa(model, data = data.frame(Y, Xcov), Kmax = 5, iter = 1000)
#' mcmc <- post.column.switch(mcmc)
#' mcmc <- post.sign.switch(mcmc)
#' mcmc.sum <- summary(mcmc)
#' mcmc.sum
#'
#' # compare posterior mean of regression coefficients to true values
#' beta.comp <- cbind(beta[beta != 0], mcmc.sum$beta[, 'mean'])
#' colnames(beta.comp) <- c('true', 'mcmc')
#' print(beta.comp, digits = 3)
#' }
#'
#' @export befa
#' @import checkmate
#' @importFrom stats rnorm runif
#' @useDynLib BayesFM, .registration = TRUE, .fixes = "F_"

befa <- function(model, data, burnin = 1000, iter = 10000, Nid = 3, Kmax,
                 A0 = 10, B0 = 10, c0 = 2, C0 = 1, HW.prior = TRUE,
                 nu0 = Kmax + 1, S0 = 1, kappa0 = 2, xi0 = 1, kappa = 1/Kmax,
                 indp.tau0 = TRUE, rnd.step = TRUE, n.step = 5,
                 search.delay = min(burnin, 10), R.delay = min(burnin, 100),
                 dedic.start, alpha.start, sigma.start, beta.start, R.start,
                 verbose = TRUE)
{

  checkArgs <- makeAssertCollection()

  ##############################################################################
  ## data and model specification

  if (missing(data))
    data <- parent.frame()
  else
    assertDataFrame(data)

  if (is.matrix(model))
    model <- as.data.frame(model)

  if (is.data.frame(model)) {

    assertDataFrame(model, types = c('double', 'logical'), all.missing = FALSE)
    Ytype <- sapply(model, typeof)
    Yobs  <- as.matrix(model)
    Xobs  <- nX <- 0
    YXloc <- FALSE

  } else if (is.list(model) &
             all(sapply(model, is.formula))) {

    tmp <- extract.data(model, data)
    if (!is.null(tmp$errmsg)) {
      for (w in tmp$errmsg) checkArgs$push(w)
    } else {
      for (w in tmp$warnmsg) warning(w, immediate. = TRUE)
      Yobs  <- tmp$Yobs
      Ytype <- tmp$Ytype
      Xobs  <- tmp$Xobs
      Xlab  <- colnames(Xobs)
      YXloc <- tmp$YXloc
      nX    <- ncol(Xobs)
    }

  } else {

    checkArgs$push('Y should be a matrix, a data frame or a list of formulas.')

  }

  # if any errors, report and stop
  reportAssertions(checkArgs)

  # check manifest variables are either continuous or dichotomous
  Ycat <- rep(0, length(Ytype))
  Ycat[Ytype == 'logical'] <- 2
  Yind <- Ytype %in% c('double', 'numeric', 'logical')
  if (any(!Yind)) {
    checkArgs$push(paste('following variables not continuous nor dichotomous:',
                   paste0(Ylab[!Yind], collapse = ', ')))
  }

  Ylab  <- colnames(Yobs)
  nobs  <- nrow(Yobs)
  nvar  <- ncol(Yobs)
  nbeta <- sum(YXloc)
  Ymiss <- is.na(Yobs)
  Yobs[Ymiss] <- -99999  # flag for NA (not used in Fortran subroutine)


  ##############################################################################
  ## number of latent factors and identification restrictions

  # minimum number of dedicated variables per factor
  assertInt(Nid, lower = 1, upper = nvar, add = checkArgs)

  # check maximum number of latent factors and Ledermann bound
  Ledermann.bound <- 0.5 * (2 * nvar + 1 - sqrt(8 * nvar + 1))
  if (missing(Kmax)) {
    Kmax <- floor(min(nvar/Nid, Ledermann.bound))
  } else {
    assertInt(Kmax, lower = 1, upper = nvar, add = checkArgs)
  }
  if (Kmax > Ledermann.bound) {
    warning("Check identification! (Kmax exceeds Ledermann bound)",
            immediate. = TRUE)
  }

  # check consistency of Nid and Kmax
  if (Kmax > floor(nvar/Nid)) {
    msg <- paste("Too many latent factors specified given identification",
                 "restriction. Check arguments Nid and Kmax.")
    checkArgs$push(msg)
  }

  # throw warning in case of single-factor model
  if (Kmax == 1) {
    checkArgs$push("Single-factor model not allowed.")
  }

  reportAssertions(checkArgs)


  ##############################################################################
  ## prior parameters

  # function checking prior parameter values
  check.prior <- function(x, n, min, name) {
    if (length(x) == 1)
      x <- rep(x, n)
    assertNumeric(x, len = n, lower = min, finite = TRUE, any.missing = FALSE,
      .var.name = name, add = checkArgs)
    return(x)
  }

  tiny <- 10^-9
  A0 <- check.prior(A0, nvar, tiny, "A0")
  B0 <- check.prior(B0, nvar, tiny, "B0")
  c0 <- check.prior(c0, nvar, tiny, "c0")
  C0 <- check.prior(C0, nvar, tiny, "C0")
  S0 <- check.prior(S0, Kmax, tiny, "S0")
  nu0 <- check.prior(nu0, 1, Kmax, "nu0")
  xi0 <- check.prior(xi0, 1, tiny, "xi0")
  kappa <- check.prior(kappa, Kmax, tiny, "kappa")
  kappa0 <- check.prior(kappa0, 1, tiny, "kappa0")

  # use Huang-Wand (2013) prior?
  assertFlag(HW.prior, add = checkArgs)

  # use specific tau0 parameters for each manifest variable? [see CFSHP, p.36]
  assertFlag(indp.tau0, add = checkArgs)

  # if any errors, stop here
  reportAssertions(checkArgs)

  # prior values for Fortran subroutine
  prior.idiovar <- cbind(c0, C0)
  prior.loading <- 1/A0  # precision passed to Fortran subroutine
  prior.beta <- 1/B0  # precision passed to Fortran subroutine
  prior.dedic <- c(indp.tau0, 1/A0, c0, C0, xi0, kappa0, kappa)
  prior.facdist <- c(HW.prior, nu0, S0)


  ##############################################################################
  ## starting values

  ### idiosyncratic variances
  if (missing(sigma.start)) {
    sigma.start <- 1/rgamma(nvar, shape = c0, rate = C0)
  } else {
    assertNumeric(sigma.start, len = nvar, lower = tiny, any.missing = FALSE,
      add = checkArgs)
  }
  sigma.start[Ycat > 0] <- 1  # fix variance to 1 for binary variables

  ### factor loadings
  if (missing(alpha.start)) {
    alpha.start <- rnorm(nvar, mean = 0, sd = sqrt(A0))
  } else {
    assertNumeric(alpha.start, len = nvar, any.missing = FALSE,
                  add = checkArgs)
  }

  ### regression coefficients
  if (missing(beta.start)) {
    if (nbeta > 0) {
      beta.start <- rnorm(nbeta, mean = 0, sd = rep(sqrt(B0), rowSums(YXloc)))
    } else {
      beta.start <- double()
    }
  }
  assertNumeric(beta.start, len = nbeta, any.missing = FALSE, add = checkArgs)
  # prepare matrix to be passed to Fortran subroutine
  beta.start.1 <- matrix(-99999, nX, nvar)
  if (length(beta.start) == nbeta) {
    beta.start.1[t(YXloc)] <- beta.start
    beta.start.1 <- t(beta.start.1)
  }

  ### correlation matrix of latent factors
  if (missing(R.start)) {
    R.start <- diag(Kmax)
  }
  # check matrix is a correlation matrix, positive semi-definite, and invertible
  assertMatrix(R.start, mode = "double", nrows = Kmax, ncols = Kmax,
               any.missing = FALSE, add = checkArgs)
  if (!all(diag(R.start) == 1)) {
    checkArgs$push("R.start should be a correlation matrix.")
  }
  if (!is.pos.semidefinite.matrix(R.start)) {
    checkArgs$push("R.start should be a positive semi-definite matrix.")
  }
  if (!is.invertible.matrix(R.start)) {
    checkArgs$push("R.start is not invertible (singular matrix).")
  }

  ### indicators - default: maximum number of factors, random allocation
  if (missing(dedic.start)) {
    dedic.start <- rep(0, nvar)
    ind <- matrix(sample(Nid * Kmax), ncol = Kmax)
    for (k in 1:Kmax) dedic.start[ind[, k]] <- k
    dedic.start[dedic.start == 0] <- sample(Kmax, nvar - Nid * Kmax,
                                            replace = TRUE)
  }
  assertIntegerish(dedic.start, len = nvar, lower = 0, upper = Kmax,
                   any.missing = FALSE, add = checkArgs)
  # check identification constraint
  if (any(table(dedic.start[dedic.start != 0]) < Nid)) {
    checkArgs$push("dedic.start does not correspond to an identified model.")
  }

  ### latent factors
  start.factor <- replicate(Kmax, rnorm(nobs))


  ##############################################################################
  ## MCMC tuning

  assertCount(iter, positive = TRUE, add = checkArgs)
  assertCount(burnin, positive = FALSE, add = checkArgs)
  assertInt(search.delay, lower = 0, upper = burnin + iter)
  assertInt(R.delay, lower = 0, upper = burnin + iter)
  assertFlag(rnd.step, add = checkArgs)
  if (rnd.step) {
    assertNumber(n.step, lower = 1.1, add = checkArgs)
    n.step <- n.step - 1  # so that number of steps ~ 1 + Poisson(n.step-1)
  } else {
    assertCount(n.step, positive = TRUE, add = checkArgs)
  }


  ##############################################################################
  ## if any errors in arguments, report and stop execution

  assertFlag(verbose, add = checkArgs)
  reportAssertions(checkArgs)


  ##############################################################################
  ## MCMC sampling

  # total number of model parameters
  npar <- c(nvar, nvar, Kmax*(Kmax - 1)/2, nbeta)
  npar.all <- sum(npar)

  # seed for RNG in Fortran subroutine
  seed <- round(runif(1) * 10^9)

  # call Fortran subroutine
  if (verbose) cat("starting MCMC sampling...\n")
  mcmc <- .Fortran(F_befa,
                   as.integer(nvar),
                   as.integer(nobs),
                   as.integer(Kmax),
                   as.integer(Nid),
                   as.double(Yobs),
                   as.integer(Ycat),
                   as.logical(Ymiss),
                   as.integer(nX),
                   as.double(Xobs),
                   as.logical(YXloc),
                   as.integer(burnin + iter),
                   as.integer(burnin),
                   as.integer(search.delay),
                   as.integer(R.delay),
                   as.logical(rnd.step),
                   as.double(n.step),
                   as.integer(seed),
                   as.double(cbind(prior.loading, prior.idiovar)),
                   as.double(prior.beta),
                   as.double(prior.dedic),
                   as.double(prior.facdist),
                   as.double(cbind(alpha.start, sigma.start)),
                   as.double(beta.start.1),
                   as.integer(dedic.start),
                   as.double(start.factor),
                   as.double(R.start),
                   as.logical(verbose),
                   as.integer(npar.all),
                   MCMCdraws = double(iter * npar.all),
                   MCMCdedic = integer(iter * nvar),
                   MHacc = logical(iter))
  if (verbose) cat("done with sampling!\n")


  ##############################################################################
  ## label MCMC draws and return output

  # extract MCMC draws
  par.mcmc <- split(mcmc$MCMCdraws, rep(1:4, times = iter * npar))
  par.mcmc <- lapply(par.mcmc, matrix, nrow = iter)

  # label parameters
  names(par.mcmc)[1:3] <- c("alpha", "sigma", "R")
  colnames(par.mcmc$R) <- paste("R", rep(1:(Kmax-1), times = (Kmax-1):1),
                                 unlist(mapply(seq, 2:Kmax, Kmax)), sep = ":")
  colnames(par.mcmc$alpha) <- paste0("alpha:", Ylab)
  colnames(par.mcmc$sigma) <- paste0("sigma:", Ylab)

  iter.lab <- burnin + 1:iter
  rownames(par.mcmc$alpha) <- iter.lab
  rownames(par.mcmc$sigma) <- iter.lab
  rownames(par.mcmc$R)     <- iter.lab

  if (nbeta > 0) {
    names(par.mcmc)[4] <- "beta"
    beta.lab <- c()
    for (i in 1:nvar) {
      if (!any(YXloc[i, ])) next
      beta.lab <- c(beta.lab, paste(Ylab[i], Xlab[YXloc[i, ]], sep = ":"))
    }
    colnames(par.mcmc$beta) <- beta.lab
    rownames(par.mcmc$beta) <- iter.lab
  }

  # indicators
  dedic.mcmc <- as.integer(mcmc$MCMCdedic)
  dedic.mcmc <- matrix(dedic.mcmc, nrow = iter, ncol = nvar)
  colnames(dedic.mcmc) <- paste0("dedic:", Ylab)
  rownames(dedic.mcmc) <- iter.lab

  # number of active latent factors across MCMC iterations
  nfac.mcmc <- apply(dedic.mcmc, 1, count.unique.nonzero)

  # prepare and return output
  output       <- par.mcmc
  output$dedic <- dedic.mcmc
  output$nfac  <- nfac.mcmc
  output$MHacc <- mcmc$MHacc
  attr(output, "call")   <- match.call()
  attr(output, "title")  <- "BEFA posterior sample"
  attr(output, "Kmax")   <- Kmax
  attr(output, "Nid")    <- Nid
  attr(output, "iter")   <- iter
  attr(output, "burnin") <- burnin
  attr(output, "post.column.switch") <- FALSE
  attr(output, "post.sign.switch")   <- FALSE
  class(output) <- "befa"
  return(output)

}
