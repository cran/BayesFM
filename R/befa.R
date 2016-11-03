#'
#' Bayesian Exploratory Factor Analysis
#'
#' This function implements the Bayesian Exploratory Factor Analysis
#' (\code{befa}) developed in Conti et al. (CFSHP, 2014). It runs a MCMC sampler
#' for a factor model with dedicated factors, where each manifest variable is
#' allowed to load on at most one latent factor. The allocation of the manifest
#' variables to the latent factors is not fixed \emph{a priori} but determined
#' stochastically during sampling. The minimum number of variables dedicated to
#' each factor can be controlled by the user to achieve the desired level of
#' identification. The manifest variables can be continuous or dichotomous, and
#' control variables can be introduced as covariates.
#'
#' @param model
#'        This argument can be specified in two different ways:
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
#'        Binary variables should be specified as logical vectors in the data
#'        frame to be treated as dichotomous. \code{NA} values are accepted in
#'        manifest variables only.
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
#' @param loading.start
#'        Starting values for the factor loadings. Numeric vector of length
#'        equal to the number of manifest variables. If missing, sampled from a
#'        Normal distribution with zero mean and variance \code{A0}.
#' @param dedic.start
#'        Starting values for the indicators. Vector of integers of length equal
#'        to the number of manifest variables. Each element takes a value among
#'        0, 1, ..., \code{Kmax}. If missing, random allocation of the manifest
#'        variables to the maximum number of factors \code{Kmax}, with a minimum
#'        of \code{Nid} manifest variables loading on each factor.
#' @param idiovar.start
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
#'        If \code{TRUE}, display information such as the progression of the
#'        MCMC sampler.
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
#' regression coefficients \eqn{\beta}.
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
#' \strong{Prior specification.}
#' The indicators are assumed to have the following probabilities,
#' for \eqn{k = 1, ..., K}:
#' \deqn{Prob(\Delta_m = e_k \mid \tau_k) = \tau_k}{
#'       Prob(\Delta_m = e_k | \tau_k) = \tau_k}
#' \deqn{\tau = (\tau_0, \tau_1, ..., \tau_K)}
#' where the probabilities, if \code{indp.tau0 = FALSE}, are specified as:
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
#' of formulas as \code{model} (see example below).
#'
#' To sample the correlation matrix \eqn{R} of the latent factors, marginal data
#' augmentation (van Dyk and Meng, 2001) is implemented, see CFSHP section 2.2.
#' Using the transformation \eqn{\Omega = \Lambda^{1/2} R \Lambda^{1/2}}, the
#' parameters \eqn{\Lambda = diag(\lambda_1, ..., \lambda_K)} are used as
#' \emph{working parameters}. These parameters correspond to the variances of
#' the latent factors in an expanded version of the model where the factors do
#' not have unit variances. Two prior distributions can be specified on the
#' covariance matrix \eqn{\Omega} in the expanded model:
#' \itemize{
#'   \item Inverse-Wishart distribution if \code{HW.prior = FALSE}:
#'   \deqn{\Omega \sim \mathcal{I}nv-\mathcal{W}ishart(\nu_0, diag(S_0))}{
#'         \Omega ~ Inv-Wishart(\nu_0, diag(S0))}
#'   with \eqn{\nu_0} = \code{nu0} and \eqn{S_0} = \code{S0}.
#'   \item Huang-Wand (2013) prior if \code{HW.prior = TRUE}:
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
#' @return The function returns an object of class '\code{befa}' containing:
#' \itemize{
#'   \item \code{draws}: A matrix containing the MCMC draws of the model
#'   parameters in its columns, stored in the following order: factor loadings,
#'   idiosyncratic variances, regression coefficients (if any), off-diagonal
#'   elements of the correlation matrix of the factors.
#'   \item \code{dedic}: A matrix containing the MCMC draws of the indicators.
#'   \item \code{MHacc}: Acceptance rate of the Metropolis-Hastings algorithm.
#' }
#' The parameters \code{Kmax} and \code{Nid} are saved as object attributes.
#'
#' Note that identification is achieved only with respect to the scale of the
#' latent factors. Non-identifiability problems may affect the posterior sample
#' because of column switching and sign switching of the factor loadings.
#' These issues can be addressed \emph{a posteriori} with the functions
#' \code{\link{post.column.switch}} and \code{\link{post.sign.switch}}.
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
#' @seealso \code{\link{HPPmodel}} to find the highest posterior probability
#' model.
#'
#' @examples
#' #### model without covariates
#'
#' set.seed(6)
#'
#' # generate fake data with 15 manifest variables and 3 factors
#' N <- 200    # number of observations
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
#'
#' # summarize highest posterior probability (HPP) model
#' hppm <- HPPmodel(mcmc)
#' hppm$prob                                # HPP model probability
#' hppm$dedic                               # indicators in HPP model
#' colMeans(hppm$draws)                     # posterior means in HPP model
#' attributes(Y)[c('alpha', 'sigma', 'R')]  # true model parameters
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
#' model <- c('~ X1',                           # X1 covariate in all equations
#'            paste0('Y', 1:5,   ' ~ X2'),      # X2 covariate for Y1-Y5 only
#'            paste0('Y', 6:10,  ' ~ X3'),      # X3 covariate for Y6-Y10 only
#'            paste0('Y', 11:15, ' ~ X4'))      # X4 covariate for Y11-Y15 only
#' model <- lapply(as.list(model), as.formula)  # make list of formulas
#'
#' # run MCMC sampler
#' mcmc <- befa(model, data = data.frame(Y, Xcov), Kmax = 5, iter = 1000)
#'
#' # compare posterior mean of regression coefficients to true values
#' post <- cbind(beta[beta != 0], colMeans(mcmc$draws)[31:75])
#' colnames(post) <- c('true', 'mcmc')
#' post
#'
#' @export befa
#' @import checkmate
#' @useDynLib BayesFM
#' @importFrom stats is.empty.model model.frame terms
#' @importFrom stats rWishart rgamma rnorm runif

befa <- function(model, data, burnin = 1000, iter = 10000, Nid = 3, Kmax,
                 A0 = 10, B0 = 10, c0 = 2, C0 = 1, HW.prior = TRUE,
                 nu0 = Kmax + 1, S0 = 1, kappa0 = 2, xi0 = 1, kappa = 1/Kmax,
                 indp.tau0 = TRUE, rnd.step = TRUE, n.step = 5,
                 search.delay = min(burnin, 10), R.delay = min(burnin, 100),
                 dedic.start, loading.start, idiovar.start, beta.start, R.start,
                 verbose = FALSE)
{

  checkArgs <- makeAssertCollection()

  ##############################################################################
  ## data and model specification

  if (missing(data))
    data <- parent.frame() else assertDataFrame(data)

  if (is.matrix(model))
    model <- as.data.frame(model)

  if (is.data.frame(model)) {
    assertDataFrame(model, types = c("double", "logical"), all.missing = FALSE)
    Ycat <- sapply(model, typeof)
    Yobs <- as.matrix(model)
    Xobs <- nX <- 0
    YXloc <- FALSE
  } else if (is.list(model) & all(sapply(model, is.formula))) {
    all.var <- unique(unlist(lapply(model, all.vars)))
    ind <- all.var %in% names(data)
    if (!all(ind)) {
      checkArgs$push(paste0("Following variable(s) not in data frame: ",
                            paste0(all.var[!ind],
        collapse = ", "), "."))
    } else {
      tmp <- extract.data(model, data)
      for (w in tmp$warnmsg) warning(w, immediate. = TRUE)
      for (w in tmp$errmsg) checkArgs$push(w)
      Yobs <- tmp$Yobs
      Ycat <- tmp$Ycat
      Xobs <- tmp$Xobs
      Xlab <- colnames(Xobs)
      YXloc <- tmp$YXloc
      nX <- ncol(Xobs)
    }
  } else {
    checkArgs$push("Y should be a matrix, a data frame or a list of formulas.")
  }

  # if any errors, report and stop
  reportAssertions(checkArgs)

  # check manifest variables are either continuous or dichotomous
  Ycat <- ifelse(Ycat == "double" | Ycat == "numeric", 0,
                 ifelse(Ycat == "logical", 2, 1))
  Yind <- Ycat == 0 | Ycat == 2
  if (any(!Yind)) {
    checkArgs$push("Following variables not continuous nor dichotomous:",
      paste0(Ylab[!Yind], collapse = ", "), ".")
  }

  Ylab <- colnames(Yobs)
  nobs <- nrow(Yobs)
  nmeas <- ncol(Yobs)
  nbeta <- sum(YXloc)
  Ymiss <- is.na(Yobs)
  Yobs[Ymiss] <- -99999  # flag for NA (not used in Fortran subroutine)


  ##############################################################################
  ## number of latent factors and identification restrictions

  # minimum number of dedicated variables per factor
  assertInt(Nid, lower = 1, upper = nmeas, add = checkArgs)

  # check maximum number of latent factors and Ledermann bound
  Ledermann.bound <- 0.5 * (2 * nmeas + 1 - sqrt(8 * nmeas + 1))
  if (missing(Kmax)) {
    Kmax <- floor(min(nmeas/Nid, Ledermann.bound))
  } else {
    assertInt(Kmax, lower = 1, upper = nmeas, add = checkArgs)
  }
  if (Kmax > Ledermann.bound) {
    warning("Check identification! (Kmax exceeds Ledermann bound)",
            immediate. = TRUE)
  }

  # check consistency of Nid and Kmax
  if (Kmax > floor(nmeas/Nid)) {
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
  A0 <- check.prior(A0, nmeas, tiny, "A0")
  B0 <- check.prior(B0, nmeas, tiny, "B0")
  c0 <- check.prior(c0, nmeas, tiny, "c0")
  C0 <- check.prior(C0, nmeas, tiny, "C0")
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
  if (missing(idiovar.start)) {
    idiovar.start <- 1/rgamma(nmeas, shape = c0, rate = C0)
  } else {
    assertNumeric(idiovar.start, len = nmeas, lower = tiny, any.missing = FALSE,
      add = checkArgs)
  }
  idiovar.start[Ycat > 0] <- 1  # fix variance to 1 for binary variables

  ### factor loadings
  if (missing(loading.start)) {
    loading.start <- rnorm(nmeas, mean = 0, sd = sqrt(A0))
  } else {
    assertNumeric(loading.start, len = nmeas, any.missing = FALSE,
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
  beta.start.1 <- matrix(-99999, nX, nmeas)
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
    dedic.start <- rep(0, nmeas)
    ind <- matrix(sample(Nid * Kmax), ncol = Kmax)
    for (k in 1:Kmax) dedic.start[ind[, k]] <- k
    dedic.start[dedic.start == 0] <- sample(Kmax, nmeas - Nid * Kmax,
                                            replace = TRUE)
  }
  assertIntegerish(dedic.start, len = nmeas, lower = 0, upper = Kmax,
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
  npar <- 2 * nmeas + nbeta + Kmax * (Kmax - 1)/2

  # seed for RNG in Fortran subroutine
  seed <- round(runif(1) * 10^9)

  # call Fortran subroutine
  if (verbose) cat("Starting MCMC sampling...\n")
  mcmc <- .Fortran("befa",
                   as.integer(nmeas),
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
                   as.double(cbind(loading.start, idiovar.start)),
                   as.double(beta.start.1),
                   as.integer(dedic.start),
                   as.double(start.factor),
                   as.double(R.start),
                   as.logical(verbose),
                   as.integer(npar),
                   MCMCdraws = double(iter * npar),
                   MCMCgroup = integer(iter * nmeas),
                   MHacc = logical(iter),
                   PACKAGE = "BayesFM")
  if (verbose) cat("done with sampling!\n")


  ##############################################################################
  ## label MCMC draws and return

  draws <- matrix(mcmc$MCMCdraws, nrow = iter, ncol = npar)

  alpha.lab <- paste0("alpha:", 1:nmeas)
  var.lab <- paste0("sigma:", 1:nmeas)
  beta.lab <- c()
  if (nbeta > 0) {
    for (i in 1:nmeas) {
      if (!any(YXloc[i, ]))
        next
      beta.lab <- c(beta.lab, paste(Ylab[i], Xlab[YXloc[i, ]], sep = ":"))
    }
  }
  ind <- matrix(1:Kmax, Kmax, Kmax)
  li <- lower.tri(ind)
  R.lab <- paste("R", t(ind)[li], ind[li], sep = ":")

  rownames(draws) <- seq(burnin + 1, burnin + iter)
  colnames(draws) <- c(alpha.lab, var.lab, beta.lab, R.lab)

  # indicators
  dedic <- as.integer(mcmc$MCMCgroup)
  dedic <- matrix(dedic, nrow = iter, ncol = nmeas)
  colnames(dedic) <- Ylab
  rownames(dedic) <- seq(burnin + 1, burnin + iter)

  # prepare and return output
  output <- list(draws = draws, dedic = dedic, MHacc = mean(mcmc$MHacc))
  attr(output, "title") <- "BEFA posterior sample"
  attr(output, "Kmax") <- Kmax
  attr(output, "Nid") <- Nid
  attr(output, "post.column.switch") <- FALSE
  attr(output, "post.sign.switch") <- FALSE
  class(output) <- "befa"
  return(output)

}
