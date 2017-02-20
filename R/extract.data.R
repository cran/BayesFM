#' @importFrom stats complete.cases cor terms

extract.data <- function(model, data)
{

  errmsg <- warnmsg <- NULL

  # get equation terms
  neq <- length(model)
  all.terms <- lapply(model, terms)

  # which equations have a manifest variable?
  resp <- sapply(all.terms, attr, 'response')
  resp <- as.logical(resp)

  # manifest variables
  allY <- sapply(lapply(model, all.vars), '[', 1)
  Ylab <- unique(allY[resp])
  nY   <- length(Ylab)
  allY[!resp] <- NA

  # intercept terms
  const <- sapply(all.terms, attr, 'intercept')
  const <- as.logical(const)

  # covariates
  allX <- lapply(all.terms, attr, 'term.labels')
  Xlab <- unique(unlist(allX))
  Xlab <- c('(Intercept)', Xlab)
  nX   <- length(Xlab)

  YXloc <- matrix(FALSE, nY, nX, dimnames = list(Ylab, Xlab))
  for (i in 1:neq) {
    if (is.na(allY[i])) {
      YXloc[, unlist(allX[i])] <- TRUE
      YXloc[, '(Intercept)'] <- const[i]
    } else {
      YXloc[allY[i], unlist(allX[i])] <- TRUE
      YXloc[allY[i], '(Intercept)'] <- const[i]
    }
  }

  # check that all variables are in data frame
  YXlab <- unique(unlist(c(Ylab, Xlab)))
  YXlab <- YXlab[YXlab != '(Intercept)']   # no error message for intercept term
  indat <- YXlab %in% names(data)
  if (any(!indat)) {
    errmsg <- paste('following variables not in data:\n   ',
                    paste(YXlab[!indat], collapse = ', '))
    return(list(errmsg))
  }
  YXdata <- data[YXlab]
  nobs   <- nrow(YXdata)
  YXdata[["(Intercept)"]] <- 1    # add vector of ones for intercept terms

  # type of manifest variables
  Ytype <- sapply(YXdata[Ylab], typeof)

  # discard missing values in covariates
  nomiss <- complete.cases(YXdata[Xlab])
  YXdata <- YXdata[nomiss, ]
  if (!all(nomiss)) {
    warnmsg <- paste(sum(!nomiss), 'observations discarded because of NAs',
                     'in at least one covariate')
    nobs <- nobs - sum(!nomiss)
  }

  # check for multicollinearity in covariates specified in each equation
  Xcor <- cor(data[Xlab[Xlab != '(Intercept)']])  # exclude intercept
  for (i in 1:nY) {
    Xcori <- Xcor[YXloc[i, -1], YXloc[i, -1]]
    Xcori <- Xcori[lower.tri(Xcori)]
    if (any(abs(Xcori - 1) < 1e-12)) {
      errmsg <- c(errmsg, paste('perfect multicollinearity between covariates',
                                 'of manifest variable', Ylab[i]))
    } else if (any(abs(Xcori) > 0.95)) {
      warnmsg <- c(warnmsg, paste('possible multicollinearity problem between',
                                  'covariates of manifest variable', Ylab[i]))
    }
  }

  # return
  return(list(Ytype   = Ytype,
              Yobs    = as.matrix(YXdata[Ylab]),
              Xobs    = as.matrix(YXdata[Xlab]),
              YXloc   = YXloc,
              errmsg  = errmsg,
              warnmsg = warnmsg))

}
