#' @importFrom stats na.omit cor

extract.data <- function(model, data)
{

  errmsg <- warnmsg <- NULL

  eqIng <- list()
  i <- 1
  for (eq in model) {
    eqf <- model.frame(eq, data, na.action = NULL)
    eqt <- terms(eqf)
    eqIng[[i]] <- list(dat = eqf)
    if (attr(eqt, "response") == 1) {
      eqIng[[i]]$Ylab <- colnames(eqf)[1]
      eqIng[[i]]$Ycat <- attr(eqt, "dataClasses")[[1]]
      eqIng[[i]]$Xlab <- colnames(eqf)[-1]
    } else {
      eqIng[[i]]$Ylab <- NA
      eqIng[[i]]$Ycat <- NA
      eqIng[[i]]$Xlab <- colnames(eqf)
    }
    eqIng[[i]]$cst <- attr(eqt, "intercept") == 1
    i <- i + 1
  }

  # prepare objects to be returned

  Ylab <- na.omit(sapply(eqIng, function(x) return(x$Ylab)))
  Xlab <- lapply(eqIng, function(x) return(x$Xlab))
  Xlab <- unique(unlist(Xlab))

  YXdat <- lapply(eqIng, function(x) return(x$dat))
  YXdat <- do.call("cbind", YXdat)
  YXdat <- YXdat[!duplicated(t(YXdat))]

  Yobs <- as.matrix(YXdat[Ylab])
  Xobs <- as.matrix(YXdat[Xlab])
  Xobs <- cbind("(Intercept)" = 1, Xobs)
  Xlab <- colnames(Xobs)
  nmeas <- length(Ylab)
  nX <- length(Xlab)

  YXloc <- matrix(FALSE, nmeas, nX, dimnames = list(Ylab, Xlab))
  YXloc[, 1] <- TRUE  # intercept terms by default
  Ycat <- rep(NA, nmeas)
  names(Ycat) <- Ylab
  for (eq in eqIng) {
    if (is.na(eq$Ylab)) {
      YXloc[, eq$Xlab] <- TRUE
      # YXloc[, '(Intercept)'] <- eq$cst
    } else {
      YXloc[eq$Ylab, eq$Xlab] <- TRUE
      YXloc[eq$Ylab, "(Intercept)"] <- eq$cst
      Ycat[eq$Ylab] <- eq$Ycat
    }
  }

  # missing values in covariates
  if (any(is.na(Xobs))) {
    nm <- nrow(Yobs)
    Xobs <- na.omit(Xobs)
    Yobs <- Yobs[rownames(Xobs), ]
    nm <- nm - nrow(Yobs)
    warnmsg <- c(warnmsg, paste0("Observations with NA values in at least ",
      "one covariate were discarded (", nm, " cases)."))
  }

  # check for multicollinearity
  Xcor <- cor(Xobs[, -1])  # ignore vector of ones for intercept
  for (i in 1:nmeas) {
    Xcori <- Xcor[YXloc[i, -1], YXloc[i, -1]]
    Xcori <- Xcori[lower.tri(Xcori)]
    if (any(abs(Xcori - 1) < 1e-12)) {
      errmsg <- c(errmsg, paste0("Perfect multicollinearity between covariates",
        " of measurement ", Ylab[i], "."))
    } else if (any(abs(Xcori) > 0.95)) {
      warnmsg <- c(warnmsg, paste0("Possible multicollinearity problem ",
                                   "between covariates for measurement ",
                                   Ylab[i], "."))
    }
  }

  # return objects
  return(list(Yobs = Yobs, Ycat = Ycat, Xobs = Xobs, YXloc = YXloc,
              errmsg = errmsg, warnmsg = warnmsg))

}
