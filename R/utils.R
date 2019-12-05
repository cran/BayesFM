
### check if object is a formula

is.formula <- function(x) return(class(x) == "formula")


### check if matrix is positive semidefinite

is.pos.semidefinite.matrix <- function(x) {
  if (!is.matrix(x))
    return(FALSE)
  if (!isSymmetric(x, tol = 10^-9))
    return(FALSE)
  return(all(eigen(x, only.values = TRUE)$values >= 0))
}


### check if matrix can be inverted

is.invertible.matrix <- function(x) {
  return(class(try(solve(x), silent = TRUE))[1L] == "matrix")
}


### relabel indicators for dedicated factor model
### example:
###   > d <- c(6, 6, 6, 2, 2, 2, 0, 0, 4, 4, 4)
###   > relabel(d) [1] 1 1 1 2 2 2 0 0 3 3 3

relabel.dedic <- function(d) {
  u <- unique(d[d != 0])
  t <- rep(0, max(d))
  t[u] <- 1:length(u)
  d[d != 0] <- t[d]
  return(d)
}


### count number of unique nonzero elements in x

count.unique.nonzero <- function(x) return(length(unique(x[x != 0])))


### convert matrix to list

mat2list <- function(x, byrow = FALSE) {
  if (byrow) {
    y <- split(x, c(row(x)))
    names(y) <- rownames(x)
  } else {
    y <- split(x, c(col(x)))
    names(y) <- colnames(x)
  }
  return(y)
}


### get maximum correlation (in absolute value) from correlation matrix

get.maxcor <- function(x) max(abs(x[lower.tri(x)]))


### get minimum eigenvalue of symmetric matrix

get.mineig <- function(x) min(eigen(x, symmetric = TRUE,
                                    only.values = TRUE)$values)
