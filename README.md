# BayesFM: Bayesian Inference for Factor Modeling

<!-- badges: start -->
[![CRAN](https://www.r-pkg.org/badges/version/BayesFM)](https://cran.r-project.org/package=BayesFM)
<!-- badges: end -->

## Installation

This package can be installed in R using:
```{r}
install.packages("BayesFM")
```

Alternatively, it can be installed directly from Github:
```
# install.packages("devtools")
devtools::install_github("piatek/BayesFM")
```

### Notes:

- This package requires a Fortran compiler and a C compiler.
- GNU Fortran version 4.6.3 or later is recommended.
- F95 and earlier versions of GNU Fortran may not work because
  of unsupported Fortran 2003 features used in this package.

Windows users:

- It is recommended to install 
  [Rtools](https://cran.r-project.org/bin/windows/Rtools/)
  to compile this package from source, see
  [R for Windows FAQ](https://cran.r-project.org/bin/windows/base/rw-FAQ.html).
