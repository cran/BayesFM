# Changes in version 0.1.7
- Refactored Fortran code: Removed polymorphism to allow compilation with flang 
  versions 17 and 18.

# Changes in version 0.1.6
- Changed extension of Fortran files (.f95 -> .f90) to accommodate Intel Fortran 
  compiler.

# Changes in version 0.1.5
- Added missing AC_OUTPUT to configure.ac to remove warning on CRAN.
- Minor # Changes to documentation (now using \doi{} for references).

# Changes in version 0.1.4
- Updated maintainer's email address.

# Changes in version 0.1.3
- Fixed error occurring when "class(matrix(...))" returns a vector of length 
  two.

# Changes in version 0.1.2
- Patched "configure.ac" file to solve compilation problems on some platforms.
- Fortran native routines now registered and symbol search disabled.
- README file added to provide information about installation requirements.
- Added empty "configure.win" file for Windows users.

# Changes in version 0.1.1
- Added simul.R.prior() and simul.nfac.prior() functions to simulate the prior
  distributions of the number of latent factors, and of the correlation matrix
  of the factors.
- Improved plot() generic function for "befa" object to show more posterior
  results, such as heatmaps for indicator posterior probabilities, factor
  loadings and correlation matrix of the latent factors.
- Now using ggplot2 package for nice-looking graphs.
- Added "configure" script to address cross-compilation problems.
- Minor improvements to documentation.

# Changes in version 0.1.0
- Fixed dependencies in Fortran code (Makevars file) and added gfortran as
  system requirement in DESCRIPTION file to solve cross-compilation problems.
- Added print(), summary() and plot() functions to "befa" object.
- Removed HPPmodel(), as highest posterior probability models can now be
  summarized with generic function summary().
- Renamed two arguments passed to befa() function for consistency
  (loading.start -> alpha.start, idiovar.start -> sigma.start).
- Removed deprecated cleanup file.
- Added NEWS file.

# Changes in version 0.0.2
- First version released on CRAN.
