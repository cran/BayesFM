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