.onLoad <- function(libname = find.package("mr.raps"), pkgname = "mr.raps"){

  # CRAN Note avoidance
  if (getRversion() >= "3.1.0")
    utils::globalVariables(c("w", "pval.selection", "beta.hat", "p", "w", "shape"))
  invisible()
}
