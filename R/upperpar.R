upperpar <- function(object, ...) UseMethod("upperpar")

print.upperpar <- function(x, digits = max(3, getOption("digits") - 3),
  ...)
{
  cat("Item response upper asymptote parameters (", attr(x, "model"),
    "):\n", sep = "")
  print(coef(x), digits = digits, ...)
  invisible(x)
}

coef.upperpar <- function(object, ...) {
  lbs <- names(object)
  object <- as.vector(object)
  names(object) <- lbs
  return(object)
}

vcov.upperpar <- function(object, ...) {
  attr(object, "vcov")
}
