guesspar <- function(object, ...) UseMethod("guesspar")

print.guesspar <- function(x, digits = max(3, getOption("digits") - 3),
  ...)
{
  cat("Item response guessing parameters (", attr(x, "model"), "):\n",
    sep = "")
  print(coef(x), digits = digits, ...)
  invisible(x)
}

coef.guesspar <- function(object, ...) {
  lbs <- names(object)
  object <- as.vector(object)
  names(object) <- lbs
  return(object)
}

vcov.guesspar <- function(object, ...) {
  attr(object, "vcov")
}
