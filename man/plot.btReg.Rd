\name{plot.btReg}
\alias{plot.btReg}

\title{Visualizing Bradley-Terry Models}

\description{
  Base graphics plotting function for Bradley-Terry models.
}

\usage{
\S3method{plot}{btReg}(x, worth = TRUE, index = TRUE, names = TRUE,
  ref = TRUE, abbreviate = FALSE, type = NULL, lty = NULL,
  xlab = "Objects", ylab = NULL, \dots)
}

\arguments{
  \item{x}{an object of class \code{"btReg"}.}
  \item{worth}{logical. Should worth parameters (or alternatively coefficients
    on log-scale) be displayed?}
  \item{index}{logical. Should different indexes for different items be used?}
  \item{names}{logical. Should the names for the objects be displayed?}
  \item{ref}{logical. Should a horizontal line for the reference level be drawn?}
  \item{abbreviate}{logical or numeric. Should object names be abbreviated?
    If numeric this controls the length of the abbreviation.}
  \item{type}{plot type. Default is \code{"b"} if \code{index} is \code{TRUE}.}
  \item{lty}{line type.}
  \item{xlab, ylab}{x and y axis labels.}
  \item{\dots}{further arguments passed to \code{\link{plot}}.}
}

\seealso{\code{\link{btReg.fit}}}

\examples{
## data
data("GermanParties2009", package = "psychotools")

## Bradley-Terry model
bt <- btReg.fit(GermanParties2009$preference)
plot(bt)
plot(bt, worth = FALSE)
plot(bt, index = FALSE)
}

\keyword{hplot}
