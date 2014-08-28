\name{piplot}
\alias{piplot}

\title{Person-Item Plots for IRT Models}

\description{
  Base graphics plotting function for person-item plot visualization of IRT models.
}

\usage{
  piplot(object, ref = NULL, items = NULL, xlim = NULL, names = NULL,
    labels = TRUE, main = "Person-Item Plot", xlab = "Latent trait",
    abbreviate = FALSE, cex.axis = 0.8, cex.text = 0.5, cex.points = 1.5, ...)
}

\arguments{
  \item{object}{a fitted model object of class
    \code{"raschmodel"}, \code{"rsmodel"} or \code{"pcmodel"}.}
  \item{ref}{argument passed over to internal calls of
    \code{\link{threshpar}} and \code{\link{itempar}}.}
  \item{items}{character or numeric, specifying the items which
    should be visualized in the person-item plot.}
  \item{xlim}{numeric, specifying the x axis limits.}
  \item{names}{character, specifying labels for the items.}
  \item{labels}{logical, whether to draw the number of the threshold as
    text below the threshold.}
  \item{main}{character, specifying the overall title of the plot.}
  \item{xlab}{character, specifying the x axis labels.}
  \item{abbreviate}{logical or numeric, specifying whether object names
    are to be abbreviated. If numeric this controls the length of the abbreviation.}
  \item{cex.axis}{numeric, the magnification to be used for the axis
    notation relative to the current setting of \code{cex}.}
  \item{cex.text}{numeric, the magnification to be used for the symbols
    relative to the current setting of \code{cex}.}
  \item{cex.points}{numeric, the magnification to be used for the points
    relative to the current setting of \code{cex}.}
  \item{\dots}{further arguments passed to internal calls of
    \code{\link{lines}}, \code{\link{points}} and \code{\link{text}}}
}

\details{
  The person-item plot visualization illustrates the distribution of the
  person parameters against the absolute item threshold parameters under
  a certain data set and IRT model.
}

\seealso{\code{\link{curveplot}}, \code{\link{regionplot}},
  \code{\link{profileplot}}, \code{\link{infoplot}}}


\examples{
## Verbal aggression data
data("VerbalAggression", package = "psychotools")

## Fit a partial credit model to VerbalAggression data
pcmod <- pcmodel(VerbalAggression$resp)

## Create a person-item plot visualization of the fitted PCM
plot(pcmod, type = "piplot")

## Just visualize the first six items
plot(pcmod, type = "piplot", ref = 1:6)
}

\keyword{aplot}
