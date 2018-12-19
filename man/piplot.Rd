\name{piplot}
\alias{piplot}

\title{Person-Item Plots for IRT Models}

\description{
  Base graphics plotting function for person-item plot visualization of IRT models.
}

\usage{
  piplot(object, pcol = NULL, histogram = TRUE, ref = NULL, items = NULL,
    xlim = NULL, names = NULL, labels = TRUE, main = "Person-Item Plot",
    xlab = "Latent trait", abbreviate = FALSE, cex.axis = 0.8, cex.text = 0.5,
    cex.points = 1.5, grid = TRUE, \dots)
}

\arguments{
  \item{object}{a fitted model object of class \code{"raschmodel"},
    \code{"rsmodel"}, \code{"pcmodel"}, \code{"plmodel"} or \code{"gpcmodel"}.}
  \item{pcol}{optional character (vector), specifying the color(s) used for the
    person parameter plot.}
  \item{histogram}{logical. For models estimated via MML (\code{plmodel}s and
    \code{gpcmodel}s), should a histogram of the person-wise (individual) person
    parameters be drawn additionally to the normal distribution density of the
    person parameters?}
  \item{ref}{argument passed over to internal calls of \code{\link{threshpar}}
    and \code{\link{itempar}}.}
  \item{items}{character or numeric, specifying the items which should be
    visualized in the person-item plot.}
  \item{xlim}{numeric, specifying the x axis limits.}
  \item{names}{character, specifying labels for the items.}
  \item{labels}{logical, whether to draw the number of the threshold as text
    below the threshold.}
  \item{main}{character, specifying the overall title of the plot.}
  \item{xlab}{character, specifying the x axis labels.}
  \item{abbreviate}{logical or numeric, specifying whether object names are to
    be abbreviated. If numeric, this controls the length of the abbreviation.}
  \item{cex.axis}{numeric, the magnification to be used for the axis notation
    relative to the current setting of \code{cex}.}
  \item{cex.text}{numeric, the magnification to be used for the symbols
    relative to the current setting of \code{cex}.}
  \item{cex.points}{numeric, the magnification to be used for the points
    relative to the current setting of \code{cex}.}
  \item{grid}{logical or color specification of horizontal grid lines. If set to
    \code{FALSE} or \code{"transparent"} grid lines can be suppressed.}
  \item{\dots}{further arguments passed to internal calls of
    \code{\link{lines}}, \code{\link{points}} and \code{\link{text}}.}
}

\details{
  The person-item plot visualization illustrates the distribution of the person
  parameters against the absolute item threshold parameters under a certain data
  set and IRT model. For models estimated via MML (\code{plmodel}s and
  \code{gpcmodel}s), the normal distribution density of the person parameters is
  drawn. If \code{histogram} is set to \code{TRUE} (the default), a histogram of
  the person-wise (individual) person parameters is drawn additionally. If a
  multiple group model has been fitted by supplying an \code{impact} variable,
  multiple person parameter plots are drawn, each corresponding to a specific
  level of this variable.
}

\seealso{\code{\link{curveplot}}, \code{\link{regionplot}},
  \code{\link{profileplot}}, \code{\link{infoplot}}}

\examples{
## load verbal agression data
data("VerbalAggression", package = "psychotools")

## fit partial credit model to verbal aggression data
pcmod <- pcmodel(VerbalAggression$resp)

## create a person-item plot visualization of the fitted PCM
plot(pcmod, type = "piplot")

## just visualize the first six items and the person parameter plot
plot(pcmod, type = "piplot", items = 1:6, pcol = "lightblue")

if(requireNamespace("mirt")) {
## fit generalized partial credit model to verbal aggression data
gpcmod <- gpcmodel(VerbalAggression$resp)

## create a person-item plot visualization of the fitted GPCM
plot(gpcmod, type = "piplot")

## turn off the histogram and grid
plot(gpcmod, type = "piplot", histogram = FALSE, grid = FALSE)

## fit GPCM to verbal aggression data accounting for gender impact
mgpcmod <- gpcmodel(VerbalAggression$resp, impact = VerbalAggression$gender)

## create a person-item plot visualization of the fitted GPCM
plot(mgpcmod, type = "piplot", pcol = c("darkgreen", "darkorange"))
}
}

\keyword{aplot}

