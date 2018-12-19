\name{curveplot}
\alias{curveplot}

\title{Response Curve Plots for IRT Models}

\description{
  Base graphics plotting function for response curve plot visualization of IRT
  models.
}

\usage{
  curveplot(object, ref = NULL, items = NULL, names = NULL,
    layout = NULL, xlim = NULL, ylim = c(0, 1), col = NULL,
    lty = NULL, main = NULL, xlab = "Latent trait",
    ylab = "Probability", add = FALSE, \dots)
}

\arguments{
  \item{object}{a fitted model object of class \code{"raschmodel"},
    \code{"rsmodel"}, \code{"pcmodel"}, \code{"plmodel"} or \code{"gpcmodel"}.}
  \item{ref}{argument passed over to internal calls of \code{\link{predict}}.}
  \item{items}{character or numeric, specifying the items for which response
    curves should be visualized.}
  \item{names}{character, specifying labels for the items.}
  \item{layout}{matrix, specifying how the response curve plots of different
    items should be arranged.}
  \item{xlim, ylim}{numeric, specifying the x and y axis limits.}
  \item{col}{character, specifying the colors of the response curve lines. The
    length of \code{col} should be the maximum number of available categories.}
  \item{lty}{numeric, specifying the line type of the response curve lines. The
    length of \code{lty} should either be one or the maximum number of
    available categories. In the first case, a single line type is used for all
    category response curves. In the latter case, separate line types for each
    category response curve are used.}
  \item{main}{character, specifying the overall title of the plot.}
  \item{xlab, ylab}{character, specifying the x and y axis labels.}
  \item{add}{logical. If \code{TRUE}, new response curves are added to an
    existing plot. Only possible when a single item is visualized.}
  \item{\dots}{further arguments passed to internal calls of
    \code{\link{matplot}}.}
}

\details{
  The response curve plot visualization illustrates the predicted probabilities
  as a function of the ability parameter \eqn{\theta} under a certain IRT model.
  This type of visualization is sometimes also called item/category operating
  curves or item/category characteristic curves.
}

\seealso{\code{\link{regionplot}}, \code{\link{profileplot}},
  \code{\link{infoplot}}, \code{\link{piplot}}}

\examples{
## load verbal aggression data
data("VerbalAggression", package = "psychotools")

## fit Rasch, rating scale and partial credit model to verbal aggression data
rmmod <- raschmodel(VerbalAggression$resp2)
rsmod <- rsmodel(VerbalAggression$resp)
pcmod <- pcmodel(VerbalAggression$resp)

## curve plots of the dichotomous RM
plot(rmmod, type = "curves")

## curve plots under the RSM for the first six items of the data set
plot(rsmod, type = "curves", items = 1:6)

## curve plots under the PCM for the first six items of the data set with
## custom labels
plot(pcmod, type = "curves", items = 1:6, names = paste("Item", 1:6))

## compare the predicted probabilities under the RSM and the PCM for a single
## item
plot(rsmod, type = "curves", item = 1)
plot(pcmod, type = "curves", item = 1, lty = 2, add = TRUE)
legend(x = "topleft", y = 1.0, legend = c("RSM", "PCM"), lty = 1:2, bty = "n")

if(requireNamespace("mirt")) {
## fit 2PL and generaliced partial credit model to verbal aggression data
twoplmod <- plmodel(VerbalAggression$resp2)
gpcmod <- gpcmodel(VerbalAggression$resp)

## curve plots of the dichotomous 2PL
plot(twoplmod, type = "curves", xlim = c(-6, 6))

## curve plots under the GPCM for the first six items of the data set
plot(gpcmod, type = "curves", items = 1:6, xlim = c(-6, 6))
}
}

\keyword{aplot}

