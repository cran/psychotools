\name{infoplot}
\alias{infoplot}

\title{Information Plots for IRT Models}

\description{
  Base graphics plotting function for information plot visualization of IRT models.
}

\usage{
  infoplot(object, what = c("categories", "items", "test"),
    ref = NULL, items = NULL, names = NULL, layout = NULL, xlim = NULL,
    ylim = NULL, col = NULL, lty = NULL, lwd = NULL, main = NULL, legend = TRUE,
    xlab = "Latent trait", ylab = "Information", add = FALSE, \dots)
}

\arguments{
  \item{object}{a fitted model object of class \code{"raschmodel"},
    \code{"rsmodel"}, \code{"pcmodel"}, \code{"plmodel"} or \code{"gpcmodel"}.}
  \item{what}{character, specifying the type of information to visualize.}
  \item{ref}{argument passed over to internal calls of \code{\link{predict}}.}
  \item{items}{character or numeric, specifying the items for which information
    curves should be visualized.}
  \item{names}{character, specifying labels for the items.}
  \item{layout}{matrix, specifying how the item or category information curves
    of different items should be arranged. If \code{null} and \code{what} is set
    to \code{"items"}, the item information curves are overlayed within a single
    plot.}
  \item{xlim, ylim}{numeric, specifying the x and y axis limits.}
  \item{col}{character, specifying the colors of the test, item or category
    information curves.}
  \item{lty}{numeric, specifying the line type of the information curves.}
  \item{lwd}{numeric, specifying the line width of the information curves.}
  \item{main}{character, specifying the overall title of the plot.}
  \item{legend}{logical, specifying if a legend is drawn when multiple item
    information curves are overlayed. The labels in the legend correspond to
    the item names (which can be specified in the argument \code{names}).}
  \item{xlab, ylab}{character, specifying the x and y axis labels.}
  \item{add}{logical. If \code{TRUE}, new information curves are added to
    an existing plot. Only possible for a test or a single item information
    curve.}
  \item{\dots}{further arguments passed to internal calls of
    \code{\link{matplot}}.}
}

\details{
  The information plot visualization illustrates the test, item or category
  information as a function of the ability parameter \eqn{\theta} under a
  certain IRT model. Further details on the computation of the displayed
  information can be found on the help page of the function
  \code{\link{predict.pcmodel}}.
}

\seealso{\code{\link{curveplot}}, \code{\link{regionplot}},
  \code{\link{profileplot}}, \code{\link{piplot}}}

\examples{
## load verbal aggression data
data("VerbalAggression", package = "psychotools")

## fit Rasch and partial credit model to verbal aggression data
rmmod <- raschmodel(VerbalAggression$resp2)
pcmod <- pcmodel(VerbalAggression$resp)

## category information plots for all items under the dichotomous RM
plot(rmmod, type = "information", what = "categories")

## category information plots for all items under the PCM
plot(pcmod, type = "information", what = "categories")

## overlayed item information plots for the first six items of the
## data set under the PCM
plot(pcmod, type = "information", what = "items", items = 1:6)

## a comparison of the item information for the first six items under the
## dichotomous RM and the PCM
plot(pcmod, type = "information", what = "items", items = 1:6,
  xlim = c(-5, 5))
plot(rmmod, type = "information", what = "items", items = 1:6,
  lty = 2, add = TRUE)
legend(x = "topright", legend = c("PCM", "RM"), lty = 1:2, bty = "n")

## a comparison of the test information based on all items of the
## data set under the dichotomous RM and the PCM
plot(pcmod, type = "information", what = "test", items = 1:6, xlim = c(-5, 5))
plot(rmmod, type = "information", what = "test", items = 1:6, lty = 2,
  add = TRUE)
legend(x = "topright", legend = c("PCM", "RM"), lty = 1:2, bty = "n")

if(requireNamespace("mirt")) {
## fit 2PL to verbal aggression data
twoplmod <- plmodel(VerbalAggression$resp2)

## category information plots for all items under the dichotomous 2PL
plot(twoplmod, type = "information", what = "categories")
}
}

\keyword{aplot}

