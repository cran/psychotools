\name{regionplot}
\alias{regionplot}

\title{Region Plots for IRT Models}

\description{
  Base graphics plotting function for region plot visualization of IRT models.
}

\usage{
  regionplot(object, parg = list(type = NULL, ref = NULL, alias = TRUE),
    names = TRUE, main = NULL, xlab = "", ylab = "Latent trait", ylim = NULL,
    off = 0.1, col = NULL, linecol = 2, srt = 45, adj = c(1.1, 1.1),
    axes = TRUE, \dots)
}

\arguments{
  \item{object}{a fitted model object of class \code{"raschmodel"},
    \code{"rsmodel"}, \code{"pcmodel"}, \code{"nplmodel"} or \code{"gpcmodel"}.}
  \item{parg}{list of arguments passed over to internal calls of
    \code{\link{threshpar}}. See the help page of \code{\link{threshpar}}
    for more details.}
  \item{names}{logical or character. If \code{TRUE}, the names of the items are
    displayed on the x-axis. If \code{FALSE}, numbers of items are shown.
    Alternatively a character vector of the same length as the number of items
    can be supplied.}
  \item{main}{character, specifying the overall title of the plot.}
  \item{xlab, ylab}{character, specifying the x and y axis labels.}
  \item{ylim}{numeric, specifying the y axis limits.}
  \item{off}{numeric, the distance (in scale units) between two item
    rectangles.}
  \item{col}{character, list or function, specifying the colors of the regions.
    Either a single vector with \eqn{k} color names, a list with \eqn{m}
    elements and each element is a character vector with color names for the
    regions of item \eqn{j} or a color-generating function like, e.g.,
    \code{gray.colors}, which is then directly used to create the color names.}
  \item{linecol}{color for lines indicating \dQuote{hidden} categories.}
  \item{srt, adj}{numeric. Angle (\code{srt}) and adjustment (\code{adj}) in
    case names (rather than numbers) are used as x-axis labels. These are
    passed to \code{\link[graphics]{text}}.}
  \item{axes}{logical. Should axes be drawn?}
  \item{\dots}{further arguments passed to \code{\link{plot}}.}
}

\details{
  The region plot visualization implemented here was already used by Van der
  Linden and Hambleton (1997) in the context of IRT and has been called "effect
  plots" by Fox & Hong (2009). In our implementation, these plots show,
  dependent on the chosen type of threshold parameters, different regions for
  the categories of an item over the theta axis. If \code{type} is set to
  \code{"modus"}, the cutpoints correspond to the threshold parameters and the
  rectangles mark the theta regions where a category is the single most probable
  category chosen with a certain value of the latent trait. If \code{type} is
  set to \code{"median"}, the cutpoints correspond to the point on the theta
  axis, where the cumulative probability to score in category \eqn{k} or higher
  is 0.5, i.e., \eqn{P(X_{ij} \geq k) = 0.5}. If set to \code{"mean"}, the
  cutpoints correspond to the point on the theta axis where the expected score
  \eqn{E(X_{ij})} is exactly between two categories, e.g., 0.5 for a dichotomous
  item.

  If \code{type} is set to \code{"mode"} and there are unordered threshold
  parameters, the location of the original threshold parameters are indicated by
  red dashed lines.
}

\references{
  Fox J, Hong J (2009).
    Effect Displays in R for Multinomial and Proportional-Odds Logit Models:
    Extensions to the effects Package.
    \emph{Journal of Statistical Software}, \bold{32}(1), 1--24.

  Van der Linden WJ, Hambleton RK (1997).
    \emph{Handbook of Modern Item Response Theory}.
    Springer, New York.
}

\seealso{\code{\link{curveplot}}, \code{\link{profileplot}},
  \code{\link{infoplot}}, \code{\link{piplot}}}

\examples{
## load verbal aggression data
data("VerbalAggression", package = "psychotools")

## fit a Partial credit model to the items of the first other-to-blame
## situation: "A bus fails to stop for me"
pcm <- pcmodel(VerbalAggression$resp[, 1:6])

## a region plot with modus as cutpoint and custom labels
lab <- paste(rep(c("Curse", "Scold", "Shout"), each = 2),
  rep(c("Want", "Do"), 3 ), sep = "-")
plot(pcm, type = "regions", names = lab)

## compare the cutpoints (with ylim specified manually)
opar <- par(no.readonly = TRUE)
ylim <- c(-2, 2)
layout(matrix(1:3, ncol = 1))
plot(pcm, type = "regions", parg = list(type = "mode"),
  main = "Modus as Cutpoint", ylim = ylim) 
plot(pcm, type = "regions", parg = list(type = "median"),
  main = "Median as Cutpoint", ylim = ylim)
plot(pcm, type = "regions", parg = list(type = "mean"),
  main = "Mean as Cutpoint", ylim = ylim)
par(opar)

## PCM for full verbal aggression data set
pcm_va <- pcmodel(VerbalAggression$resp)
plot(pcm_va, type = "regions")

if(requireNamespace("mirt")) {
## generalized partial credit model for full verbal aggression data set
gpcm_va <- gpcmodel(VerbalAggression$resp)
plot(gpcm_va, type = "regions")
}
}

\keyword{aplot}

