\name{profileplot}
\alias{profileplot}

\title{Profile Plots for IRT Models}

\description{
  Base graphics plotting function for profile plot visualization of IRT models.
}

\usage{
  profileplot(object,
    what = c("items", "thresholds", "discriminations", "guessings", "uppers"),
    parg = list(type = NULL, ref = NULL, alias = TRUE, logit = FALSE), index = TRUE,
    names = TRUE, main = NULL, abbreviate = FALSE, ref = TRUE,
    col = "lightgray", border = "black", pch = NULL, cex = 1,
    refcol = "lightgray", linecol = "black", lty = 2, ylim = NULL,
    xlab = NULL, ylab = NULL, add = FALSE, srt = 45, adj = c(1.1, 1.1),
    axes = TRUE, \dots)
}

\arguments{
  \item{object}{a fitted model object of class \code{"raschmodel"},
    \code{"rsmodel"}, \code{"pcmodel"}, \code{"nplmodel"} or \code{"gpcmodel"}.}
  \item{what}{character, specifying the type of parameters to be plotted.}
  \item{parg}{list of arguments passed over to internal calls of
    \code{\link{itempar}}, \code{\link{threshpar}}, \code{\link{discrpar}},
    \code{\link{guesspar}}, or \code{\link{upperpar}}.}
  \item{index}{logical, should different indexes for different items be used?}
  \item{names}{logical or character. If \code{TRUE}, the names of the items are
    displayed on the x-axis. If \code{FALSE}, numbers of items are shown.
    Alternatively a character vector of the same length as the number of items
    can be supplied.}
  \item{main}{character, specifying the overall title of the plot.}
  \item{abbreviate}{logical or numeric, specifying whether object names are to
    be abbreviated. If numeric this controls the length of the abbreviation.}
  \item{ref}{logical, whether to draw a horizontal line for the reference level.
    Only takes effect if argument \code{what} is \code{"items"} or
    \code{"discriminations".}}
  \item{col, border, pch, cex}{graphical appearance of plotting symbols. Can be
    of the same length as the number of items, i.e., a different graphical
    appearance is used for each item. If \code{what = "thresholds"}, \code{col}
    and \code{pch} can be matrices with a number of columns equal to the number
    of threshold parameters per item resulting in different symbols and colors
    used for different threshold parameter profiles.}
  \item{refcol}{character, specifying the line color for the reference line
    (if \code{ref} is set to \code{TRUE}).}
  \item{linecol}{character or numeric, specifying the line color to be used for
    the profiles.}
  \item{lty}{numeric, specifying the line type for the profiles.}
  \item{ylim}{numeric, specifying the y axis limits.}
  \item{xlab, ylab}{character, specifying the x and y axis labels.}
  \item{add}{logical. If \code{TRUE}, new plotted profiles are added to an
    existing plot.}
  \item{srt, adj}{numeric. Angle (\code{srt}) and adjustment (\code{adj}) in
    case names (rather than numbers) are used as x-axis labels. These are
    passed to \code{\link[graphics]{text}}.}
  \item{axes}{logical. Should axes be drawn?}
  \item{\dots}{further arguments passed over to \code{\link{plot}}.}
}

\details{
  The profile plot visualization illustrates profiles of specific estimated
  parameters under a certain IRT model.
}


\seealso{\code{\link{curveplot}}, \code{\link{regionplot}},
  \code{\link{infoplot}}, \code{\link{piplot}}}

\examples{
## load verbal aggression data
data("VerbalAggression", package = "psychotools")

## fit Rasch, rating scale and partial credit model to verbal aggression data
rmmod <- raschmodel(VerbalAggression$resp2)
rsmod <- rsmodel(VerbalAggression$resp)
pcmod <- pcmodel(VerbalAggression$resp)

## profile plots of the item parameters of the three fitted IRT models
plot(rmmod, type = "profile", what = "items", col = 4)
plot(rsmod, type = "profile", what = "items", col = 2, add = TRUE)
plot(pcmod, type = "profile", what = "items", col = 3, add = TRUE)
legend(x = "topleft", legend = c("RM", "RSM", "PCM"), col = 1,
  bg = c(4, 2, 3), pch = 21, bty = "n")

## profile plots of the threshold parameters of type "mode"
plot(rmmod, type = "profile", what = "thresholds", parg = list(type = "mode"))
plot(rsmod, type = "profile", what = "thresholds", parg = list(type = "mode"))
plot(pcmod, type = "profile", what = "thresholds", parg = list(type = "mode"))

## profile plot of the discrimination parameters of the dichotomous RM
plot(rmmod, type = "profile", what = "discrimination")

\donttest{
if(requireNamespace("mirt")) {
## fit 2PL and generalized partial credit model to verbal aggression data
twoplmod <- nplmodel(VerbalAggression$resp2)
gpcmod <- gpcmodel(VerbalAggression$resp)

## profile plot of the discrimination parameters of a dichotomous 2PL
plot(twoplmod, type = "profile", what = "discrimination")

## profile plot of the item parameters of the 2PL and GPCM
plot(twoplmod, type = "profile", what = "items", col = 4)
plot(gpcmod, type = "profile", what = "items", col = 2, add = TRUE)
}
}
}

\keyword{aplot}

