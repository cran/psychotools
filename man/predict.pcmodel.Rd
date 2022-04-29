\name{predict.pcmodel}
\alias{predict.pcmodel}

\alias{predict.rsmodel}
\alias{predict.raschmodel}
\alias{predict.gpcmodel}
\alias{predict.nplmodel}

\title{Predict Methods for Item Response Models}

\description{
  Prediction of (cumulated) response probabilities and responses based on fitted
  item response models.
}

\usage{
\S3method{predict}{pcmodel}(object, newdata = NULL, type = c("probability",
  "cumprobability", "mode", "median", "mean", "category-information",
  "item-information", "test-information"), ref = NULL, \dots)
}

\arguments{
  \item{object}{a fitted model object whose item parameters should be used for
    prediction.}
  \item{newdata}{an optional (possibly named) vector of person parameters
    used for prediction. If \code{NULL} (the default), the person parameters of
    the subjects used to fit the model in \code{object} are used.}
  \item{type}{character of length one which determines the type of
    prediction (see details below).}
  \item{ref}{arguments passed over to internal calls of \code{itempar} or
    \code{threshpar}. Not used for models estimated via MML.}
  \item{\dots}{further arguments which are currently not used.}
}

\details{
  Depending on the value of \code{type} either probabilities, responses or
  some form of information under the model specified in \code{object} are
  returned:

  If \code{type} is \code{"probability"}, the category response probabilities
  are returned.
  
  If \code{type} is \code{"cumprobability"}, the cumulated category response
  probabilities are returned, i.e., \eqn{P(X_{ij} \geq k)} with \eqn{k}
  corresponding to the categories of item \eqn{j}.

  If \code{type} is \code{"mode"}, the most probable category response for a
  given subject and item is returned.

  If \code{type} is \code{"median"}, the first category \eqn{k} where
  \eqn{P(X_{ij} = k) \geq 0.5} is returned.

  If \code{type} is \code{"mean"}, the rounded expected category response,
  i.e., \eqn{E(X_{ij}|\theta_{i})}, is returned.

  If \code{type} is \code{"category-information"}, the item-category
  information as suggested by Bock (1972) is returned.
  
  If \code{type} is \code{"item-information"}, the item information as
  suggested by Samejima (1974) is returned.
  
  If \code{type} is \code{"test-information"}, the sum over the individual
  item information values is returned.
}
 
\value{
  A (possibly named) numeric matrix with rows corresponding to subjects and
  columns corresponding to the whole test, the single items or categories. The
  exact content depends on the value of \code{type} (see details above).
}

\references{
  Bock RD (1972).
    Estimating Item Parameters and Latent Ability When Responses Are Scored in
    Two or More Nominal Categories.
    \emph{Psychometrika}, \bold{37}(1), 29--51.

  Samejima F (1974).
    Normal Ogive Model on the Continuous Response Level in the Multidimensional
    Latent Space.
    \emph{Psychometrika}, \bold{39}(1), 111--121.
}

\seealso{
  The help page of the generic function \code{\link{predict}} and other
  predict methods (e.g., \code{\link{predict.lm}}, \code{\link{predict.glm}},
  \dots)
}


\examples{
o <- options(digits = 4)

## load verbal aggression data
data("VerbalAggression", package = "psychotools")

## fit a partial credit model to first ten items
pcmod <- pcmodel(VerbalAggression$resp[, 1:10])

## predicted response probabilities for each subject and category (the default)
head(predict(pcmod), 3)

## predicted mode (most probable category) for certain subjects whose person
## parameters are given via argument "newdata"
predict(pcmod, type = "mode",
  newdata = c("Sarah" = 1.2, "Michael" = 0.1, "Arnd" = -0.8))

## rounded expected category value for the same subjects
predict(pcmod, type = "mean",
  newdata = c("Sarah" = 1.2, "Michael" = 0.1, "Arnd" = -0.8))

## in the Rasch model mode, mean and median are the same
raschmod <- raschmodel(VerbalAggression$resp2[, 1:10])
med <- predict(raschmod, type = "median")
mn <- predict(raschmod, type = "mean")
mod <- predict(raschmod, type = "mode")

head(med, 3)

all.equal(med, mn)
all.equal(mod, mn)

options(digits = o$digits)
}

\keyword{regression}

