\name{SoundQuality}
\alias{SoundQuality}

\title{Quality of Multichannel Reproduced Sound}

\description{
  Paired comparison judgments of 40 selected listeners with respect to
  eight audio reproduction modes and four types of music. 
}

\usage{data("SoundQuality")}

\format{
A data frame containing 783 observations on 6 variables.
\describe{
  \item{id}{Factor. Listener ID.}
  \item{time}{Factor. Listening experiment before or after elicitation and
    scaling of more specific auditory attributes.}
  \item{progmat}{Factor. The program material: Beethoven, Rachmaninov, Steely
    Dan, Sting.}
  \item{repet}{The repetition within each time point.}
  \item{session}{The experimental session coding the presentation order of
    the program material.}
  \item{preference}{Paired comparison of class \code{\link{paircomp}}.
    Preferences for all 28 paired comparisons from 8 audio reproduction modes:
    Mono, Phantom Mono, Stereo, Wide-Angle Stereo, 4-channel Matrix,
    5-channel Upmix 1, 5-channel Upmix 2, and 5-channel Original.}
}
}

\details{
  The data were collected within a series of experiments conducted at the
  Sound Quality Research Unit (SQRU), Department of Acoustics, Aalborg
  University, Denmark, between September 2004 and March 2005.

  The results of scaling listener preference and spatial and timbral auditory
  attributes are reported in Choisel and Wickelmaier (2007).

  Details about the loudspeaker setup and calibration are given in Choisel and
  Wickelmaier (2006).

  The attribute elicitation procedure is described in Wickelmaier and
  Ellermeier (2007) and in Choisel and Wickelmaier (2006).

  The selection of listeners for the experiments is described in Wickelmaier
  and Choisel (2005).

  An extended version of this data set, including judgments on spatial and
  timbral auditory attributes and including listener variables, is available
  via \code{data("soundquality", package = "eba")}.
}

\references{
  Choisel S, Wickelmaier F (2006).
    Extraction of Auditory Features and Elicitation of Attributes for the Assessment of Multichannel Reproduced Sound.
    \emph{Journal of the Audio Engineering Society}, \bold{54}(9), 815--826.

  Choisel S, Wickelmaier F (2007).
    Evaluation of Multichannel Reproduced Sound: Scaling Auditory Attributes Underlying Listener Preference.
    \emph{Journal of the Acoustical Society of America}, \bold{121}(1), 388--400. 
    \doi{10.1121/1.2385043}

  Wickelmaier F, Choisel S (2005).
    Selecting Participants for Listening Tests of Multichannel Reproduced Sound.
    Presented at the AES 118th Convention, May 28--31, Barcelona, Spain, convention paper 6483.

  Wickelmaier F, Ellermeier W (2007).
    Deriving Auditory Features from Triadic Comparisons.
    \emph{Perception & Psychophysics}, \bold{69}(2), 287--297.
    \doi{10.3758/BF03193750}
}

\seealso{\code{\link{paircomp}}}

\examples{
data("SoundQuality", package = "psychotools")
summary(SoundQuality$preference)
ftable(xtabs(~ time + repet + progmat, data = SoundQuality))
}

\keyword{datasets}
