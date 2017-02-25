\name{dilutionCurve}
\alias{dilutionCurve}
\title{
Compute dilution curve
}
\description{
Computes the curve of expected VAFs given total copy number, number of mutated copies and
dilution with normal cells harboring two non mutated copies
}
\usage{
fsq = dilutionCurve(sq, qt)
fsq = dilutionCurve(sq, qt, alpha)
}
\arguments{
\item{sq}{integer number of mutated copies}
\item{qt}{integer total copy number}
\item{alpha}{numeric vector of cancer cell fractions over which fsq will be computed}
}
\value{
numeric vector of expected VAFs where fsq = 0 \eqn{\ldots} 100
}
\author{
\cr
David Brown\cr
\email{dbrown@ulb.ac.be}\cr
\url{http://homepages.ulb.ac.be/~dbrown/}
}
\seealso{
\code{\link{setupCurves}}
}
\keyword{}