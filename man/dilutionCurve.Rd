\name{dilutionCurve}
\alias{dilutionCurve}
\title{
Compute dilution curve
}
\description{
Computes the curve of expected VAF given total copy number, number of mutated copies and
dilution with normal cells harboring two non mutated copies
}
\usage{
fsq = dilutionCurve(sq, qt, alpha)
}
\arguments{
\item{sq}{integer number of mutated copies}
\item{qt}{integer total copy number}
\item{alpha}{numeric vector of cellular fractions over which fsq will be computed}
}
\value{
Numeric vector of expected VAF where fsq = 0 \eqn{\ldots} 1
}
\author{
\cr
David Brown\cr
\email{brownd7@mskcc.org}\cr
}
\keyword{}