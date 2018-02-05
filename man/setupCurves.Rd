\name{setupCurves}
\alias{setupCurves}
\title{
Plot dilution curves at given multiplicities
}
\description{
Plots dilution curves of expected VAFs at all multiplicites i.e. 2/1, 2/2, 3/1, 3/2, etc
given contamination with normal cells harboring two non mutated copies
}
\usage{
setupCurves(nqt, col, lty, f_hat, qt, sq)
}
\arguments{
\item{nqt}{integer maximum number of total copies}
\item{col}{character vector of colors for the different copy numbers where \eqn{||col|| = nqt}}
\item{lty}{numeric vector of line types for the different multiplicities where \eqn{||lty|| = nqt}}
\item{f_hat}{optional numeric vector of observed VAFs where f_hat = 0 \eqn{\ldots} 1}
\item{qt}{optional integer vector total copy numbers where \eqn{||qt|| = ||f_hat||}}
\item{sq}{optional integer vector of number of mutated copies where \eqn{||sq|| = ||f_hat||}}
}
\value{
None
}
\note{
If f_hat, qt and sq are provided, the function \code{cellularFraction(...)} is called
to compute the maximum likelihood of cellular fractions. These are displayed as a
scatter plot above the dilution curves. For clarity, user should not attempt to go beyond
\eqn{nqt = 3} and similarly remove points in f_hat, qt and sq where \eqn{qt>3}.
}
\author{
\cr
David Brown\cr
\email{brownd7@mskcc.org}\cr
}
\keyword{}