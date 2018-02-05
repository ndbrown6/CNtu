\name{cellularFraction}
\alias{cellularFraction}
\title{
Cellular fraction
}
\description{
Computes the maximum likelihood estimate of cancer cell fractions given observed VAF
}
\usage{
alpha_hat = cellularFraction(f_hat, qt, sq)
}
\arguments{
\item{f_hat}{numeric vector of variant allele fractions where f_hat = 0 \eqn{\ldots} 1}
\item{qt}{numeric vector of integer estimates of total copy numbers}
\item{sq}{numeric vector of integer estimates of mutated copy numbers}
}
\value{
Numeric vector of cancer cell fractions where alpha_hat = 0 \eqn{\ldots} 1
}
\author{
\cr
David Brown\cr
\email{brownd7@mskcc.org}\cr
}
\keyword{}