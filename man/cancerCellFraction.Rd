\name{cancerCellFraction}
\alias{cancerCellFraction}
\title{
Cancer cell fraction
}
\description{
Computes the maximum likelihood estimate of cancer cell fractions given observed VAFs
}
\usage{
alpha_hat = cancerCellFraction(f_hat, qt, sq)
}
\arguments{
\item{f_hat}{numeric vector of variant allele fractions where f_hat = 0 \eqn{\ldots} 1}
\item{qt}{numeric vector of integer estimates of total copy numbers}
\item{sq}{numeric vector of integer estimates of mutated copy numbers}
}
\value{
numeric vector of cancer cell fractions where alpha_hat = 0 \eqn{\ldots} 1
}
\author{
\cr
David Brown\cr
\email{dbrown@ulb.ac.be}\cr
\url{http://homepages.ulb.ac.be/~dbrown/}
}
\seealso{
\code{\link{mutatedCopies}}
}
\keyword{}