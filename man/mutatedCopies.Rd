\name{mutatedCopies}
\alias{mutatedCopies}
\title{
Predict number of mutated copies
}
\description{
Computes the maximum likelihood estimate of the number of mutated copies given the
observed VAFs
}
\usage{
sq_hat = mutatedCopies(f_hat, n, qt, q2, alpha, e)
sq_hat = mutatedCopies(f_hat, n, qt, q2, alpha, e, dist)
}
\arguments{
\item{f_hat}{numeric vector of observed VAF where f_hat = 0 \eqn{\ldots} 1}
\item{n}{integer vector of total number of reads covering mutated loci}
\item{qt}{integer vector of total copy numbers}
\item{q2}{integer vector of major allele counts}
\item{alpha}{numeric for the estimate of cancer cell fraction}
\item{e}{numeric for expected sequencing error}
\item{dist}{character string for the underlying distribution 'beta' or 'binomial'}
}
\value{
integer vector of estimates of mutated copies
}
\author{
\cr
David Brown\cr
\email{dbrown@ulb.ac.be}\cr
\url{http://homepages.ulb.ac.be/~dbrown/}
}
\seealso{
}
\keyword{}