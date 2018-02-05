\name{mutatedCopies}
\alias{mutatedCopies}
\title{
Predicts the number of mutated copies
}
\description{
Computes the maximum likelihood estimate of the number of mutated copies given the
observed VAF
}
\usage{
sq_hat = mutatedCopies(f_hat, n, qt, q2, alpha, e, d)
}
\arguments{
\item{f_hat}{numeric vector of observed VAF where f_hat = 0 \eqn{\ldots} 1}
\item{n}{integer vector of total number of reads covering mutated loci}
\item{qt}{integer vector of total copy numbers}
\item{q2}{integer vector of major allele counts}
\item{alpha}{numeric for the global cellular fraction where alpha = 0 \eqn{\ldots} 1}
\item{e}{numeric for expected sequencing error where e = 0 \eqn{\ldots} 1}
\item{d}{character string for the underlying distribution 'beta' or 'binomial'}
}
\value{
Integer vector of estimates of mutated copies
}
\author{
\cr
David Brown\cr
\email{brownd7@mskcc.org}\cr
}
\keyword{}