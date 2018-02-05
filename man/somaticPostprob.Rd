\name{somaticPostprob}
\alias{somaticPostprob}
\title{
Posterior probability of clonal SNV
}
\description{
Computes posterior probabilities for SNV given observed VAF, total
copy numbers, allelic information and cellular fraction
}
\usage{
pr = somaticPostprob(f_hat, n, qt, q2, alpha, w, sthetaq, e)
}
\arguments{
\item{f_hat}{numeric vector of observed VAFs where f_hat = 0 \eqn{\ldots} 1}
\item{n}{integer vector of total number of reads covering mutated loci}
\item{qt}{integer vector of total copy numbers}
\item{q2}{integer vector of major allele counts}
\item{alpha}{numeric estimate of cancer cell fraction}
\item{w}{list of numeric mixture weights}
\item{sthetaq}{numeric vector}
\item{e}{numeric error rate where e = 0 \eqn{\ldots} 1}
}
\value{
list of two as follows:
\item{ppr}{data frame of posterior probabilities with following column description:\cr
		   \enumerate{
    			\item{Pr_somatic_clonal:	}{probability that SNV is somatic and strictly clonal}
    			\item{Pr_germline:	}{probability that SNV is a germline mutation}
    			\item{Pr_subclonal:	}{probability that SNV is somatic and subclonal}
    			\item{Pr_subclonal_wt0:	}{probability that SNV is subclonal when qt = 1, 0 otherwise}
    			\item{Pr_wt0:	}{probability that q2 copies are mutated if locus has uniparental allele (qt = q2), 0 otherwise}
    			\item{Pr_ge2:	}{warning: undocumented \& not in use}
    			\item{Pr_GL_som_HZ_alt:	}{warning: undocumented \& not in use}
    			\item{Pr_GL_som_HZ_ref:	}{warning: undocumented \& not in use}
    			\item{Pr_cryptic_SCNA:	}{warning: undocumented}
    			\item{LL:	}{log likelihood}
    			\item{sq:	}{most probable number of mutated copies where sq = 1 \eqn{ldots} q2}
    			\item{Pr_somatic:	}{probability that SNV is somatic i.e. Pr_somatic_clonal + Pr_subclonal}
    		}
}
\item{som_mut_Q_tab}{list of numeric mixture weights}
}
\note{
This functionality has been lifted straight from ABSOLUTE with the unfortunate caveat that
most of it is undocumented
}
\author{
\cr
David Brown\cr
\email{brownd7@mskcc.org}\cr
}
\keyword{}