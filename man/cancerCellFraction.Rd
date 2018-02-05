\name{cancercellFraction}
\alias{cancercellFraction}
\title{
Cancer cell fraction
}
\description{
Computes the maximum likelihood estimate of cancer cell fractions given observed VAF
}
\usage{
alpha_hat = cancercellFraction(f_hat, n, qt, q2, alpha, e=0.01)
}
\arguments{
\item{f_hat}{numeric vector of variant allele fractions where f_hat = 0 \eqn{\ldots} 1}
\item{n}{integer vector of total depth of coverage}
\item{qt}{integer vector of total copy numbers}
\item{q2}{integer vector of major allele}
\item{alpha}{numeric value of global purity where alpha = 0 \eqn{\ldots} 1}
\item{e}{numeric value of sequencing error rate where e = 0 \eqn{\ldots} 1}
}
\value{
Data frame with columns as follows:
	\enumerate{
		\item{cell_multiplicity:	}{}
		\item{scaled_cell_frac:	}{}
		\item{cancer_cell_frac:	}{}
		\item{ccf_95CI_low:	}{}
		\item{ccf_95CI_high:	}{}
		\item{Pr_somatic_clonal:	}{probability that SNV is somatic and clonal}
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
\author{
\cr
David Brown\cr
\email{brownd7@mskcc.org}\cr
}
\keyword{}