'mutatedCopies' <- function (f_hat, n, qt, q2, alpha, e, d = "beta")
{
	sq_hat = vector(mode="numeric", length=length(f_hat))
	if (d=="beta") {
		for (i in 1:length(f_hat)) {
			sq_hat[i] = .backend.beta(f_hat[i], n[i], qt[i], q2[i], alpha, e)
		}
	} else if (d=="binomial") {
		for (i in 1:length(f_hat)) {
			sq_hat[i] = .backend.binomial(f_hat[i], n[i], qt[i], q2[i], alpha, e)
		}
	}
	return(invisible(sq_hat))
}

