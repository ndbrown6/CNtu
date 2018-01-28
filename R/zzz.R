'.backend.beta' <- function (f_hat, n, qt, q2, alpha, e)
{
	sq = 1:q2
	fsq = sq*alpha/(2*(1-alpha)+alpha*qt)
	fe = (1/(3-2*e))*((3-4*e)*fsq + e)
	p = dbeta(x=fe, round(n*f_hat)+1, round(n*(1-f_hat))+1, log=TRUE)
	sq_hat = sq[which.max(p)]
	return(invisible(sq_hat))
}

'.backend.binomial' <- function (f_hat, n, qt, q2, alpha, e)
{
	sq = 1:q2
	fsq = sq*alpha/(2*(1-alpha)+alpha*qt)
	fe = (1/(3-2*e))*((3-4*e)*fsq + e)
	p = dbinom(x=round(n*f_hat), n, fe, log=TRUE)
	sq_hat = sq[which.max(p)]
	return(invisible(sq_hat))
}

'.backend.addPoints' <- function (f_hat, qt, sq, col)
{
	cols = col[qt]
	ccf = cancerCellFraction(f_hat, qt, sq)
	points(ccf*100, f_hat*100, type="p", pch="X", col=cols, cex=1.2, lwd=1.5)
}

'.backend.GetHscnSomaticMutComb' <- function(alpha, qt, q2)
{
	sq = 1:q2
	fsq = alpha*sq/(2*(1-alpha) + alpha*qt)
	fsq = cbind(sq, fsq)
	return(invisible(fsq))
}

'.backend.SomaticMutPrior' <- function(W, N, thetaq)
{
	w = thetaq[c(1:N)]
	w = w/sum(w)
	w = W[["SM"]]*w
	return(invisible(w))
}

'.backend.FhatCombPost' <- function(f_hat, n, fsq, e=0.01)
{
	fe = (1/(3-2*e))*((3-4*e)*fsq + e)
	logpp = matrix(NA, nrow=length(f_hat), ncol=length(fe))
	for (i in 1:length(fe)) {
		logpp[,i] = dbeta(fe[i], f_hat*n + 1, (1-f_hat)*n + 1, log=TRUE)
	}
	return(invisible(logpp))
}

'.backend.GermlineMutPrior' <- function(W)
{
	w = c(49, 49, 1)
	w = w/sum(w)
	w = W[["GL"]]*w
	return(invisible(w))
}

'.backend.GetHscnGermlineMutComb' <- function(alpha, qt, q2)
{
	eps = 0.001
	sq = c(qt-q2, q2, qt)
	fsq = ((1-alpha) + alpha*sq)/(2*(1-alpha) + alpha*qt)
	fsq[3] = 1 - eps
	fsq = cbind(sq, fsq)
	return(invisible(fsq))
}

'.backend.UnifSubclonalPost' <- function(f_hat, sdelta, n, e=0.01)
{
	fe = (1/(3-2*e))*((3-4*e)*sdelta + e)
	beta.init = pbeta(fe, f_hat*n + 1, (1-f_hat)*n + 1, lower.tail=TRUE, log.p=TRUE)
	loglik = beta.init + log(1/fe)
	return(invisible(loglik))
}

'.backend.CrypticScnaPost' <- function (f_hat, sdelta, n, e=0.01)
{
	fe = (1/(3-2*e))*((3-4*e)*sdelta + e)
	beta.init = pbeta(fe, f_hat*n + 1, (1 - f_hat) * n + 1, lower.tail=FALSE, log.p=TRUE)
	loglik = beta.init + log(1/fe)
	return(loglik)
}

'.backend.LogAdd' <- function(X)
{
    if (is.vector(X)) {
        mix = which.max(X)
        max = X[mix]
        res = max + log(sum(exp(X - max)))
    } else if (is.matrix(X)) {
        mv = apply(X, 1, max)
        res = mv + log(rowSums(exp(X - mv)))
    }
    return(invisible(res))
}

'.backend.jitter' <- function (x, y)
{
	qx = seq(from=min(x, na.rm=TRUE), to=max(x, na.rm=TRUE), length=10)
	nx = vector(mode="numeric", length=(length(qx)-1))
	for (i in 1:length(nx)) {
		nx[i] = sum(x>=qx[i] & x<qx[i+1], na.rm=TRUE)
	}
	nx = .25*(nx - min(nx))/(max(nx)-min(nx))
	z = rep(y, length(x))
	for (i in 1:length(nx)) {
		index = which(x>=qx[i] & x<qx[i+1])
		z[index] = jitter(z[index], amount=nx[i])
	}
	return(invisible(z))	
}

'.backend.ccfposteriorgrid' <- function (alt, ref, alpha, q, grid)
{
	f = (alpha * grid)/(2 * (1 - alpha) + alpha * q)
	N = alt + ref
	ll_grid = dbinom(alt, N, f, log=TRUE)
	pr_grid = exp(ll_grid - .backend.LogAdd(ll_grid))
	return(invisible(pr_grid))
}

.CntuEnv <- new.env()
kQ = 15
pi_som_theta_q=c(100, 50, rep(2, kQ-2))
mut_class_w=list("SM"=.5, "GL"=0, "SC"=.5, "OL"=1e-3, "Pi_SM"=15, "Pi_SC"=15)
tM = matrix(c(1,0,0,0,0, 1,1,1,0,0, 0,1,1,1,0, 0,0,1,1,1, 0,0,0,1,1), nrow=5, ncol=5, byrow=TRUE, dimnames=list(c(0:4), c(0:4)))
ColNames = c("Pr_somatic_clonal", "Pr_germline", "Pr_subclonal", "Pr_subclonal_wt0", "Pr_wt0", "Pr_ge2", "Pr_GL_som_HZ_alt", "Pr_GL_som_HZ_ref", "Pr_cryptic_SCNA", "LL")

assign("ColNames", ColNames, envir=.CntuEnv)
assign("tM", tM, envir=.CntuEnv)
assign("mut_class_w", mut_class_w, envir=.CntuEnv)
assign("pi_som_theta_q", pi_som_theta_q, envir=.CntuEnv)
assign("kQ", kQ, envir=.CntuEnv)

rm(ColNames, kQ, pi_som_theta_q, mut_class_w, tM)
