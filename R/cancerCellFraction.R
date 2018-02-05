'cancercellFraction' <- function (f_hat, n, qt, q2, alpha, e=0.01)
{
	q_s = rep(NA, length(f_hat))
    s_ppr = somaticPostprob(f_hat, n, qt, q2, w=NULL, sthetaq=NULL, alpha)
    modal_q_s = mutatedCopies(f_hat, n, qt, q2, alpha, e=0.01, d="beta")
    pr_somatic_clonal = s_ppr$ppr[, "Pr_somatic_clonal"]
    c.ix = pr_somatic_clonal > 0.5
    q_s[c.ix] = modal_q_s[c.ix]
    q_s[!c.ix] = 1
  	som_delta = alpha/(2*(1-alpha) + alpha*qt)
  	f_s = q_s * som_delta    
	mut_scale = 1/f_s
	cell_frac = f_hat*mut_scale
	cell_mult = f_hat/som_delta
	ccf_grid = seq(0, 1, by=0.001)
	ccf_dens = matrix(NA, nrow=length(f_hat), ncol=length(ccf_grid))
	ccf_ci95 = matrix(NA, nrow=nrow(ccf_dens), ncol=2)
	ccf_hat = rep(NA, length(f_hat))
  	for (i in 1:length(f_hat)) {
  		if (qt[i]==0) {
  			next
  		}
  		ccf_dens[i,] = CNtu:::.backend.ccfposteriorgrid(round(f_hat[i]*n[i]), round((1-f_hat[i])*n[i]), alpha, qt[i], ccf_grid)
    	ccf_hat[i] = ccf_grid[which.max(ccf_dens[i,])]
    	ecdf = cumsum(ccf_dens[i, ])
    	ccf_ci95[i, ] = approx(x=ecdf, y=ccf_grid, xout=c(0.025, 0.975))$y
  	}
  	nix1 = is.na(ccf_ci95[, 1])
  	ccf_ci95[nix1, 1] = min(ccf_grid)
  	nix2 = is.na(ccf_ci95[,2])
  	ccf_ci95[nix2, 2] = max(ccf_grid)
  	ix = ccf_ci95[, 2] > ccf_grid[length(ccf_grid) - 1]
  	ccf_ci95[ix, 2] = 1.0
  	ix = qt[i] == 0
  	ccf_ci95[ix, ] = NA
  	res = cbind(cell_mult, cell_frac, ccf_hat, ccf_ci95)
  	colnames(res) = c("cell_multiplicity", "scaled_cell_frac", "cancer_cell_frac", "ccf_95CI_low", "ccf_95CI_high")
  	res = cbind(res, s_ppr$ppr)
  	return(invisible(res))
}
