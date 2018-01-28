'somaticPostPr' <- function(f_hat, n, qt, q2, alpha, w=NULL, sthetaq=NULL, e=0.01)
{
	if (is.null(w)) {
		w = CNtu::.CntuEnv$mut_class_w
	}
	if (is.null(sthetaq)) {
		sthetaq = CNtu::.CntuEnv$pi_som_theta_q
	}
	PPr = matrix(NA, nrow=length(qt), ncol=length(CNtu::.CntuEnv$ColNames))
	colnames(PPr) = CNtu::.CntuEnv$ColNames
	somq = matrix(NA, nrow=length(qt), ncol=length(sthetaq))
	hkeyz = paste(qt, q2, sep="/")
	nkeyz = table(hkeyz)
	ukeyz = names(nkeyz)
	ukeyz = ukeyz[order(nkeyz, decreasing=TRUE)]
	nkeyz = sort(nkeyz, decreasing=TRUE)
	for (i in 1:length(nkeyz)) {
		tmpk = strsplit(ukeyz[i], "/" )[[1]]
		qtL = as.integer(tmpk[1])
		q2L = as.integer(tmpk[2])
		indx = (qt==qtL) & (q2==q2L)
		if (qtL==0) {
			PPr[indx,] = 0
      		PPr[indx,"LL"] = log(w[["OL"]])
      		somq[indx,] = 0
      		somq[indx,1] = 1
			next
		}
    	sdelta = alpha/(2*(1-alpha) + alpha*(qtL))
    	nL = n[indx]
    	f_hatL = f_hat[indx]
		fsqL = .backend.GetHscnSomaticMutComb(alpha=alpha, qt=qtL, q2L)
		somW = .backend.SomaticMutPrior(W=w, N=nrow(fsqL), thetaq=sthetaq)
		somW =  matrix(somW, ncol=nrow(fsqL), nrow=length(f_hatL), byrow=TRUE)
    	somLogPr = log(somW) + .backend.FhatCombPost(f_hat=f_hatL, n=nL, fsq=fsqL[,"fsq"], e=e)
    	somPr = cbind(somLogPr, matrix(-Inf, nrow=nrow(somLogPr), ncol=(length(sthetaq)-ncol(somLogPr))))
   		somPr = somPr - .backend.LogAdd(somPr)
   		somq[indx,] = exp(somPr)
    	gmlW = .backend.GermlineMutPrior(w)
		if (any(gmlW>0)) {
    		fsqL = .backend.GetHscnGermlineMutComb(alpha=alpha, qt=qtL, q2=q2L)
      		gmlW =  matrix(gmlW, ncol=nrow(fsqL), nrow=length(f_hatL), byrow=TRUE)
      		gmlLogPr = log(gmlW) + .backend.FhatCombPost(f_hat=f_hatL, n=nL, fsq=fsqL[,"fsq"], e=e)
    	} else { 
      		gmlLogPr = matrix(log(0), ncol=length(gmlW), nrow=length(f_hatL))
		}
		sclLogPr = log(w[["SC"]]) + .backend.UnifSubclonalPost(f_hat=f_hatL, sdelta=sdelta, n=nL, e=e)
    	cscnaLogPr = rep(log(0.01), length(f_hatL)) + .backend.CrypticScnaPost(f_hat=f_hatL, sdelta=sdelta, n=nL, e=e)
    	LogLik = cbind(gmlLogPr, cscnaLogPr, sclLogPr, somLogPr)
    	Z = .backend.LogAdd(LogLik)
    	Pr = matrix(exp(LogLik-Z), nrow=nrow(LogLik), ncol=ncol(LogLik))
		PPr[indx,"LL"] = Z
		if (any(!is.finite(Z))) {
      		stop("Non-finite log likelihood")
      	}
    	PPr[indx,"Pr_somatic_clonal"] = 1 - rowSums(Pr[,1:5,drop=FALSE])
    	PPr[indx,"Pr_germline"] = rowSums(Pr[,1:3,drop=FALSE])
    	PPr[indx,"Pr_subclonal"] = Pr[,5,drop=FALSE]
    	PPr[indx,"Pr_cryptic_SCNA"] = Pr[,4,drop=FALSE]
		if (qtL==q2L) {
    		PPr[indx,"Pr_wt0"] = Pr[,ncol(Pr),drop=FALSE]
      		PPr[indx,"Pr_GL_som_HZ_alt"] = Pr[,2]/rowSums(Pr[,1:3,drop=FALSE])
			PPr[indx,"Pr_GL_som_HZ_ref"] = Pr[,1]/rowSums(Pr[,1:3,drop=FALSE])
    	} else {
			PPr[indx,"Pr_wt0"] = 0
		}
    	if (qtL==1) {
    		PPr[indx,"Pr_subclonal_wt0"] = PPr[indx, "Pr_subclonal"]
		} else {
    		PPr[indx,"Pr_subclonal_wt0"] = 0
		}
		if (any(c(qtL-q2L,q2L) > 1)) {
			PPr[indx,"Pr_ge2"] = rowSums(Pr[,7:ncol(Pr),drop=FALSE])
    	}
	}
	sq = rep(NA, nrow(somq))
	if (w[["SM"]] > 0) {
		indx = PPr[,"Pr_somatic_clonal"] < 0.5
		indx[is.na(indx)] = TRUE
		if (any(!indx)) {
			sq[!indx] = apply(somq[!indx,,drop=FALSE], 1, which.max)
		}
		sq[indx] = 1
	}
	PPs = cbind(PPr, sq=sq)
	somPr = rowSums(PPs[,c("Pr_somatic_clonal", "Pr_subclonal"),drop=FALSE], na.rm=TRUE)
	PPs = cbind(PPs, Pr_somatic=somPr)
	return(invisible(list(PPr=PPs, som_mut_Q_tab=somq)))
}
