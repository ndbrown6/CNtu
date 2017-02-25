'setupCurves' <- function (nqt = 3, col = c("steelblue", "orange", "salmon"), lty = 1:3, f_hat = NULL, qt = NULL, sq = NULL)
{
	plot(0,0, type="n", axes=FALSE, xlab="", ylab="", frame.plot=FALSE, xlim=c(0,100), ylim=c(0,100))
    box(lwd=2)
    for (i in 1:nqt) {
    	for (j in 1:i) {
    		fsq = dilutionCurve(sq=j, qt=i)
    		lines(1:99, fsq, col=col[i], lwd=2, lty=lty[j])
    	}
    }
    axis(1, at=NULL, cex.axis = 1.5, cex.lab=1.5, las=1)
   	axis(2, at=NULL, cex.axis=1.5, cex.lab=1.5, las=1)
    title(xlab="Cancer cell fraction (%)", ylab="Variant allele fraction (%)", cex.lab=1.75, cex.main=1.75)
    if (!is.null(f_hat)) {
    	if (!is.null(qt) & length(qt)==length(f_hat)) {
    		if (!is.null(sq) & length(sq)==length(qt)) {
    			.backend.addPoints(f_hat, qt, sq, col)
    		}
    	}
    }
}
