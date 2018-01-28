'boxPlot' <- function (x, main = "", sub  = "", xlab = "", ylab = "", col, file = "",
					   lwd  = 2, ...)
{
	if ("list" %in% is(x)) {
		nboxes = length(x)
		if (nboxes>10) {
			stop("no more than 10 boxes allowed")
		}
	} else if ("matrix" %in% is(x)) {
		nboxes = ncol(x)
		if (nboxes>10) {
			stop("no more than 10 boxes allowed")
		}
	}
	if (missing(col)) {
		col = vector(mode="character", length=nboxes)
		ucol = rainbow_hcl(nboxes, start=30, end=300)
		for (i in 1:nboxes) {
			col[i] = transparentRgb(col=ucol[i], alpha=126)
		}
	}
	if (file!="") {
		pdf(file)
	}
    par(mar=c(6.1, 6.5, 4.1, 1.1))
    boxplot(x, outline=FALSE, main="", xlab="", ylab="", axes=FALSE, frame=FALSE, lwd=1, ...)
    if (any(xlab=="")) {
    	axis(1, at=1:nboxes, labels=rep(xlab, nboxes), cex.axis=1.5, padj=0.25)
    } else {
    	axis(1, at=1:nboxes, labels=xlab, cex.axis=1.5, padj=0.25)
	}
    axis(2, at=NULL, cex.axis=1.5, las=1)
    mtext(side=2, text=ylab, line=4, cex=1.5)
    if (sub=="") {
    	title(main=main, cex.main=2.0)
    } else {
    	title(main=paste(main, "\n", sep=""), cex.main=2.0)
    	title(main=paste("\n", sub, sep=""), cex.main=1.5)
    }
    if ("list" %in% is(x)) {
		for (i in 1:nboxes) {
			points(jitter(rep(i, length(x[[i]])), amount=.25), x[[i]], pch=19, col=col[[i]], cex=2.0)
		}
    } else if ("matrix" %in% is(x)) {
		for (i in 1:nboxes) {
	    		points(jitter(rep(i, nrow(x)), amount=.25), x[,i], pch=19, col=col[i], cex=2.0)
	    }
	}
    if (nboxes==2) {
    	if ("list" %in% is(x)) {
    		p = wilcox.test(x=x[[1]], y=x[[2]], alternative="two.sided", exact=FALSE)$p.value
    	} else if ("matrix" %in% is(x)) {
    		p = wilcox.test(x=x[,1], y=x[,2], alternative="two.sided", exact=FALSE)$p.value
    	}
    } else {
    	if ("list" %in% is(x)) {
    		p = kruskal.test(x=x)$p.value
    	} else if ("matrix" %in% is(x)) {
    		y = list()
    		for (i in 1:nboxes) {
    			y[[i]] = x[,i]
    		}
    		p = kruskal.test(x=y)$p.value
    	}
    }
    mtext(text=paste("P =", toupper(signif(p,3))), side=1, line=4, cex=1.5)
    box(lwd=lwd)
    if (file!="") {
    	dev.off()
    }
}
