'scatterPlot' <- function (x, y, z, main = "", sub = "", xlab = "", ylab = "",
						   col = transparentRgb(col = "black", alpha = 85),
						   fit = list(add = TRUE, type = "lm"),
						   cor = list(add = TRUE, pos = "bottomright"),
						   legend = list(add = FALSE, pos = "topleft", text = ""),
						   file = "", cex = 2, lwd = 2, ...)
{
	if (file!="") {
		pdf(file)
	}
    par(mar=c(6.1, 6.5, 4.1, 1.1))
    if (!missing(z)) {
    	uz = unique(z)
    	nz = length(uz)
    	if (nz>4) {
    		stop("no more than 4 groups allowed")
    	}
    	col = vector(mode="character", length=length(z))
		ucol = c("steelblue", "orange", "salmon", "black")
		for (i in 1:nz) {
			col[z==uz[i]] = transparentRgb(col=ucol[i], alpha=126)
		}
    }
    plot(x, y, col=col, pch=19, cex=cex, axes=FALSE, frame.plot=FALSE, main="", xlab="", ylab="", ...)
    axis(1, at=NULL, cex.axis=1.5, padj=0.25)
    axis(2, at=NULL, cex.axis=1.5, las=1)
    mtext(side=1, text=xlab, line=4, cex=1.5)
    mtext(side=2, text=ylab, line=4, cex=1.5)
    if (sub=="") {
    	title(main=main, cex.main=2.0)
    } else {
    	title(main=paste(main, "\n", sep=""), cex.main=2.0)
    	title(main=paste("\n", sub, sep=""), cex.main=1.5)
    }
    if (fit$add) {
    	if (fit$type=="lm") {
    		abline(lm(y ~ x), col="red")
    	} else if (fit$type=="lqs") {
    		abline(MASS::lqs(y ~ x), col="red")
    	}
    }
    if (cor$add) {
    	z = cor.test(x, y, method="spearman", alternative="two.sided", exact=FALSE)
    	if (cor$pos=="bottomright") {
        	pos.x = min(x, na.rm=TRUE) + .675*(max(x, na.rm=TRUE)-min(x, na.rm=TRUE))
        	pos.y = min(y, na.rm=TRUE) + .1*(max(y, na.rm=TRUE)-min(y, na.rm=TRUE))
        	p.pos.y = pos.y - 0.075*(max(y, na.rm=TRUE)-min(y, na.rm=TRUE))
        	text(pos.x, pos.y, cex=1.5, labels=bquote(rho== .(signif(z$estimate,3))), pos=4)
        	text(pos.x, p.pos.y, cex=1.5, labels=paste("(P = ", toupper(signif(z$p.value,3)), ")", sep=""), pos=4)
        } else if (cor$pos=="topleft") {
        	pos.x = min(x, na.rm=TRUE) + .2*(max(x, na.rm=TRUE)-min(x, na.rm=TRUE))
        	pos.y = min(y, na.rm=TRUE) + .9*(max(y, na.rm=TRUE)-min(y, na.rm=TRUE))
        	p.pos.y = pos.y - 0.075*(max(y, na.rm=TRUE)-min(y, na.rm=TRUE))
        	text(pos.x, pos.y, cex=1.5, labels=bquote(rho== .(signif(z$estimate,3))), pos=4)
        	text(pos.x, p.pos.y, cex=1.5, labels=paste("(P = ", toupper(signif(z$p.value,3)), ")", sep=""), pos=4)
        }	
    }
    if (legend$add) {
    	if (legend$pos=="topleft") {
    		legend("topleft", legend=legend$text, cex=1, pch=19, pt.cex=2.0, col=unique(col), bty="n")
    	} else if (legend$pos=="bottomright") {
    		legend("bottomright", legend=legend$text, cex=1, pch=19, pt.cex=2.0, col=unique(col), bty="n")
    	}
    }
    box(lwd=lwd)
    if (file!="") {
    	dev.off()
    }
}

