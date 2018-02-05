'cellularFraction' <- function (f_hat, qt, sq)
{
	alpha = signif(2*(f_hat)/(sq-(f_hat)*(qt-2)), 2)
	return(invisible(alpha))
}
