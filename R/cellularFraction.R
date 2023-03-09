'cellularFraction' <- function (f_hat, qt, sq)
{
	alpha = signif( (2*f_hat)/(sq + (2*f_hat - qt*f_hat)), 3)
	return(invisible(alpha))
}
