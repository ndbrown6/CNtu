'dilutionCurve' <- function (sq, qt, alpha = seq(from=.1, to=.99, length=99))
{
	fsq = sq*alpha/(alpha*qt + 2*(1-alpha))
	return(invisible(fsq))
}
