'dilutionCurve' <- function (sq, qt, alpha = 1:99)
{
	fsq = 100*sq*alpha/(alpha*qt + 2*(100-alpha))
	return(invisible(fsq))
}
