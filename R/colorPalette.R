'colorPalette' <- function (palette = "ubuntu orange")
{
	data('tints')
	tmp = unlist(strsplit(x=palette, split=" ", fixed=TRUE))
	tint = tmp[1]
	palette = tmp[2]
	tint = unlist(strsplit(x=tint, split="", fixed=TRUE))[1]
	palette = unlist(strsplit(x=palette, split="", fixed=TRUE))[1]
	if (paste0(tint, palette)%in%colnames(tints)) {
		return(invisible(as.vector(tints[,paste0(tint, palette)])))
	} else {
		return(invisible(NA))
	}
}

