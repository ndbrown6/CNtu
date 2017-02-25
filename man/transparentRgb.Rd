\name{transparentRgb}
\alias{transparentRgb}
\title{
RGB to HTML color code
}
\description{
Converts a character string specifying a color into HTML color code with alpha transparency
}
\usage{
htmlcol = transparentRgb()
htmlcol = transparentRgb(col, alpha)
}
\arguments{
\item{col}{character string specifying a known color in \code{colors(...)}}
\item{alpha}{numeric between 1 and 255 specifying the alpha transparency}
}
\value{
character string of specified HTML color
}
\author{
\cr
David Brown\cr
\email{dbrown@ulb.ac.be}\cr
\url{http://homepages.ulb.ac.be/~dbrown/}
}
\seealso{
\code{\link{scatterPlot}}
}
\keyword{}