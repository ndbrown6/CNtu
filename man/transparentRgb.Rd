\name{transparentRgb}
\alias{transparentRgb}
\title{
RGB to HTML color code
}
\description{
Converts a character string specifying a color into HTML color code with alpha transparency
}
\usage{
html = transparentRgb(col, alpha)
}
\arguments{
\item{col}{character string specifying a known color in \code{colors(...)}}
\item{alpha}{numeric between 1 and 255 specifying the alpha transparency}
}
\value{
Character string of specified HTML color
}
\author{
\cr
David Brown\cr
\email{brownd7@mskcc.org}\cr
}
\keyword{}