\name{boxPlot}
\alias{boxPlot}
\title{
Box plot
}
\description{
Draws a box plots overlaid with data points
}
\usage{
boxPlot(x, main, sub, xlab, ylab, col, file, ...)
}
\arguments{
\item{x}{list or matrix with values to plot}
\item{main}{character string for main title}
\item{sub}{character string for sub title}
\item{xlab}{vector of character strings for labels of individual boxes}
\item{ylab}{character string for y axis label}
\item{col}{character string giving the colors of the dots}
\item{file}{character string specifying the path and file name of the output pdf}
\item{...}{additional arguments that may be passed to \code{boxplot}}
}
\value{
none
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