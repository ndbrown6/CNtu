\name{scatterPlot}
\alias{scatterPlot}
\title{
Scatter plot
}
\description{
Draws a scatterplot of data points with regression fit
}
\usage{
scatterPlot(x, y, z, main, sub, xlab, ylab, col, fit, cor, legend, file, ...)
}
\arguments{
\item{x}{numeric vector of values}
\item{y}{numeric vector of values}
\item{z}{character or numeric vector giving the group memberships of elements in x an y}
\item{main}{character string for main title}
\item{sub}{character string for sub title}
\item{xlab}{character string for x axis label}
\item{ylab}{character string for y axis label}
\item{col}{character string giving the colors of the dots if z is missing}
\item{fit}{list of two giving the type of fitting\cr
		   \enumerate{
    			\item{add:	}{logical to add or not the regression line}
    			\item{type:	}{character string specifying the type of regression 'lm' or 'lqs'}
  			}
}
\item{cor}{list of two giving the placement of correlation coefficient and P-value\cr
		   \enumerate{
    			\item{add:	}{logical to add or not the correlation coefficient}
    			\item{pos:	}{character string for placement 'bottomright' or 'topleft'}
  			}
}
\item{legend}{list of three giving the necessary information for the legend\cr
		   \enumerate{
    			\item{add:	}{logical to add or not the legend}
    			\item{pos:	}{character string for placement 'bottomright' or 'topleft'}
    			\item{text:	}{vector of characters for display}
  			}
}
\item{file}{character string specifying the path and file name of the output pdf}
\item{...}{addition arguments that may be passed to \code{plot}}
}
\value{
None
}
\author{
\cr
David Brown\cr
\email{brownd7@mskcc.org}\cr
}
\keyword{}