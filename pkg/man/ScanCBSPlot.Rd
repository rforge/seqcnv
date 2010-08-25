\name{ScanCBSPlot}
\alias{ScanCBSPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Main Plotting of the scan statistic segmentation
}
\description{
This is an overall plotting function to display the segmentation for a chromosome
}
\usage{
ScanCBSPlot(cases, controls, CBSObj, filename, mainTitle, CIObj=NULL, length.out=10000, localWindow=0.5*10^5, localSeparatePlot=TRUE, smoothF=0.025, xlabScale=10^6, width=12, height=18)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cases}{
The case read positions (should be restricted to a chromosome)
}
  \item{controls}{
The control read positions (should be restricted to a chromosome)
}
  \item{CBSObj}{
The output object of the \code{\link{ScanCBS}} function
}
  \item{filename}{
The output file names of the plot
}
  \item{mainTitle}{
The title of the plot
}
  \item{CIObj}{
Optional; the Bayesian CI computed by \code{BayesCptCI} function
}
  \item{length.out}{
The number of windows to use for the display of smoothed rate estimates
}
  \item{localWindow}{
The number of genome locations to show around each of the called change points
}
	\item{localSeparatePlot}{
Whether to show the local behavior of each change point in a seperate PDF file. Default to TRUE. The output file are the given filename attached with the index and actual location of the change point.
}
  \item{smoothF}{
The lowess smoothing factor. The proportion of windows around the current window that affects its smoothed rate estimate
}
  \item{xlabScale}{
The scaling factor of the read positions, often in 10^6, or Mb
}
  \item{width}{
The width of the output graph in inches
}
  \item{height}{
The height of the output graph in inches
}
}
\details{
This function produces three sub-graphs, showing the segmentation calls, the smoothed rate estimate, and the inferred relative copy number. It is crucial that one seperates the plot for each chromosome. It also makes a zoom-in plot for a region around each of the called change points.
}
\value{
No return object
}
\author{
Jeremy J. shen
}

\seealso{
\code{\link{ScanCBS}}, \code{\link{ScanCBSSimPlot}}, \code{\link{relCNComp}}
}