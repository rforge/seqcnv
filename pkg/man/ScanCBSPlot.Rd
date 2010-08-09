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
ScanCBSPlot(cases, controls, CBSObj, filename, mainTitle, length.out = 1000, smoothF = 25, xlabScale = 10^6, width = 12, height = 18)
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
  \item{length.out}{
The number of windows to use for the display of smoothed rate estimates
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
This function produces three sub-graphs, showing the segmentation calls, the smoothed rate estimate, and the inferred relative copy number. It is crucial that one seperates the plot for each chromosome
}
\value{
No return object
}
\author{
Jeremy J. shen
}

\seealso{
\code{\link{ScanCBS}}, \code{\link{relCNComp}}
}