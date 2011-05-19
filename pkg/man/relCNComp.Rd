\name{relCNComp}
\alias{relCNComp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute the Relative Copy Number
}
\description{
This computes the relative copy number by each of the segment called
}
\usage{
relCNComp(combX, combZ, tauHatInd, p, alpha)
}
\arguments{
  \item{combX}{
The number of reads at each unique read position
}
  \item{combZ}{
The number of case/tumor reads at each unique read position
}
  \item{tauHatInd}{
The index of change points called
}
  \item{p}{
The overall proportion of case reads
}
	\item{alpha}{
Significance level for testing whether each segment is a gain (relative CN > 1) or loss (relative CN < 1). The method internally corrects for multiple testing.
}
}
\details{
The relative CN is defined as the number of case reads divided by the number of control reads in a window, adjusted for overall proportion of case reads (divided by the overall relative CN).
}
\value{
Returns a vector of relative CN for each of the segment between two change points
}

\author{
Jeremy J. Shen
}
