\name{nhppSimulate}
\alias{nhppSimulate}
\title{
Simulate a non-homogeneous Poisson Process
}
\description{
This function simulates an NHPP by blocked thinning
}
\usage{
nhppSimulate(smoothRates)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{smoothRates}{
A list containing x and y, which are the mid-points of the window and the smoothed number of events in this window
}
}
\details{
The list component y of the argument represents the smoothed number of events in the window, namely, they represent the window rate
}
\value{
Returns a vector of events of a realization of the NHPP
}
\author{
Jeremy J. Shen
}

\seealso{
\code{\link{nhppRateEstimate}}, \code{\link{nhppSpike}}
}