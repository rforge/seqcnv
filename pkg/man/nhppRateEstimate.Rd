\name{nhppRateEstimate}
\alias{nhppRateEstimate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate the rate of non-homogeneous PP with data
}
\description{
Given a vector of point events, give a rough estimate of the rate of underlying non-homogeneous Poisson process by window and smoothing
}
\usage{
nhppRateEstimate(controls, length.out = floor(length(controls)/20), lowessF = 0.1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{controls}{
A vector of point locations (read positions) of a control sample for which the rate is wanted
}
  \item{length.out}{
The number of windows to be used for the rate estimate vector; default to be number of observations/100
}
  \item{lowessF}{
Smoothing factor for the lowess smoothing of the windowed rates, describes the proportion of windows around a particular window that has influence on its smoothed rate estimate
}
}
\details{
This is used to give a realistic estimate of the rate nhpp of control samples
}
\value{
Returns a vector of length \code{length.out} that contains the smoothed rate estimate of each window
}

\author{
Jeremy J. Shen
}

\seealso{
\code{\link{nhppSimulate}}
}