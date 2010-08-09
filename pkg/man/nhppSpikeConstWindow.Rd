\name{nhppSpikeConstWindow}
\alias{nhppSpikeConstWindow}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Spike NHPP rate with constant window width
}
\description{
Randomly spike the smoothed control rate of an NHPP according to the parameters, with constant window width.
}
\usage{
nhppSpikeConstWindow(smoothRates, nSpike = 25, cptLen = 5, minGain = 1.5, maxGain = 10, minLoss = 0.01, maxLoss = 0.5, pGain = 0.6)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{smoothRates}{
The smoothed rate estimate of the control process
}
  \item{nSpike}{
Number of signal spikes
}
  \item{cptLen}{
Window width of each signal region
}
  \item{minGain}{
Minimal Gain relative CN
}
  \item{maxGain}{
Maximal Gain relative CN
}
  \item{minLoss}{
Minimal Loss relative CN
}
  \item{maxLoss}{
Maximal Loss relative CN
}
  \item{pGain}{
Proportion of signal regions that are CN gain
}
}
\details{
The signal strength are uniformly drawn between the two limits.
}
\value{
  \item{spikeMat }{A matrix containing the actual signal spike information}
  \item{caseRates }{The rate of the case PP to be simulated after signal spike}
}
\author{
Jeremy J. Shen
}

\seealso{
\code{\link{nhppSimulate}}
}