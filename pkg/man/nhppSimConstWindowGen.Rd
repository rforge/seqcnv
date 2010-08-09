\name{nhppSimConstWindowGen}
\alias{nhppSimConstWindowGen}
\title{
Simulate a Non-Homogeneous PP with constant window spike
}
\description{
Simulate non-homogeneous Poisson processes with a number of constant-widths windows of signal spike, and output the data and meta file
}
\usage{
nhppSimConstWindowGen(controlRates, filename, chromosomeN, nSpike=25, cptLen=c(3,5,8,12,20,30,50,75,100), nPair=2, nRepeat=10, minGain=1.5, maxGain=4, minLoss=0.01, maxLoss=0.5, pGain=0.6)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{controlRates}{
The estimated rate of nhpp for the control
}
	\item{filename}{
The prefix of all the output files from this simulation
}
	\item{chromosomeN}{
The chromosome number. Should be the number from which the samples are emulated
}
  \item{nSpike}{
Number of signal spikes
}
  \item{cptLen}{
Constant width of the signal (CN gain or loss) region to simulate, can be a vector of different values for which to test
}
  \item{nPair}{
Number of tumor samples to generate for each choice of the width of the signal; number of normal samples to generate
}
  \item{nRepeat}{
Number of times to repeat the simulation data generation
}
  \item{minGain}{
Minimal signal gain
}
  \item{maxGain}{
Maximal signal gain
}
  \item{minLoss}{
Minimal signal loss
}
  \item{maxLoss}{
Maximal signal loss
}
  \item{pGain}{
Proportion of the signal regions that are CN gain
}
}
\details{
This function is used in conjunction with a modified, windowed rate vector to simulate non-homogeneous Poisson processes with a number of constant-widths windows of signal spike.
One should use the \code{\link{nhppRateEstimate}} function to estimate the rate of a control sample one wishes to mimic. This function randomly choose windows of a specified constant width, and spike in signals (change points) which can be either gain or loss of copy numbers.
}
\value{
No return value. Generates a number of .txt files, one for each normal/tumor sample as raw data, one input meta file and a file with the true change points for each choice of cptLen.
}

\author{
Jeremy J. Shen
}

\seealso{
\code{\link{nhppSimulate}}
}