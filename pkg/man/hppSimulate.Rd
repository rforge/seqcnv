\name{hppSimulate}
\alias{hppSimulate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulate a homogeneous Poisson Process
}
\description{
Simulation of a homogeneous poisson process using poisson and uniform distribution
}
\usage{
hppSimulate(lambda, maxVal)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lambda}{
	The rate of the poisson process
}
  \item{maxVal}{
	The maximum length of the process to be observed
}
}
\details{
This is a very simple simulation function meant to be used in the NHPP generation.
}
\value{
	Returns a vector of point events generated for this process
}

\author{
Jeremy J. Shen
}

\seealso{
\code{\link{nhppSimulate}}

}