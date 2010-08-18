\name{BayesCptCI}
\alias{BayesCptCI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bayesian Point-wise Confidence Interval for Change-Point Model
}
\description{
This algorithm computes a point-wise Bayesian CI for the p parameter in change-point model on binomial process.
}
\usage{
BayesCptCI(cases, controls, CBSRes, stepSize="adaptive", adaptMaxMix=80, alpha=0.05, epsilon=10^-4, epsCDF=10^-4, verbose=FALSE)
}
\arguments{
  \item{cases}{
A numeric vector of the case/tumor reads
}
  \item{controls}{
A numeric vector of the control/normal reads
}
  \item{CBSRes}{
Output from the ScanCBS algorithm on this case/control data
}
  \item{stepSize}{
An actual point-wise computation is time-consuming; by using stepSize = n, a Bayesian CI is computed at every n reads. The adaptive option gives good computational speed by choosing stepSize based on the data.
}
	\item{adaptMaxMix}{
An upper bound for the number of unique weights calculated at each change point under the adaptive method. The default is 80 for an average of approximately 5000 mixture components at each point.
}
  \item{alpha}{
Defaults to 0.05 for the usual CI.
}
  \item{epsilon}{
The cutoff for the likelihood ratio of a model with shifted change point compared to the ScanCBS estimated change-point. The likelihood decreases exponentially around the true point or a 'good' estimate of it. Only alternatives with above the cutoff likelihood ratio are considered plausible and integrated in following computations.
}
  \item{epsCDF}{
This is an error tolerance value for finding the quantile of the posterior, which is a Beta mixture distribution.
}
  \item{verbose}{
If \code{TRUE}, then will print much information on each segmentation. For diagnostics only.
}
}
\details{
This method is a Bayesian point-wise CI for our change-point method. It takes model complexity (number of change points) as a given. With the ScanCBS-estimated change points, it evaluates alternatives for the change points around the estimated value and computes the likelihood of the alternative models. Through theoretical derivation, we then have the estimated probability of a case read, p, at a given read index, to have a posterior density given by Beta Mixture. We then compute the quantiles of this distribution using a safe version of Newton-Raphson implemented in C as the CI at this read.
}
\value{
  \item{CIRes }{A matrix containing the location and its CI of p, each column is a location, or a strech of location if stepSize>1}
  \item{wkRes }{The likelihood ratio for alternatives around each estimated change point}
  \item{mixStruct }{The Beta mixture components for each unique location, containing a collection of two shape parameters and the weight wk}
  \item{timeCIRes }{A list containing the result of the timing of this algorithm}
}
\author{
Jeremy J. Shen
}
\seealso{
\code{\link{ScanBIC}}
}