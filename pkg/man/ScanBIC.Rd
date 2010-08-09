\name{ScanBIC}
\alias{ScanBIC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute the modified BIC for change-point models
}
\description{
This computes mBIC for the current change point model. We then use this to determine the appropriate model complexity.
}
\usage{
ScanBIC(combX, combZ, tauHat, lik0, nTotal)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{combX}{
The number of reads at each unique read position
}
  \item{combZ}{
The number of case/tumor reads at each unique read position
}
  \item{tauHat}{
The change points called
}
  \item{lik0}{
The null likelihood. Computed in the main routine.
}
  \item{nTotal}{
The total number of reads
}
}
\details{
This is meanted to be called as a subrountine of \code{\link{ScanCBS}}
}
\value{
Returns a numerical value of mBIC for the current model
}
\author{
Jeremy J. Shen
}

\seealso{
\code{\link{ScanCBS}}
}