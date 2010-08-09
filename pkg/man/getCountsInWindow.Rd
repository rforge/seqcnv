\name{getCountsInWindow}
\alias{getCountsInWindow}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get number of reads in fixed-width window
}
\description{
Computes the number of reads for each fixed-width window between two limits
}
\usage{
getCountsInWindow(events, startE, endE, windowSize = 10000, sorted = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{events}{
	A vector of the read positions
}
  \item{startE}{
	Left limit
}
  \item{endE}{
	Right Limit
}
  \item{windowSize}{
	Size of the window
}
  \item{sorted}{
	Whether events is sorted, default F
}
}
\details{
	Uses hist() function
}
\value{
A vector of counts for each window
}

\author{
Jeremy J. Shen
}

\examples{
	getCountsInWindow(sample(1:10000, 3000, replace=TRUE), 0, 10000, 100, FALSE)
}

\seealso{
\code{\link{ScanCBS}}
}