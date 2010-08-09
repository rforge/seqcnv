\name{getAutoGridSize}
\alias{getAutoGridSize}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get Automatic Grid Sizes
}
\description{
This produces a default set of grid sizes to be used in Interative Grid Scan
}
\usage{
getAutoGridSize(nL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nL}{
	Number of unique read positions
}
}
\details{
The default grid sizes are powers of 10
}
\value{
 numeric vector of grid sizes
}
\author{
Jeremy J. Shen
}

\examples{
	## Should produce a vector of power of ten up to 10000
	getAutoGridSize(2*10^5)
}
\seealso{
\code{\link{ScanCBS}}
}