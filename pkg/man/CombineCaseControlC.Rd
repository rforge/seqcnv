\name{CombineCaseControlC}
\alias{CombineCaseControlC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Combine case and control reads
}
\description{
Combine the case and control reads; finds the unique read positions and count the number of case and control reads.
}
\usage{
CombineCaseControlC(cases, controls)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cases}{
	A vector of numeric read positions from case sample
}
  \item{controls}{
	A vector of numeric read positions from control sample
}
}
\details{
A few C functions are used for efficient implementation
}
\value{
  \item{combX}{Number of total reads at read position}
  \item{combZ}{Number of case reads at read position}
 	\item{combL}{Vector of unique read positions}

}
\author{
Jeremy J. Shen
}
\seealso{
\code{\link{ScanCBS}}
}