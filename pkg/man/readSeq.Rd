\name{readSeq}
\alias{readSeq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Wrapper for managing the reading of different raw data formats
}
\description{
This is a wrapper function. It calls one of the subroutines to reads in a datafile, depending on the format
}
\usage{
readSeq(filename, formatName)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{
The file name of the data file to be read
}
  \item{formatName}{
The format the file is in. Can be either 'Chiang' or 'ELANDPaired'. We recommend using Chiang since this is the minimal required format.
}
}
\details{
We recommend using the 'Chiang' format, as used by the datasets of Chiang (2009). This format requires minimal memory and contains all relevant information for this program. It is a table with two columns, first being the chromosome of the mapped read, and the second being the position of the read in the chromosome. One line for each observation. If one has paired read, please use only one of the reads and the mapped location should be the 5'-end.
}
\value{
  \item{seqF }{Read position for each read}
  \item{seqChr }{Chromosome of each mapped read}
}
\references{
Chiang et al., Nature Methods, 2009, Vol.6 No.1
}
\author{
Jeremy J. Shen
}
\examples{
# This shows the recommended format, the Chiang data format
data(JSSim_NormalSim1)
print(head(JSSim_NormalSim1))
}
\seealso{
\code{\link{readSeq}}, \code{\link{readSeqChiang}}, \code{\link{readSeqELANDPaired}}
}