\name{readSeqChiang}
\alias{readSeqChiang}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Read data formatted as in Chiang (2009)
}
\description{
Read data formatted as in Chiang (2009), which we recommend using.
}
\usage{
readSeqChiang(filename)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{
The file name of the data set
}
}
\details{
This format requires minimal memory and contains all relevant information for this program. It is a table with two columns, first being the chromosome of the mapped read, and the second being the position of the read in the chromosome. One line for each observation. In case of paired read, we only use the front read (whichever has a smaller position label) and ask that you use only that for input.
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
# This shows the format of this type of data
data(JSSim_NormalSim1)
print(head(JSSim_NormalSim1))
}
\seealso{
\code{\link{readSeq}}, \code{\link{readSeqChiang}}, \code{\link{readSeqELANDPaired}}
}