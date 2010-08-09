\name{SegSeqResProcess}
\alias{SegSeqResProcess}
\title{
Read and Process result of SegSeq
}
\description{
Read the segmentation results of SegSeq and returns the change points called
}
\usage{
SegSeqResProcess(filename)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{
The filename of the SegSeq output file to be processed
}
}
\details{
This function is used to read in the SegSeq results and use for performance evaluation and comparison
}
\value{
Return a list the length of unique chromosomes in the result file. For each entry, the label is the chromosome label; and there is a vector of the change point locations called by SegSeq
}

\author{
Jeremy J. Shen
}

\seealso{
\code{\link{nhppSimulate}}
}