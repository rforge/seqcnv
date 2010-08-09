\name{readSeqELANDPaired}
\alias{readSeqELANDPaired}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Read raw data formatted as in paired ELAND output
}
\description{
Read datasets with paired-end format, possible output format of ELAND
}
\usage{
readSeqELANDPaired(filename)
}
\arguments{
  \item{filename}{
The file name of the data set
}
}
\details{
This format has two reads per line, each looking like 'NACGATGAAACCCCGTCTCTACTAACCATACAAAAA hs_ref_chr17.fa 12091150 R	TGTCGCCCAGGCTGCAATGCAGTGGCGCGATCTCGG hs_ref_chr17.fa 12091018 F'. There are 8 columns, 4 for each of the paired read. The first is the actual read sequence, which we discard; the second is the chromosome of the mapped read; the third is the read position; and the last is indicating whether it is a front- or rear- end read. We only use the reads with the same mapped chromosome and only the front read. This contains more information than needed; the Chiang format is prefered.
}
\value{
  \item{seqF }{Read position for each read}
  \item{seqChr }{Chromosome of each mapped read}
}
\author{
Jeremy J. Shen
}

\seealso{
\code{\link{readSeq}}, \code{\link{readSeqChiang}}
}