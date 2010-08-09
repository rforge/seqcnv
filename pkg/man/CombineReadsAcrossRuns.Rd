\name{CombineReadsAcrossRuns}
\alias{CombineReadsAcrossRuns}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Combine multiple read lists
}
\description{
Combines multiple lists in the same format of the same sample into one list of the said format.
}
\usage{
CombineReadsAcrossRuns(seqs)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{seqs}{
	A list of lists, each containing equal number of numeric vectors that can be concatenated together. Both number of lists and number of variables can be arbitrary.
}
}
\value{
Returns a list of the same format as the input lists
}
\examples{
data(JSSim_NormalSim1)
data(JSSim_NormalSim2)
write.table(JSSim_NormalSim1, file="JSSim_NormalSim1.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(JSSim_NormalSim2, file="JSSim_NormalSim2.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
JSSim_Normal1 = readSeqChiang("JSSim_NormalSim1.txt")
JSSim_Normal2 = readSeqChiang("JSSim_NormalSim2.txt")
file.remove(c("JSSim_NormalSim1.txt", "JSSim_NormalSim2.txt"))
combJSNormal = CombineReadsAcrossRuns(list(JSSim_Normal1, JSSim_Normal2))
print(c(length(JSSim_Normal1$seqF), length(JSSim_Normal2$seqF), length(combJSNormal$seqF)))
}
\author{
Jeremy J. Shen
}
\seealso{
\code{\link{ScanCBS}}
}