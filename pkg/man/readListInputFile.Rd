\name{readListInputFile}
\alias{readListInputFile}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Read meta file containing list of raw data files
}
\description{
Reads a meta file that contains the file names and type of the data files. See details for the format.
}
\usage{
readListInputFile(inputFilename, sep = "\t")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inputFilename}{
The name of file, containing relevant information of all input files
}
  \item{sep}{
Delimiter of the meta input file, default is tab-delimited
}
}
\details{
The meta input file should be organized in a table format with 2 columns, one of which is 'file' and the other is 'type', indicating the data file names and whether the data is from 'normal' or 'tumor'.
}
\value{
  \item{normalFiles }{A character vector containing the names of files with the normal reads}
  \item{tumorFiles }{A character vector containing the names of files with the tumor reads}
}
\author{
Jeremy J. Shen
}
\examples{
# This shows the format of the meta file
data(JSSim_Meta)
print(JSSim_Meta)
}
\seealso{
\code{\link{readInput}}, \code{\link{readSeq}}
}