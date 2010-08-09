\name{ScanIterateGrid}
\alias{ScanIterateGrid}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Main Scan with Iterative Grid Search
}
\description{
This is a computational speed-up to prevent a quadratic order computation.
}
\usage{
ScanIterateGrid(combX, combZ, combL, statistic, grid.size, nGridSize, timeIGSBreakDown, takeN, verbose, timing)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{combX}{
The number of reads at the unique read positions
}
  \item{combZ}{
The number of case reads at the unique read positions
}
  \item{combL}{
The set of the labels for the unique read positions
}
  \item{statistic}{
The type of statistic to be used. Can be 'binomial','rabinowitz', or 'normal'
}
  \item{grid.size}{
The set of grid sizes for the iterative search. An automatic default can be given
}
  \item{nGridSize}{
The number of grid sizes
}
  \item{timeIGSBreakDown}{
Cumulative timing of IGS, in a broken down fashion
}
  \item{takeN}{
The number of candidate change points to be added to a temporary set at each grid size
}
  \item{verbose}{
If \code{TRUE}, then will print much information on each segmentation. For diagnostics only.
}
  \item{timing}{
If \code{TRUE}, perform a timing of this algorithm, include in the output data file.
}
}
\details{
This algorithm is a computational speed-up tool. It computes the statistic on coarse grids, and refine to finer grids. Also, at each refinement, it computes all new smaller windows on the finer grid that would not have been captured by the coarse grid. Hence it has a New Scan step and a Refine Scan step, both implemented in C for speed. The three types of statistics are by the use of exact binomial likelihood ('binomial'), score statistic ('rabinowitz') or using normal approximation to the binomial ('normal').
}
\value{
  \item{cptsRet }{The current set of change points called after the IGS scan of the current region}
  \item{timeIGSBreakDown }{A break-down of the time used at the stages of the IGS}
}
\author{
Jeremy J. Shen
}

\seealso{
\code{\link{ScanCBS}}, \code{\link{ScanStatNewComp}}, \code{\link{ScanStatRefineComp}}
}