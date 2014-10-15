\name{ScanStatRefineComp}
\alias{ScanStatRefineComp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Main refining window scan statistics computation
}
\description{
This is a wrapper function to call the C routines for the scan statistic to refine current candidate segmentations computing from the IGS
}
\usage{
ScanStatRefineComp(combZCumSum, combXCumSum, combZPoint, combXPoint, 
  p, nTotal, grid.cur, grid.LR, max.win, statistic)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{combZCumSum}{
A cumulative sum of the number of case reads
}
  \item{combXCumSum}{
A cumulative sum of the number of reads
}
  \item{combZPoint}{
The number of case reads at the grid points
}
  \item{combXPoint}{
The number of reads at the grid points
}
  \item{p}{
The proportion of case reads in the current region
}
  \item{nTotal}{
The total number of reads in the current region
}
  \item{grid.cur}{
The current grid to be computed on
}
  \item{grid.LR}{
The left and right limits of the existing candidate segmentations that will be refined, indexed by the current grid
}
  \item{max.win}{
The maximum inter-window to be considered for new scan
}
  \item{statistic}{
The type of statistic. Can be 'binomial','rabinowitz' or 'normal'.
}
}
\details{
The computations are done in C for speed. The three types of statistics are by the use of exact binomial likelihood ('binomial'), score statistic ('rabinowitz') or using normal approximation to the binomial ('normal').
}
\value{
Returns a matrix containing the refined candidate change points
}
\author{
Jeremy J. Shen
}

\seealso{
\code{\link{ScanCBS}}, \code{\link{ScanIterateGrid}}
}