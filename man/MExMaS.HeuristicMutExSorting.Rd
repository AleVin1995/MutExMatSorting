\name{MExMaS.HeuristicMutExSorting}
\alias{MExMaS.HeuristicMutExSorting}
\title{Minimal overlap sorting}
\usage{
MExMaS.HeuristicMutExSorting(mutPatterns, display = TRUE)
}
\description{
This function implements an heuristic algrorithm that takes in input a sparse binary matrix and sorts its rows and column in a way that the patterns of non null entries have a minimal overalp across rows.
}
\arguments{
 \item{mutPatterns}{numeric binary matrix of the values to be sorted.}
 \item{display}{boolean, whether to display the original and sorted matrices. Default is true.}
}
\examples{
library(pheatmap)

# Generating a random binary matrix with row and column names
r <- 100
c <- 100
dens <- 0.10
mutPatterns <- matrix(0, r, c,dimnames = list(paste('row',1:r,sep=''),paste('col',1:c,sep='')))
mutPatterns[sample(r*c,round(r*c*dens))] <- 1

# Executing mutual exclusivity sorting
sortedMat <- MExMaS.HeuristicMutExSorting(mutPatterns)
}