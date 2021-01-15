\name{MExMaS.HeuristicMutExSorting}
\alias{MExMaS.HeuristicMutExSorting}
\title{Minimal overlap sorting}
\usage{
MExMaS.HeuristicMutExSorting(mutPatterns)
}
\description{
This function implements an heuristic algrorithm that takes in input a sparse binary matrix and sorts its rows and column in a way that the patterns of non null entries have a minimal overalp across rows.
}
\arguments{
 \item{mutPatterns}{numeric binary matrix of the values to be sorted.}
}
\examples{
library(pheatmap)

#Generating a random binary matrix with row and column names
r <- 100
c <- 100
dens<-0.10
mutPatterns <- matrix(0, r, c,dimnames = list(paste('row',1:r,sep=''),paste('col',1:c,sep='')))
mutPatterns[sample(r*c,round(r*c*dens))]<-1

#Executing mutual exclusivity sorting
sortedMat<-MExMaS.HeuristicMutExSorting(mutPatterns)

#visualising original matrix
pheatmap(mutPatterns,cluster_rows = FALSE,cluster_cols = FALSE,legend = FALSE,
  show_colnames = FALSE,show_rownames = FALSE,main='Original Matrix',col=c('white','blue'))
#visualising original matrix
pheatmap(sortedMat,cluster_rows = FALSE,cluster_cols = FALSE,legend = FALSE,
  show_colnames = FALSE,show_rownames = FALSE,main='Sorted Matrix',col=c('white','blue'))
}
