\name{MExMaS.Memo}
\alias{MExMaS.Memo}
\title{OncoPrint sorting}
\usage{
MExMaS.Memo(mutPatterns, display = TRUE)
}
\description{
This function implements the sorting algorithm derived from Memo. It takes in input a sparse binary matrix and sorts its rows according to the numbers of non-null entries, whereas columns are sorted through a weighted scoring based on rows ordering.
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
sortedMat <- MExMaS.Memo(mutPatterns)
