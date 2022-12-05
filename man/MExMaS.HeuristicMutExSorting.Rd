\name{MExMaS.HeuristicMutExSorting}
\alias{MExMaS.HeuristicMutExSorting}
\title{Minimal overlap sorting}
\usage{
MExMaS.HeuristicMutExSorting(mutPatterns, 
                             display = TRUE, 
                             cluster_rows = FALSE, 
                             cluster_cols = FALSE, 
                             legend = FALSE, 
                             show_rownames = FALSE, 
                             show_colnames = FALSE,
                             col = c('white','blue'))
}
\description{
This function implements an heuristic algrorithm that takes in input a sparse binary matrix and sorts its rows and column in a way that the patterns of non null entries have a minimal overalp across rows.
}
\arguments{
 \item{mutPatterns}{numeric binary matrix of the values to be sorted in the format (features x samples).}
 \item{display}{boolean, whether to display the original and sorted matrices. Default is true.}
 \item{cluster_rows}{boolean values determining if rows should be clustered.}
 \item{cluster_cols}{boolean values determining if columns should be clustered.}
 \item{legend}{logical to determine if legend should be drawn or not.}
 \item{show_rownames}{boolean specifying if column names are be shown.}
 \item{show_colnames}{boolean specifying if column names are be shown.}
 \item{col}{vector of colors used in heatmap.}
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