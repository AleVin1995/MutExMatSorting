## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, include = FALSE---------------------------------------------------
library(MutExMatSorting)

## -----------------------------------------------------------------------------
library(pheatmap)

r <- 100
c <- 100

dens<-0.1
BinMat <- matrix(0, r, c,dimnames = list(paste('row',1:r,sep=''),paste('col',1:c,sep='')))
BinMat[sample(r*c,round(r*c*dens))]<-1

#Executing mutual exclusivity sorting
sortedMat<-MExMaS.HeuristicMutExSorting(BinMat)

#visualising original matrix
pheatmap(BinMat,cluster_rows = FALSE,cluster_cols = FALSE,
         legend = FALSE,show_colnames = FALSE,show_rownames = FALSE,main='Original Matrix',
         col=c('white','blue'))

## -----------------------------------------------------------------------------
#visualising sorted matrix
pheatmap(sortedMat,cluster_rows = FALSE,cluster_cols = FALSE,
         legend = FALSE,show_colnames = FALSE,show_rownames = FALSE,main='Sorted Matrix',
         col=c('white','blue'))

