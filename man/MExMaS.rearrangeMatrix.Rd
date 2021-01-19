\name{MExMaS.rearrangeMatrix}
\alias{MExMaS.rearrangeMatrix}
\title{Rearrange matrix columns to minimise row-wise entry overlap}
\usage{
MExMaS.rearrangeMatrix(patterns,GENES)
}
\description{
This function rearranges the binary matrix columns in order to minimise row-wise entry overlap based on exclusive coverage.
}
\arguments{
 \item{patterns}{numeric binary matrix of the values to be sorted, after removing samples with no entries.}
 \item{GENES}{character vector containing rownames ordered according to exclusive coverage.}
}
\examples{
#Generating a random binary matrix with row and column names
r <- 100
c <- 100
dens<-0.10
mutPatterns <- matrix(0, r, c,dimnames = list(paste('row',1:r,sep=''),paste('col',1:c,sep='')))
mutPatterns[sample(r*c,round(r*c*dens))]<-1

#Rowwise sorting
nsamples<-ncol(mutPatterns)

coveredGenes<-NA
uncoveredGenes<-rownames(mutPatterns)

idNull<-which(colSums(mutPatterns)==0)
nullCol<-matrix(c(mutPatterns[,idNull]),nrow(mutPatterns),
                length(idNull),dimnames = list(rownames(mutPatterns),colnames(mutPatterns)[idNull]))

idNonNull<-which(colSums(mutPatterns)>0)
mutPatterns<-matrix(c(mutPatterns[,idNonNull]),
                    nrow(mutPatterns),length(idNonNull),
                    dimnames=list(rownames(mutPatterns),colnames(mutPatterns)[idNonNull]))

coveredSamples<-NA
uncoveredSamples<-colnames(mutPatterns)
BS<-NA

while(length(uncoveredGenes)>0 & length(uncoveredSamples)>0){

  patterns<-matrix(c(mutPatterns[uncoveredGenes,uncoveredSamples]),
                   nrow = length(uncoveredGenes),
                   ncol = length(uncoveredSamples),
                   dimnames = list(uncoveredGenes,uncoveredSamples))

  if(length(uncoveredGenes)>1){
    bestInClass<-MExMaS.findBestInClass(patterns)
  }else{
    bestInClass<-uncoveredGenes
  }

  if(is.na(BS[1])){
    BS<-bestInClass
  }else{
    BS<-c(BS,bestInClass)
  }

  if(is.na(coveredGenes[1])){
    coveredGenes<-bestInClass
  }else{
    coveredGenes<-c(coveredGenes,bestInClass)
  }

  uncoveredGenes<-setdiff(uncoveredGenes,coveredGenes)
  toCheck<-matrix(c(patterns[bestInClass,uncoveredSamples]),
    nrow = 1,ncol=ncol(patterns),
    dimnames = list(bestInClass,uncoveredSamples))

  if (length(coveredGenes)==1){
    coveredSamples<-names(which(colSums(toCheck)>0))
  }else{
    coveredSamples<-c(coveredSamples,names(which(colSums(toCheck)>0)))
  }

  uncoveredSamples<-setdiff(uncoveredSamples,coveredSamples)

}

GENES<-c(BS,uncoveredGenes)

## Columnwise sorting
CID<-MExMaS.rearrangeMatrix(mutPatterns,GENES)
CID
}
\keyword{MutExMatSorting}
