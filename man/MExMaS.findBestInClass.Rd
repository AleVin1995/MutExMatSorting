\name{MExMaS.findBestInClass}
\alias{MExMaS.findBestInClass}
\title{Find gene with the highest exclusive coverage}
\usage{
MExMaS.findBestInClass(patterns)
}
\arguments{
 \item{patterns}{numeric binary matrix of the values to be sorted, after removing samples with no entries.}
}
\description{
This function finds the gene (i.e. row) with the highest exclusive coverage. The exclusive coverage for a gene g is defined as the number of uncovered samples in which this gene is mutated minus the number of samples in which at least another uncovered gene is mutated.
}
\examples{
#Generating a random binary matrix with row and column names
r <- 100
c <- 100
dens<-0.10
mutPatterns <- matrix(0, r, c,dimnames = list(paste('row',1:r,sep=''),paste('col',1:c,sep='')))
mutPatterns[sample(r*c,round(r*c*dens))]<-1

#Removing samples with no entries
idNull<-which(colSums(mutPatterns)==0)
nullCol<-matrix(c(mutPatterns[,idNull]),nrow(mutPatterns),
                length(idNull),dimnames = list(rownames(mutPatterns),colnames(mutPatterns)[idNull]))

idNonNull<-which(colSums(mutPatterns)>0)
patterns<-matrix(c(mutPatterns[,idNonNull]),
                    nrow(mutPatterns),length(idNonNull),
                    dimnames=list(rownames(mutPatterns),colnames(mutPatterns)[idNonNull]))

#Find gene that maximise exclusive coverage
BS<-MExMaS.findBestInClass(patterns)
BS
}
\keyword{MutExMatSorting}
