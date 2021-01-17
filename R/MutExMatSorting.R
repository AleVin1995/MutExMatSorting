MExMaS.HeuristicMutExSorting<-function(mutPatterns){

  mutPatterns<-sign(mutPatterns)

  ngenes<-nrow(mutPatterns) 
  nsamples<-ncol(mutPatterns)
  
  if (is.null(rownames(mutPatterns))){
    rownames(mutPatterns) <- 1:ngenes
  }

  if (is.null(colnames(mutPatterns))){
    colnames(mutPatterns) <- 1:nsamples
  }
  
  coveredGenes<-NA
  uncoveredGenes<-rownames(mutPatterns)

  if (length(uncoveredGenes)>1){

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
      toCheck<-matrix(c(patterns[bestInClass,uncoveredSamples]),nrow = 1,ncol=ncol(patterns),dimnames = list(bestInClass,uncoveredSamples))

      if (length(coveredGenes)==1){
        coveredSamples<-names(which(colSums(toCheck)>0))
      }else{
        coveredSamples<-c(coveredSamples,names(which(colSums(toCheck)>0)))
      }

      uncoveredSamples<-setdiff(uncoveredSamples,coveredSamples)

    }

    BS<-c(BS,uncoveredGenes)

    CID<-MExMaS.rearrangeMatrix(mutPatterns,BS)

    FINALMAT<-mutPatterns[BS,CID]

    FINALMAT<-cbind(FINALMAT,nullCol[rownames(FINALMAT),])

    return(FINALMAT)
  }

}
MExMaS.findBestInClass<-function(patterns){

  if(nrow(patterns)==1){
    return(rownames(patterns))
  }

  if(ncol(patterns)==1){
    return(rownames(patterns)[1])
  }

  exclCov<-colSums(t(2*patterns)-colSums(patterns))
  names(exclCov)<-rownames(patterns)

  return(names(sort(exclCov,decreasing=TRUE))[1])
}
MExMaS.rearrangeMatrix<-function(patterns,GENES){

  remainingSamples<-colnames(patterns)

  toAdd<-NULL

  DD<-t(t(2*patterns)-colSums(patterns))
  colnames(DD) <- remainingSamples

  for (g in GENES){
    cols <- remainingSamples[order(DD[g,remainingSamples],decreasing = TRUE)]
    toAdd<-c(toAdd,names(which(patterns[g,cols]>0)))
    remainingSamples<-setdiff(remainingSamples,toAdd)

    if(length(remainingSamples)==0){
      break
    }
  }

  toAdd<-c(toAdd,remainingSamples)

  return(toAdd)
}
