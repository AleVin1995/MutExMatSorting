MExMaS.HeuristicMutExSorting <- function(mutPatterns, display = TRUE, cluster_rows = FALSE, cluster_cols = FALSE, legend = FALSE,
                                         show_rownames = FALSE, show_colnames = FALSE, col = c('white','blue')){
  
  mutPatterns <- as.matrix(mutPatterns)
  mutPatterns <- sign(mutPatterns)

  ngenes <- nrow(mutPatterns) 
  nsamples <- ncol(mutPatterns)
  
  if (is.null(rownames(mutPatterns))){
    rownames(mutPatterns)  <-  1:ngenes
  }

  if (is.null(colnames(mutPatterns))){
    colnames(mutPatterns)  <-  1:nsamples
  }

  if (nrow(mutPatterns) > 1 & ncol(mutPatterns) > 1){
    
    RowNull <- names(which(rowSums(mutPatterns) == 0))
    RowNonNull <- which(rowSums(mutPatterns) > 0)
    
    ColNull <- names(which(colSums(mutPatterns) == 0))
    ColNonNull <- which(colSums(mutPatterns) > 0)
    
    mutPatterns <- matrix(c(mutPatterns[RowNonNull, ColNonNull]),
                        length(RowNonNull),length(ColNonNull),
                        dimnames=list(rownames(mutPatterns)[RowNonNull],colnames(mutPatterns)[ColNonNull]))
    
    if (nrow(mutPatterns) > 1 & ncol(mutPatterns) > 1){

      coveredGenes <- NA
      uncoveredGenes <- rownames(mutPatterns)
      
      coveredSamples <- NA
      uncoveredSamples <- colnames(mutPatterns)
      BS <- NA
  
      while(length(uncoveredGenes) > 0 & length(uncoveredSamples) > 0){
  
        patterns <- matrix(c(mutPatterns[uncoveredGenes, uncoveredSamples]),
                         nrow = length(uncoveredGenes),
                         ncol = length(uncoveredSamples),
                         dimnames = list(uncoveredGenes, uncoveredSamples))
  
        if(length(uncoveredGenes) > 1){
          bestInClass <- .MExMaS.findBestInClass(patterns)
        }else{
          bestInClass <- uncoveredGenes
        }
  
        if(is.na(BS[1])){
          BS <- bestInClass
        }else{
          BS <- c(BS, bestInClass)
        }
  
        if(is.na(coveredGenes[1])){
          coveredGenes <- bestInClass
        }else{
          coveredGenes <- c(coveredGenes, bestInClass)
        }
  
        uncoveredGenes <- setdiff(uncoveredGenes, coveredGenes)
        toCheck <- matrix(c(patterns[bestInClass, uncoveredSamples]),nrow = 1, ncol=ncol(patterns),
                        dimnames = list(bestInClass, uncoveredSamples))
  
        if (length(coveredGenes) == 1){
          coveredSamples <- names(which(colSums(toCheck) > 0))
        }else{
          coveredSamples <- c(coveredSamples, names(which(colSums(toCheck) > 0)))
        }
  
        uncoveredSamples <- setdiff(uncoveredSamples, coveredSamples)
  
      }
  
      BS <- c(BS, uncoveredGenes)
  
      CID <- .MExMaS.rearrangeMatrix(mutPatterns, BS)
  
      FINALMAT <- mutPatterns[BS, CID]
  
      nullCol <- matrix(0, nrow(FINALMAT),length(ColNull),
                      dimnames = list(rownames(FINALMAT),ColNull))
      
      FINALMAT <- cbind(FINALMAT, nullCol)
      
      nullRow <- matrix(0, length(RowNull),ncol(FINALMAT),
                      dimnames = list(RowNull, colnames(FINALMAT)))
      
      FINALMAT <- rbind(FINALMAT, nullRow)
  
      if (display){
        ## original matrix
        pheatmap::pheatmap(mutPatterns, cluster_rows = cluster_rows, cluster_cols = cluster_cols, legend = legend,
                           show_rownames = show_rownames, show_colnames = show_colnames, main='Original Matrix', col = col)
        ## sorted matrix
        pheatmap::pheatmap(FINALMAT, cluster_rows = cluster_rows, cluster_cols = cluster_cols, legend = legend,
                           show_rownames = show_rownames, show_colnames = show_colnames, main='Sorted Matrix', col = col)
      }
      return(FINALMAT)
      
    } else {
      stop('Matrix must have at least 2 non-null rows and 2 non-null columns')
    }
  } else {
    stop('Matrix must have at least 2 rows and 2 columns') 
  }
}

.MExMaS.findBestInClass <- function(patterns){

  if(nrow(patterns) == 1){
    return(rownames(patterns))
  }

  if(ncol(patterns) == 1){
    return(rownames(patterns)[1])
  }

  exclCov <- colSums(t(2*patterns)-colSums(patterns))
  names(exclCov) <- rownames(patterns)

  return(names(sort(exclCov, decreasing=TRUE))[1])
}

.MExMaS.rearrangeMatrix <- function(patterns, GENES){

  remainingSamples <- colnames(patterns)

  toAdd <- NULL

  DD <- t(t(2*patterns)-colSums(patterns))
  colnames(DD)  <-  remainingSamples

  for (g in GENES){
    cols  <-  remainingSamples[order(DD[g, remainingSamples],decreasing = TRUE)]
    toAdd <- c(toAdd, names(which(patterns[g, cols] > 0)))
    remainingSamples <- setdiff(remainingSamples, toAdd)

    if(length(remainingSamples) == 0){
      break
    }
  }

  toAdd <- c(toAdd, remainingSamples)

  return(toAdd)
}

MExMaS.MEMo <- function(mutPatterns, display = TRUE, cluster_rows = FALSE, cluster_cols = FALSE, legend = FALSE,
                        show_rownames = FALSE, show_colnames = FALSE, col = c('white','blue')){
  
  mutPatterns <- as.matrix(mutPatterns)
  mutPatterns <- sign(mutPatterns)
  
  ngenes <- nrow(mutPatterns) 
  nsamples <- ncol(mutPatterns)
  
  if (is.null(rownames(mutPatterns))){
    rownames(mutPatterns)  <-  1:ngenes
  }
  
  if (is.null(colnames(mutPatterns))){
    colnames(mutPatterns)  <-  1:nsamples
  }
  
  if (nrow(mutPatterns) > 1 & ncol(mutPatterns) > 1){
    
    RowNull <- names(which(rowSums(mutPatterns) == 0))
    RowNonNull <- which(rowSums(mutPatterns) > 0)
    
    ColNull <- names(which(colSums(mutPatterns) == 0))
    ColNonNull <- which(colSums(mutPatterns) > 0)
    
    mutPatterns <- matrix(c(mutPatterns[RowNonNull, ColNonNull]),
                          length(RowNonNull),length(ColNonNull),
                          dimnames=list(rownames(mutPatterns)[RowNonNull],colnames(mutPatterns)[ColNonNull]))
    
    
    if (nrow(mutPatterns) > 1 & ncol(mutPatterns) > 1){
      
      BS <- names(sort(rowSums(mutPatterns), decreasing = TRUE))
      scores <- apply(mutPatterns[BS,], 2, .MExMaS.scoreCol)
      CID <- names(sort(scores, decreasing = TRUE))
      
      FINALMAT <- mutPatterns[BS, CID]
      
      nullCol <- matrix(0, nrow(FINALMAT),length(ColNull),
                        dimnames = list(rownames(FINALMAT),ColNull))
      
      FINALMAT <- cbind(FINALMAT, nullCol)
      
      nullRow <- matrix(0, length(RowNull),ncol(FINALMAT),
                        dimnames = list(RowNull, colnames(FINALMAT)))
      
      FINALMAT <- rbind(FINALMAT, nullRow)
      
      if (display){
        ## original matrix
        pheatmap::pheatmap(mutPatterns, cluster_rows = cluster_rows, cluster_cols = cluster_cols, legend = legend,
                           show_rownames = show_rownames, show_colnames = show_colnames, main='Original Matrix', col = col)
        ## sorted matrix
        pheatmap::pheatmap(FINALMAT, cluster_rows = cluster_rows, cluster_cols = cluster_cols, legend = legend,
                           show_rownames = show_rownames, show_colnames = show_colnames, main='Sorted Matrix', col = col)
      }
      
      return(FINALMAT)
      
    } else {
      stop('Matrix must have at least 2 non-null rows and 2 non-null columns')
    }
  } else {
    stop('Matrix must have at least 2 rows and 2 columns') 
  }
}

.MExMaS.scoreCol <- function(x){
  s <- 2^((length(x)-1):0)
  s <- sum(s*x)
  
  return(s)
}
