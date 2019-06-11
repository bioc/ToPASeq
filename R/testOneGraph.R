## Copyright 2010 Laurent Jacob, Pierre Neuvial and Sandrine Dudoit.

## This file is part of DEGraph.

## DEGraph is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.

## DEGraph is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with DEGraph.  If not, see <http://www.gnu.org/licenses/>.

## Modified by Ivana Ihnatova, 2018


testOneGraph <- function(graph, data, classes, useInteractionSigns=TRUE, ...) {
  ## 'graph': an object of class 'graph' or 'graphAM' or 'graphNEL'.
  ## 'data': a 'matrix' (size: number 'p' of genes x number 'n' of samples) of gene expression
  ## 'classes': a 'vector' (length: 'n') of class assignments
  ## 'useInteractionSigns': should the sign of interaction be taken into account ?
  ## '...': further arguments to be passed to testOneConnectedComponent
  
  
  genes <- nodes(graph)
  genes <- intersect(genes, rownames(data))
  if (length(genes) == 0) 
    stop("There is no intersection between expression feature names and the node names on the graph.")
  graph <- subGraph(genes, graph)
  data <- data[genes, , drop = FALSE]
  
  graphList <- .getConnectedComponentList(graph)
  
  gsizes <- sapply(graphList, graph::numNodes)
  w1 <- which(gsizes==1)
  if (length(w1)) {
    graphList <- graphList[-w1]
  }
  res <- lapply(graphList, .testOneConnectedComponent, data, classes, ...)
  res
}


.getConnectedComponentList <- function(graph) {
  nodeList <- RBGL::connectedComp(graph)
  sizes <- sapply(nodeList, length)
  oo <- order(sizes, decreasing=TRUE)
  graphList <- lapply(nodeList[oo], function(x) {
    g<-subGraph(x, graph)
    g@graphData$signMat<-graph@graphData$signMat[x,x]
    return(g)
  }
  )
  graphList
}

.testOneConnectedComponent <- function(graph, data, classes, ..., prop=0.2) {
  ## Check that there is only one connex component
  cc <- RBGL::connectedComp(graph)
  if (length(cc)>1) {
    stop("More than one connex component in graph")
  }
  rm(cc)
  
  p <- graph::numNodes(graph)
  geneNames <- rownames(data)
  
  data<-data[nodes(graph),]
  cls <- sort(unique(classes))
  X1 <- t(data[, classes==cls[1]])
  X2 <- t(data[, classes==cls[2]])
  
  ## get sign information (if any) and infer from its presence the type of graph
  signMat <- graph@graphData$signMat
  if (is.null(signMat)) {
    adjMat <- as(graph,"graphAM")@adjMat
    mat <- adjMat
  } else {
    mat <- signMat
  }
  
  
  ## Argument 'prop'
  if (is.null(prop)) {
    prop <- c(0.05, 0.1, 0.2, 0.4)
    ndim <- unique(c(1:5, ceiling(p*prop)))
    ndim <- sort(ndim[ndim<=p])
  } else {
    prop <- as.numeric(prop, range=c(0, 1))
    ndim <- unique(ceiling(p*prop))
  }
  
  dg <- rowSums(abs(mat))
  lfA <- .laplacianFromA(mat, ..., ltype="meanInfluence", k=1)
  U <- lfA$U

  res <- NULL
  resNames <- NULL
  
  ## T2 in original space
  T2H <- try(.hotellingIPF(t(data), classes, nperm=0, type="cov", equal=TRUE, paired=FALSE))
  if (class(T2H)=="try-error") {
    pH <- NA
  } else {
    pH <- T2H$pval
  }
  res <- c(res, pH)
  resNames <- c(resNames, "T2")
  
  ## T2 in Fourier space
  pU <- rep(NA, length=length(ndim))
  names(pU) <- ndim
  for (kk in seq(along=ndim)) {
    kR <- ndim[kk]
    T2U <- .graph.T2.test(X1, X2, G=NULL, lfA=lfA, k=kR)
    pU[kk] <- T2U$pval
  }
  
  res <- c(res, pU)
  resNames <- c(resNames, paste("T2 (", ndim, " Fourier components)", sep=""))
  
  
  names(res) <- resNames
  ret <- list(p.value=res, graph=graph, k=ndim)
  
  ret
}

.laplacianFromA <- function(A, k=1, ltype=c("meanInfluence", "normalized", "unnormalized", "totalInfluence")) {
  ltype <- match.arg(ltype)
  tol <- 1e-8
  rownames(A) <- NULL
  colnames(A) <- NULL
  
  ltype <- match.arg(ltype)
  
  p <- nrow(A)
  I <- diag(rep(1,p))
  
  if(ltype %in% c("normalized","unnormalized")) # A must be made symmetric
  {
    tsIdx <- ((A == 0) & (t(A) != 0))
    A[tsIdx] <- t(A)[tsIdx]
  }
  
  if(ltype == "normalized")
  {
    iDs <- diag(1/sqrt(rowSums(abs(A))))
    L <- I - (iDs %*% A %*% iDs)
  }
  
  if(ltype == "unnormalized")
  {
    D <- diag(rowSums(abs(A)))
    L <- D - A
  }
  
  if(ltype == "meanInfluence")
  {
    ImA <- diag(as.integer(rowSums(abs(t(A))) != 0)) - diag(1/pmax(1,rowSums(abs(t(A)))))%*%t(A)
    L <- t(ImA)%*%ImA
  }
  
  if(ltype == "totalInfluence")
  {
    ImA <- diag(as.integer(rowSums(abs(t(A))) != 0)) - t(A)
    L <- t(ImA)%*%ImA
  }
  
  edL <- eigen(L, symmetric=TRUE)
  egVal <- rev(edL$values)
  kIdx <- (egVal <= max(egVal[k], tol))
  return(list(U=edL$vectors[,ncol(edL$vectors):1], l=egVal, kIdx=kIdx))
}



.graph.T2.test <- function(X1, X2, G=NULL, lfA=NULL, ..., k=ncol(X1))
{
  tol <- 1e-8
  p <- ncol(X1)
  
  if (is.null(lfA))
  {
    if(is.null(G))
      A <- diag(rep(1,p))
    else
      A <- as(G, "graphAM")@adjMat

    
    lfA <- .laplacianFromA(A, ..., k=k)
  }
  
  U <- lfA$U
  egVal <- lfA$l
  kIdx <- (egVal <= max(egVal[k],tol)) #lfA$kIdx
  rk <- max(which(kIdx)) ## "round" the number of kept eigenvectors to take into account eigenvalue multiplicity
  ## Perform the test in the new basis
  
  X1<<-X1
  U<<-U[, 1:rk, drop=FALSE]
  lfA<<-lfA
  G<<-G
  ut <- .hotellingIPF(rbind(X1%*%U[, 1:rk, drop=FALSE], X2%*%U[, 1:rk, drop=FALSE]), c(rep(1, nrow(X1)), rep(2, nrow(X2))), nperm=0, type="cov", equal=TRUE, paired=FALSE)
  
}
