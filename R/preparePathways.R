

.preparePathways <- function(pwys, method, both.directions=TRUE,
    genes, which="proteins", edgeType=NULL, maxNodes=150, minEdges=0,
    commonTh=2, filterSPIA=FALSE, convertTo="entrez", convertBy=NULL, EdgeAttrs=NULL)
{
    stopifnot(is.list(pwys) || is(pwys, "PathwayList"))


    #konverzia indetifikatorov
    if (is.null(convertBy)) {
      if (convertTo!="none") {
        pwys<-try(sapply(pwys, convertIdentifiers, convertTo), silent=TRUE)
        if (class(pwys)=="try-error") stop ("pathway identifiers could not be converted")
      }
    } else {
      message("Converting identifiers as specified in 'convertBy'")
      pwys<-sapply(pwys, convertIdentifiersByVector, convertBy)
    }
    .CheckNames(pwys, genes)


    N <- length(pwys)
    if (!is.null(maxNodes)) pwys <- .bigPaths(pwys, maxNodes)
    if (!is.null(minEdges)) pwys <- .fewEdges(pwys, minEdges)
    if (!is.null(commonTh)) pwys <- .commonGenes(pwys, which, genes, commonTh)
    
    if (is.null(EdgeAttrs)) EdgeAttrs<-edgeAttrs
    
    pwys <- lapply(pwys, function(p)
      list(.transformPathway(p, method=method, which=which, edgeType=edgeType,
                             both.directions, EdgeAttrs=EdgeAttrs), p@title))

    if (filterSPIA) pwys <- .filterSPIA(pwys)

    nr.filtered <- N - length(pwys)
    if(nr.filtered) message(paste(nr.filtered, "pwys were filtered out"))
    if (length(pwys)==0) stop("No pathways were left for analysis")
    return(pwys)
}



.transformPathway <- function(pwy, method, which, edgeType=NULL,both.directions=TRUE, EdgeAttrs=NULL)
{
  stopifnot(class(pwy) == "Pathway")
  if (is.null(EdgeAttrs))
    EdgeAttrs<-edgeAttrs


  if (method == "TAPPA")
  {
    pwy <- .buildGraphNEL(edges(pwy, which), TRUE, edgeType)
    pwy <- as(pwy, "matrix")
    pwy <- pwy + t(pwy)
    pwy[pwy > 1] <- 1
    pwy[lower.tri(pwy)] <- 0
    diag(pwy) <- 1
    rownames(pwy)<-gsub(".*:","",rownames(pwy))
    colnames(pwy)<-gsub(".*:","",colnames(pwy))
  } else {
    if (method == "SPIA") {
      pwy<-.prepareSPIA2(pwy, which, both.directions, EdgeAttrs)
      pwy<-.getdatp(pwy, EdgeAttrs$beta$rel, EdgeAttrs$beta$beta)
    } else {
      if(method %in% c("PRS", "PWEA", "TopologyGSA","clipper","DEGraphNoSigns", "CePa", "DEGraph"))
      {
        x<-.buildGraphNEL(edges(pwy, which), both.directions, edgeType)
        nodes(x)<-gsub(".*:","",nodes(x))
        
        if(method == "PRS")
          x <- as(x, "matrix")
          if(method == "DEGraph") {
            eA<-merge(EdgeAttrs[[1]], EdgeAttrs[[2]], by.x=2, by.y=1, all=TRUE)
            pos<-as.character(eA[,2][!is.na(eA[,2]) & eA[,"beta"]==1])
            neg<-as.character(eA[,2][!is.na(eA[,2]) & eA[,"beta"]==-1])
            neu<-as.character(eA[,2][!is.na(eA[,2]) & eA[,"beta"]==0])

            signMat<-as(x,"matrix")
            signMat[,]<-0

            e<-edges(pwy,which)
            posedg<-e[,c("src","dest")][e[,"type"] %in% pos, ]
            posind<-nrow(signMat)*(match(posedg[,2], colnames(signMat))-1)+match(posedg[,1], colnames(signMat))
            signMat[posind] <-  1
            negedg<-e[,c("src","dest")][e[,"type"] %in% neg, ]
            negind<-nrow(signMat)*(match(negedg[,2], colnames(signMat))-1)+match(negedg[,1], colnames(signMat))
            signMat[negind] <- -1

            #neuedg<-e[,1:2][e[,4] %in% neu, ]
            #g<-removeEdge(from=neuedg[,1], to=neuedg[,2], g)
            x@graphData$signMat<-signMat
          }
        pwy<-x
      }
    }
  }
  return(pwy)
}

.prepareSPIA2<-function (p, which, both.directions, EdgeAttrs=edgeAttrs){

  es <- graphite::edges(p,which)
  ns <- graphite::nodes(p)
  ns <- gsub(".*:","",ns)
  
  spiaEdges<-EdgeAttrs[[1]]
  if (!all(unique(as.character(es[,6])) %in% spiaEdges[,1] )) stop("Unexpected edge type ", levels(es[,4])[!levels(es[,4]) %in% spiaEdges[,1]], " Please modify edgeAttrs")

  es <- merge(es, spiaEdges, all.x = TRUE)
  l <- sapply(EdgeAttrs[[2]][,1], simplify = FALSE, USE.NAMES = TRUE,
              function(edgeType) {
                est <- es[es[, "spiaType"] == edgeType, , drop = FALSE]
                gnl <- .buildGraphNEL(est, both.directions, NULL)
                nodes(gnl)<-gsub(".*:","",nodes(gnl))
                
                gnl <- graph::addNode(setdiff(ns,graph::nodes(gnl)), gnl, list())
                gnl <- t(as(gnl, "matrix"))
                return(gnl[ns,ns])
              })
  l$title <- p@title
  l$nodes <- ns
  l$NumberOfReactions <- 0

  return(l)

}
.getdatp<-function(x, rel, Beta){
  names(Beta)<-rel
  con<-Reduce("+",Map(function(m, b) {m*abs(sign(b))}, x[rel],as.list(Beta[rel])), init=0  )
  s<- Reduce("+",Map(function(m, b) {m*b}, x[rel],as.list(Beta[rel])), init=0  )
  z = matrix(rep(apply(con, 2, sum), dim(con)[1]), dim(con)[1], dim(con)[1], byrow = TRUE)
  z[z == 0] <- 1
  return(s/z)
}

.filterSPIA<-function(pathways){
  pathways <- Filter(function(p) sum(abs(p))==0, pathways)
  return(pathways)
}

.CheckNames<-function(pathway, expr){

  IDmatchsum<-sapply(pathway, function(x) sum(.nodes(x,"mixed") %in% expr))
  IDmatchmean<-sapply(pathway, function(x) mean(.nodes(x,"mixed") %in% expr))
  if (sum(IDmatchsum)==0)
    stop("Gene labels and node labels do not match. Please, correct your gene identifiers\n",
         paste(utils::head(.nodes(pathway[[1]],"mixed")), collapse=" "), paste(utils::head(expr), collapse=" "))
  cat(sum(IDmatchsum),"node labels mapped to the expression data\n")
  cat("Average coverage", mean(IDmatchmean,na.rm=TRUE)*100,"%\n")
  cat(sum(IDmatchsum==0)," (out of ",length(pathway),") pathways without a mapped node\n", sep="")
}
.nodes<-function(pwy, which){
  e<-graphite::edges(pwy, which)
  unique(c(e$src, e$dest))
}

.preparePerms <- function(de=NULL, all=NULL, x=NULL, group=NULL, nperm=NULL, test=NULL, method)
{
  if(method == "PRS")
  {
    ind <- as.numeric(all %in% names(de))
    perms.ind <- replicate(nperm, sample(ind))
    rownames(perms.ind) <- all
    perms <- apply(perms.ind, 2, function(x) {
      x[x == 1] <- sample(de)
      x
    })
  }
  if (method == "PWEA") {
    perms <- replicate(nperm, test(x, sample(group))$stat)
  }
  if (method=="TopologyGSA") {
    perms<-replicate(nperm, sample(group))
  }
  return(perms)
}

.fewEdges <- function(pwys, minEdges)
    Filter(function(p) nrow(edges(p)) > minEdges, pwys)

.bigPaths <- function(pwys, maxNodes)
    Filter(function(p) length(nodes(p)) <= maxNodes, pwys)

.commonGenes <- function(pwys, which, genes, threshold)
    Filter(function(p) length(intersect(.nodes(p, which), genes)) >= threshold, pwys)

.dagOnly <- function(pwys) Filter(function(p) gRbase::is.DAG(p), pwys)

