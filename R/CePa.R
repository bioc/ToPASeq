#' Centrlity-based Pathway entrichment (CePa)
#'
#' The function runs CePa method on microarray or RNA-Seq data. The implementation includes the identification of differentially expressed genes and transformation of pathways' topologies to an appropriate form. Only the ORA version of the CePa method is implemented and covers centralities: equal-weight, in-degree, out-degree, in-reach, out-reach and betweenness.
#'
#' @param x An \code{ExpressionSet} object or a gene expression data matrix or count matrix, rows refer to genes, columns to samples
#' @param group  Name or number of the phenoData column or a character vector or factor that contains required class assigments
#' @param pathways  A list of pathways in a form from \code{graphite} package or created by \code{preparePathways()}
#' @param type Type of the input data, \code{"MA"} for microarray and \code{"RNASeq"} for RNA-Seq
#' @param which Character, which type of nodes is preserved in a pathway. Possible values are \code{"proteins"},\code{"metabolites"},\code{"mixed"}
#' @param edgeType Character, which type of edges is preserved in a pathway. If \code{NULL}, all edges are kept.
#' @param preparePaths  Logical, by default the pathways are transformed with \code{preparePathways()}. Use \code{FALSE}, if you have done this transformation separately
#' @param norm.method Character, the method to normalize RNAseq data. If \code{NULL} then vst-normalization is performed. Possible values are: \code{"edgeR", "vst", "rLog", "none"}
#' @param test.method Character, the method for differentiall expression analysis of RNAseq data. If \code{NULL} then \code{"voomlimma"} is used. Possible values are: \code{"DESeq2", "voomlimma", "vstlimma", "edgeR"}. 
#' @param p.th Numeric, threshold for p-values of tests for differential expression of genes. Use \code{1} if you don't want any threshold to be applied
#' @param logFC.th Numeric, threshold for log fold-change of a gene to identify the gene as differentially expressed. Use negative if you don't want any threshold to be applied
#' @param nperm Numeric, number of permutations
#' @param both.directions,maxNodes,minEdges,commonTh,filterSPIA,convertTo,convertBy Arguments for the \code{preparePathways()}
#'
#'
#' @return A list:
#' \item{res}{A matrix, each row refers to one pathway, each column to one centrality and the value is a p-value.}
#' \item{topo.sig}{A list of weights for genes (nodes) in individual pathways}
#' \item{degtest}{A numeric vector of gene-level differential expression statistics of all genes in the dataset}
#'
#' @references Gu Z., Liu J., Cao K., Zhang J., Wang J.: Centrality-based pathway enrichment: a systematic approach for finding significant pathways dominated by key genes. BMC Systems Biology 2012, 6:56
#'
#' @author Ivana Ihnatova
#' @keywords htest
#'
#' @examples
#' if (require(breastCancerVDX)) {
#' data("vdx")
#' pathways<-pathways("hsapiens","biocarta")[1:3]
#' MAdata<-Biobase::exprs(vdx)[,1:10]
#' rownames(MAdata)<-Biobase::fData(vdx)[,"Gene.symbol"]
#' MAdata<-MAdata[!duplicated(rownames(MAdata)),]
#'
#' CePa(MAdata, Biobase::pData(vdx)[,"er"][1:10], pathways, type="MA", convertTo="SYMBOL")
#' }
#' @export
CePa<-function(x, group, pathways, type,which="proteins", edgeType=NULL, preparePaths=TRUE,  norm.method=NULL, test.method=NULL, p.th=0.05, logFC.th=2, nperm=1000,
               both.directions=TRUE, maxNodes=150, minEdges=0, commonTh=2, filterSPIA=FALSE, convertTo="none", convertBy=NULL){
  degs<-.prepareData(x, group, type, method="CePa", norm.method, test.method, p.th=p.th, logFC.th=logFC.th)
  degs[[1]]<-names(degs[[1]])
  if (preparePaths) paths<-.preparePathways(pathways, method="CePa", both.directions, degs[[2]], which=which, edgeType=edgeType, maxNodes, minEdges, commonTh, filterSPIA, convertTo, convertBy ) else paths<-pathways
  res<-.CePaORA(degs[[1]], degs[[2]], paths, nperm)

  if (type=="MA") {
    gedm<-.processMA(x, group)
    deg.table<-.testMA(gedm[[1]], gedm[[2]])
  } else {
    if (is.null(test.method)) test.method<-"voomlimma"
    if (type=="RNASeq") {
      deg.table<-.testRNAseq(x, group, test.method)
    } else
      deg.table<-NULL
  }

  cen<-c("equal.weight","in.degree","out.degree",
         "in.reach", "out.reach", "betweenness")
  topo.sig<-lapply(paths, function(x) {.CePaWeights(x[[1]], cen)})

  out<-list(res=res, topo.sig=topo.sig, degtable=deg.table)
  class(out)<-c(class(out), "topResultW","topResult")
  return(out)
}
## CePa finalna moja implementacia
# centrality, g=graphNEL


.CePaWeights<-function(x, method){
  .indegree<-function(px){
    return(graph::degree(px)$inDegree)
  }
  
  .outdegree<-function(px){
    return(graph::degree(px)$outDegree)
  }
  
  .inreach<-function(px){
    sp<-RBGL::floyd.warshall.all.pairs.sp(px)
    sp[is.infinite(sp)]<-NA
    out<-apply(sp, 2, max, na.rm=TRUE) #inreach
    return(out)
  }
  
  .outreach<-function(px){
    sp<-RBGL::floyd.warshall.all.pairs.sp(px)
    sp[is.infinite(sp)]<-NA
    out<-apply(sp, 1, max, na.rm=TRUE) #outreach
    return(out)
  }
  
  #.between<-function(px){
  #  return(betweenness(as(px,"matrix"), nodes(px)))
  #}
  
  #Brandes algorithm
  # graph is adjacency matrix
  betweennesR<-function(graph){
    nds<-rownames(graph)
    
    n<-length(nds)
    CB<-setNames(rep(0, n),nds)
    for (s in nds) {
      S<-c()
      P<-lapply(nds, function(x) character())
      names(P)<-nds
      
      sigma<-setNames(rep(0, n), nds)
      sigma[s]<-1
      d<-setNames(rep(-1, n), nds)
      d[s]<-0
      Q<-s
      while(length(Q)>0) {
        v<-Q[1]
        Q<-Q[-1]
        S<-c(v,S)
        for (w in names(which(graph[v,]==1))) {
          if (d[w]<0) {
            Q<-c(Q,w)
            d[w]<-d[v]+1
          } 
          if (d[w]==d[v]+1) {
            sigma[w]<-sigma[w]+sigma[v]
            P[[w]]<-c(P[[w]], v)
          }
        }
      }
      delta<-setNames(rep(0, n),nds)
      while (length(S)>0) {
        w<-S[length(s)]
        S<-S[-length(s)]
        for (v in unique(P[[w]])) delta[v]<- delta[v]+sigma[v]/sigma[w]*(1+delta[w])
        if (w!=s) CB[w]<-CB[w]+delta[w]
      }
      
    }
    
    return(CB)
  }
  
  out<-list()
  length(out)<-length(method)
  names(out)<-method
  if (any(method == "equal.weight")) {
    out[["equal.weight"]]<-setNames(rep(1, length(nodes(x))), graph::nodes(x))
  }
  if (any(method == "in.degree")) {
    out[["in.degree"]]<-.indegree(x)
  }
  if (any(method == "out.degree")) {
    out[["out.degree"]]<-.outdegree(x)
  }
  if (any(method == "degree")) {
    out[["degree"]]<-graph::degree(x, nodes(x))
  }
  if (any(method == "betweenness")) {
    out[["betweenness"]]<-betweennesR(as(x,"matrix"))
  }
  if (any(method == "in.reach")) {
    out[["in.reach"]]<-.inreach(x)
  }
  if (any(method == "out.reach")) {
    out[["out.reach"]]<-.outreach(x)
  }

  out<-lapply(out, function(x) {
    nz<-x>0
    if (sum(nz)>0) return(x+min(x[nz])/100) else return(x)
  })
  return(out)
}

# priprava permutacii
# globalne
#preparePermsCePa<-function(de, all, nperm){
#  pdiff<-length(de)/length(all)
#  rde<-sapply(1:nperm, function(i) all[as.logical(rbinom(length(all),1, pdiff))])
#}


.preparePermsCePap<-function(de, all, nperm, p){
  pdiff<-length(de)/length(all)
  g<-unique(unlist(strsplit(nodes(p), " ", fixed=TRUE)))
  rde<-sapply(seq_len(nperm), function(i) g[as.logical(rbinom(length(g),1, pdiff))])
}


# analyza jednej drahy
.CePaSingle<-function(de, all, weights, rde){
  nperm<-length(rde)
  nod<-names(weights[[1]])
  d<-sapply(nod, function(x) {
    g<-strsplit(x, " ", fixed=T)
    as.numeric(any(g %in% de))
  })
  
  ps<-sapply(weights, function(w) sum(w*d))

  rps<-sapply(rde, function(x) {
    dr<-sapply(nod, function(y) {
      g<-strsplit(y, " ", fixed=TRUE)
      as.numeric(any(g %in% x))
    })

    psr<-sapply(weights, function(w) sum(w*dr))
    return(psr)
  })
  #psr<<-rps
  p<-rowSums(rps>=matrix(ps, nrow=nrow(rps), ncol=ncol(rps)))/nperm
  return(p)
}

.CePaORA<-function(de, all, paths, nperm) {
  cen<-c("equal.weight","in.degree","out.degree",
         "in.reach", "out.reach", "betweenness")
CenList<-lapply(paths, function(x) {.CePaWeights(x[[1]], cen)})

RDE<-lapply(1:length(paths), function(i) .preparePermsCePap(as.character(de), as.character(all), nperm, paths[[i]][[1]]))
OUT<- t(sapply(seq_len(length(RDE)), function(i) .CePaSingle(de, all, CenList[[i]], RDE[[i]])))
rownames(OUT)<-sapply(paths, function(x) x[[2]])
return(OUT)
}
