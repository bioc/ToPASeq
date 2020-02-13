#' Differential Expression of Graph (DEGraph)
#'
#' DEGraph implements recent hypothesis testing methods which directly assess whether a particular gene network is differentially expressed between two conditions. In employs Graph Laplacian, Fourier transformation and multivariate T2-statistic
#'
#' @param x An \code{ExpressionSet} object or a gene expression data matrix or count matrix, rows refer to genes, columns to samples
#' @param group  Name or number of the phenoData column or a character vector or factor that contains required class assigments
#' @param pathways  A list of pathways in a form from \code{graphite} package or created by \code{preparePathways()}
#' @param type Type of the input data, \code{"MA"} for microarray and \code{"RNASeq"} for RNA-Seq
#' @param which Character, which type of nodes is preserved in a pathway. Possible values are \code{"proteins"},\code{"metabolites"},\code{"mixed"}
#' @param edgeType Character, which type of edges is preserved in a pathway. If \code{NULL}, all edges are kept.
#' @param preparePaths  Logical, by default the pathways are transformed with \code{preparePathways()}. Use \code{FALSE}, if you have done this transformation separately
#' @param norm.method Character, the method to normalize RNAseq data. If \code{NULL} then vst-normalization is performed. Possible values are: \code{"edgeR", "vst", "rLog", "none"}
#' @param test.method Character, the method for differentiall expression analysis of RNAseq data. If \code{NULL} then \code{"voomlimma"} is used. Possible values are: \code{"DESeq2", "voomlimma", "vstlimma", "edgeR"}. This analysis is needed only for the visualization.
#' @param overall Character, how should the overall p-value for a pathway be calculated. The possible values are: "mean", "min", "biggest". "biggest" returns the p-value of the biggest connected component.
#' @param useInteractionSigns Logical, should types of interaction be included in the analysis?
#' @param EdgeAttrs A list containing two data.frames. See \code{edgeData} for the details. The interactions are assigned signs according to the \code{beta} column of the second data.frame.
#' @param both.directions,maxNodes,minEdges,commonTh,filterSPIA,convertTo,convertBy Arguments for the \code{preparePathways()}
#'
#' @return A list:
#' \item{res}{Results from analysis of individual pathways. The first column refers to the overall p-value for a pathway. Then groups of four columns follows. One group refers to one connected component and contains a pair of p-values (without and with Fourier transformation), graph and number of Fourier componets used in the test. The number of groups is equal to the highest number of components in analysed pathways. Components are sorted in the decreasing order of their nodes number.}
#' \item{topo.sig }{\code{NULL}, present for the compatibility with outputs from other methods}
#' \item{degtest}{A data.frame of gene-level statistics of all genes in the dataset} A list:
#'
#' @references L. Jacob, P. Neuvial, and S. Dudoit. Gains in power from structured two-sample tests of means on graphs. Technical Report arXiv:q-bio/1009.5173v1, arXiv, 2010.
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
#' DEGraph(MAdata, Biobase::pData(vdx)[,"er"][1:10], pathways, type="MA", convertTo="SYMBOL")
#' }
#' @export
DEGraph<-function(x, group, pathways, type, which="proteins", edgeType=NULL,preparePaths=TRUE, norm.method=NULL, test.method=NULL, overall="biggest", useInteractionSigns=TRUE, EdgeAttrs=NULL,
                  both.directions=TRUE, maxNodes=150, minEdges=0, commonTh=2, filterSPIA=FALSE, convertTo="none", convertBy=NULL){
  gedm<-.prepareData(x, group, type, method="DEGraph", norm.method)
  if (preparePaths) {
    if (useInteractionSigns) {
      paths<-.preparePathways(pathways, method="DEGraph", both.directions,rownames(gedm[[1]]), which=which, edgeType=edgeType,maxNodes, minEdges, commonTh, filterSPIA, convertTo, convertBy, EdgeAttrs)
    } else
      paths<-.preparePathways(pathways, method="DEGraphNoSigns", both.directions, rownames(gedm[[1]]),which=which, edgeType=edgeType, maxNodes, minEdges, commonTh, filterSPIA, convertTo, convertBy)
  } else paths<-pathways
  
  res<-.degraph(gedm[[1]], gedm[[2]], paths, overall)

  if (type=="MA") deg.table<-.testMA(gedm[[1]], gedm[[2]])
  if (is.null(test.method)) test.method<-"voomlimma"
  if (type=="RNASeq") deg.table<-.testRNAseq(x, group, test.method)

  out<-list(res=res, topo.sig=NULL, degtable=deg.table)
  class(out)<-c(class(out), "topResultE","topResult")
  return(out)
}


.processDEGraph<-function(res, overall){
if (length(res)>0) {
   nc<-max(unlist(sapply(res, length)))
   p<-suppressWarnings(t(sapply(res, function(x) sapply(1:nc, function(i) as.numeric(try(x[[i]][[1]][1], TRUE))))))
   pFourier<-suppressWarnings(t(sapply(res, function(x) sapply(1:nc, function(i) as.numeric(try(x[[i]][[1]][2], TRUE))))))
   graphs<-suppressWarnings(t(sapply(res, function(x) sapply(1:nc, function(i)
   if ( any("graphNEL" == is((try(x[[i]][[2]], TRUE))))) x[[i]][[2]] else NA
 ))) )
   colnames(graphs)<-paste("Comp",1:nc,".","graph", sep="")
   k<-suppressWarnings(t(sapply(res, function(x) sapply(1:nc, function(i) as.numeric(try(x[[i]][[3]], TRUE))))))

   ord<-as.numeric(sapply(1:nc, function(i) seq(i,3*nc, nc)))

out<-cbind(p, pFourier, k, deparse.level=0)[,ord, drop=FALSE]
tmp<-c("p","pFourier", "k")

colnames(out)<-paste("Comp",rep(1:nc, each=3),".",tmp, sep="")



if (!(any(overall==c("min","mean","biggest")))) stop("Invalid value for 'overall'. Use one of the following: 'min','mean','biggest'")
if (overall=="min") f<-min
if (overall=="mean") f<-mean
if (overall=="biggest") f<-function(x, na.rm=TRUE) {x[1]}
overalp<-apply(pFourier, 1, function(x) f(x,na.rm=TRUE))

out<-cbind(Overall.p=overalp, out)
out<-list(data.frame(out), graphs=graphs)
}  else out<-res
return(out)
}

.degraph<-function(exprs, group,pathways,  overall="biggest"){

out<-.catchErr(pathways, function(p) {
testOneGraph(p[[1]], exprs, group, useInteractionSigns = FALSE)
})

if (length(out[[1]])==0) stop(print(out), print(head(rownames(exprs))), print(pathways[[1]]), "No pathway was sucessfully tested")


out[[1]]<-.processDEGraph(out[[1]], overall)
out[[1]][[1]]$Overall.q.value<-p.adjust(out[[1]][[1]][,"Overall.p"],"fdr")
out[[1]][[1]]<-out[[1]][[1]][,c(1, ncol(out[[1]][[1]]), 2:(ncol(out[[1]][[1]])-1))]
return(out)
}


