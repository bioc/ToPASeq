#' Pathway Regulation Score (PRS)
#'
#' A function runs PRS method on a gene expression data matrix or count matrix and vector dividing samples into two groups and a set of pathways from \code{graphite} package. The PRS method (please see Reference for the details) was adapted to \code{graphite}'s graphs where each node is represented only by one gene.
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
#' @param p.th Numeric, threshold for p-values of tests for differential expression of genes. Use \code{1} if you don't want any threshold to be applied
#' @param logFC.th Numeric, threshold for log fold-change of a gene to identify the gene as differentially expressed. Use negative if you don't want any threshold to be applied
#' @param nperm Numeric, number of permutations
#' @param both.directions,maxNodes,minEdges,commonTh,filterSPIA,convertTo,convertBy Arguments for the \code{preparePathways()}
#'
#' @return A list:
#' \item{res}{A data frame with normalized score, p-value and FDR-adjusted p-value for each pathway}
#' \item{topo.sig}{A list with log fold-changes and number of downstream differentially expressed nodes for nodes of individual pathways}
#' \item{degtest}{A named vector of statistics from testing the differential expression of genes}
#'
#' @references Maysson Al-Haj Ibrahim, Sabah Jassim, Michael Anthony Cawthorne, and Kenneth Langlands. A Topology-Based Score for Pathway Enrichment, Journal of Computational Biology. May 2012, 19(5): 563-573
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
#' PRS_wrapper(MAdata, Biobase::pData(vdx)[,"er"][1:10], pathways, type="MA", convertTo="SYMBOL", logFC.th=-1, nperm=100)
#' }
#' @export
PRS_wrapper<-function(x, group, pathways, type,which="proteins", edgeType=NULL,  preparePaths=TRUE, norm.method=NULL, test.method=NULL, p.th=0.05, logFC.th=2, nperm=1000,
              both.directions=TRUE, maxNodes=150, minEdges=0, commonTh=2, filterSPIA=FALSE, convertTo="none", convertBy=NULL){
  degs<-.prepareData(x, group, type, method="PRS", norm.method, test.method, p.th=p.th, logFC.th=logFC.th)
  if (preparePaths) paths<-.preparePathways(pathways, method="PRS", both.directions, degs[[2]], which=which, edgeType=edgeType,maxNodes, minEdges, commonTh, filterSPIA, convertTo, convertBy ) else paths<-pathways
  res<-.prs(abs(degs[[1]]), degs[[2]], paths, nperm)
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

  topo.sig<-.collectWeightsPRS(degs[[1]], degs[[2]], paths)
  out<-list(res=res, topo.sig=topo.sig, degtable=deg.table)
  class(out)<-c(class(out), "topResultW","topResult")
  return(out)
}
#PRS

# Priprava permutacii
# all mena vsetkych
# de mena a fold-change DEG
# out matica nahodnych fold-change
#preparePermsPRS<-function(all, de){
#ind<-as.numeric(all %in% names(de))
#perms.ind<-replicate(nperm, sample(ind))
#rownames(perms.ind)<-all
#perms<-apply(perms.ind, 2, function(x) {x[x==1]<-sample(de);x})
#return(perms)
#}

# PRS skore pre jednu drahu
.PRSSingle<-function(path, de, all, perms){


weight<-.PRSweights(path, de, all)
#g<-rownames(path)[rownames(path) %in% all]
#ind<-g %in% names(de)
#nf<-sum(ind)/length(g)
#expr<-ifelse(g %in% names(de), de[g], ifelse(g %in% all, 1, 0 ))

g<-rownames(path)
g<-lapply(g, function(x) strsplit(x, " ")[[1]])
g<-lapply(g, function(x) substr(x, regexpr(":",x)+1, nchar(x)))
ind<-sapply(g, function(x) any(as.character(x) %in% all))

indde<-sapply(g, function(x) any(as.character(x) %in% names(de)))
nf<-sum(indde)/length(g)
g<-g[ind ]
expr<-sapply(g, function(x) if (any(as.character(x) %in% names(de))) max(de[as.character(x)], na.rm=TRUE) else if (any(as.character(x) %in% all)) 1 else 0)
obs<-sum(expr*weight)*nf

#random
gs<-unlist(g)[unlist(g) %in% all]
permsub<-perms[gs,, drop=FALSE]
gn<-unlist(sapply(seq_len(length(g)), function(i) rep(i, sum(g[[i]] %in% all))))
permsub<-apply(permsub, 2, function(x) tapply(x, gn, max, na.rm=T))


weight.rn<-apply(permsub,2, function(x) {

 weight<-setNames(rep(0, length(g)),rownames(path)[ind])
 if (length(g)>0 & length(g[x!=0])>=1)  weight[rownames(path)[ind][x!=0]]<-  downstreamCpp(path[ind,ind], rownames(path)[ind], rownames(path)[ind][x!=0])+1 #set, rownames(path)[ind1], rownames(path)[ind1][ind]
 weight[x==0]<-0
 return(weight)
})

nf.rn<-colSums(permsub !=0)

rand<-colSums(weight.rn*permsub)*(nf.rn/length(g))

# normalization
obs<-(obs-mean(rand))/sd(rand)
rand<-(rand-mean(rand))/sd(rand)
p.value<-sum(rand >= obs)/length(rand)
res<-c(nPRS=obs, p.value=p.value)
return(res)
}

.PRSweights<-function(path, de, all){

 g<-rownames(path)
 g<-lapply(g, function(x) strsplit(x, " ")[[1]])
 g<-lapply(g, function(x) substr(x, regexpr(":",x)+1, nchar(x)))

 ind1<-sapply(g, function(x) any(as.character(x) %in% all))
 
 g<-g[ind1]

 if (length(g)==0) stop("Gene names in a pathway and in data don't match")
 set<-path[ind1,ind1]

  ind<-sapply(g, function(x) any(as.character(x) %in% names(de)))
# ind<-g %in% names(de)
 if (length(g)>0 & sum(ind)>=1) weight<-downstreamCpp(set, rownames(path)[ind1], rownames(path)[ind1][ind]) +1 else weight<-setNames(rep(0, length(g)),rownames(path)[ind1])
wei<-setNames(rep(0, length(g[ind1])),rownames(path)[ind1])
 wei[names(weight)]<-weight
 return(wei)
 }

.prs<-function(de, all, pathways, nperm){

perms<-.preparePerms(all=all, de=de, nperm=nperm, method="PRS")

out<-.catchErr(pathways, function(p) .PRSSingle(p[[1]], de, all, perms))

if(length(out[[1]])>0){
out[[1]]<-data.frame(t(vapply(out[[1]], function(x) x, numeric(2))))
out[[1]]$q.value<-p.adjust(out[[1]]$p.value,"fdr") }

return(out)
}

.collectWeightsPRS<-function(de, all, pathways){
out<-.catchErr(pathways, function(p) .PRSweights(p[[1]], de, all))
return(out[[1]])
}
