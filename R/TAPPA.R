#' Topological Analysis of Pathway Phenotype Association (TAPPA)
#'
#' The functions analyses the differential expression of pathways via TAPPA method. Expression is compared between two groups of samples by Mann-Whitney test. P-values are later adjusted for multiple hypothesis testing by Benjamini-Hochberg's FDR method.
#' @param x An \code{ExpressionSet} object or a gene expression data matrix or count matrix, rows refer to genes, columns to samples
#' @param group  Name or number of the phenoData column or a character vector or factor that contains required class assigments
#' @param pathways  A list of pathways in a form from \code{graphite} package or created by \code{preparePathways()}
#' @param type Type of the input data, \code{"MA"} for microarray and \code{"RNASeq"} for RNA-Seq
#' @param which Character, which type of nodes is preserved in a pathway. Possible values are \code{"proteins"},\code{"metabolites"},\code{"mixed"}
#' @param edgeType Character, which type of edges is preserved in a pathway. If \code{NULL}, all edges are kept.
#' @param preparePaths  Logical, by default the pathways are transformed with \code{preparePathways()}. Use \code{FALSE}, if you have done this transformation separately
#' @param norm.method Character, the method to normalize RNAseq data. If \code{NULL} then vst-normalization is performed. Possible values are: \code{"edgeR", "vst", "rLog", "none"}
#' @param test.method Character, the method for differentiall expression analysis of RNAseq data. If \code{NULL} then \code{"voomlimma"} is used. Possible values are: \code{"DESeq2", "voomlimma", "vstlimma", "edgeR"}. This analysis is needed only for the visualization.
#' @param test Function implementing a statistical test comparing PCI scores between groups. It is employed as \code{test(PCI~group)$p.value}, where \code{PCI} is a numeric vector of the same length as \code{group}
#' @param normalize  Logical, should data be normalized?
#' @param verbose Logical, if \code{TRUE} names of the pathways are printed as they are analysed
#' @param both.directions,maxNodes,minEdges,commonTh,filterSPIA,convertTo,convertBy Arguments for the \code{preparePathways()}
#'
#' @return A list,
#' \item{res}{A data frame, rows refer to pathways. Columns contain: number of valid PCI-scores, median, min and max of the PCI scores for each group of samples, p-value of the \code{test} (\code{p.val}) and adjusted p-value (\code{p.adj}). If less than two nodes are present in the data, the function puts \code{NA}'s in all columns. }
#' \item{topo.sig}{\code{NULL}, it is preserved for the compatibility with other methods implemented in this package}
#' \item{degtest}{A numeric vector of gene-level differential expression statistics}
#'
#' @references Gao, S. and Wang, X. (2007) TAPPA: topological analysis of pathway phenotype association. Bioinformatics, 23, pages 3100-3102
#' @author Ivana Ihnatova
#' @keywords htest
#'
#' @examples
#' if (require(breastCancerVDX)) {
#' data("vdx")
#' pathways<-pathways("hsapiens","biocarta")[1:10]
#' MAdata<-Biobase::exprs(vdx)[,1:10]
#' rownames(MAdata)<-Biobase::fData(vdx)[,"Gene.symbol"]
#' MAdata<-MAdata[!duplicated(rownames(MAdata)),]
#'
#' TAPPA(MAdata, Biobase::pData(vdx)[,"er"][1:10], pathways, type="MA", convertTo="SYMBOL")
#' }
#' @export
TAPPA<-function(x, group, pathways, type, which="proteins", edgeType=NULL, preparePaths=TRUE, norm.method=NULL, test.method=NULL, test=t.test, normalize=TRUE, verbose=FALSE, both.directions=TRUE,
                maxNodes=150, minEdges=0, commonTh=2, filterSPIA=FALSE, convertTo="none", convertBy=NULL){
  gedm<-.prepareData(x, group, type, method="TAPPA", norm.method)
  if (preparePaths) paths<-.preparePathways(pathways, method="TAPPA", both.directions, rownames(gedm[[1]]),
                                            which=which, edgeType=edgeType,
                                            maxNodes, minEdges, commonTh, filterSPIA, convertTo, convertBy ) else
                                              paths<-pathways

    res<-.tappa(gedm[[1]], gedm[[2]], paths, test, normalize, verbose)
    if (type=="MA") deg.table<-.testMA(gedm[[1]], gedm[[2]])
    if (is.null(test.method)) test.method<-"voomlimma"
    if (type=="RNASeq") deg.table<-.testRNAseq(x, group, test.method)

    out<-list(res=res, topo.sig=NULL, degtable=deg.table)
    class(out)<-c(class(out), "topResultE","topResult")
    return(out)
}

#x GEDM, rows=genes, col=samples
.normalizationTAPPA<-function(x){
if (is.null(dim(x))) stop("Only matrices can be normalized")
x<-apply(x, 2, function(y) (y-mean(y))/sd(y)) #colum to zero mean and same scope
x<- 1/(1+exp(-x))-0.5
return(x)
}

#x normalized GEDM, rows=genes col=samples
#A adjacency matrix
.calculatePCI<-function(x, A){
sgn<-function(x){ifelse(x>0, 1, ifelse(x<0, -1, 0))}
# je normalizovana matica?
# je A spravna matica - s 1 na diagonale

#
A.ind<-which(A==1, arr.ind=TRUE)
if (nrow(A.ind>0)) {
PCI<-apply(A.ind, 1, function(ind) {
i<-ind[1]
j<-ind[2]
PCI<-sgn(x[i,]+x[j,])*sqrt(abs(x[i,]))*A[i,j]*sqrt(abs(x[j,]))
})
PCI<-rowSums(PCI)
norm.PCI=PCI/nrow(A)
} else norm.PCI<-NA
return(norm.PCI)
}

.extractsubset<-function(x,A){
if (is.matrix(A)) {
valid.genes<-rownames(A) %in% rownames(x)
if (sum(valid.genes)==0) warning("Expression data do not match node labels")
return(list(x=x[rownames(A)[valid.genes],], A=A[valid.genes, valid.genes]))
}
if (any(is(A)=="graphNEL")) {
valid.genes<-nodes(A) %in% rownames(x)
if (sum(valid.genes)==0) warning("Expression data do not match node labels")
return(list(x=x[nodes(A)[valid.genes],], A=A))
}

}

.tappaSingle<-function(path, group,x, test,  verbose) {
if (verbose) message(path[[2]])
A<-path[[1]]
valid.data<-.extractsubset(x,A)
PCI<-.calculatePCI(valid.data$x,valid.data$A)
desc<-unlist(tapply(PCI, group, function(x) c(N=length(x[!is.na(x)]), summary(x, na.rm=TRUE))))
out<-c(desc, p.value=test(PCI~group)$p.value)
return(out)
}

.tappa<-function(x, group, pathways, test, normalize=TRUE, verbose=FALSE ){

if (normalize) x<-.normalizationTAPPA(x)


out<-.catchErr(pathways, function(p) .tappaSingle(p, group, x,test, verbose))

if (length(out[[1]])>0) {
out[[1]]<-data.frame(t(vapply(out[[1]], function(x) x, numeric(15))))
out[[1]]$q.value<-p.adjust(out[[1]]$p.value,"fdr")
}
return(out)
}

