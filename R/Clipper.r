#' clipper
#'
#' clipper is a method for topological gene set analysis. It implements a two-step empirical approach based on the exploitation of graph decomposition into a junction tree to reconstruct the most relevant signal path. In the first step clipper selects significant pathways according to statistical tests on the means and the concentration matrices of the graphs derived from pathway topologies. Then, it "clips" the whole pathway identifying the signal paths having the greatest association with a specific phenotype.
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
#' @param method Character, \code{"mean"}  or \code{"var"}, the kind of test to perform on the cliques
#' @param testCliques Logical, if \code{TRUE} then the test is applied also on the cliques of the each pathway. It is a very time consuming calculation, especially for many or big pathways
#' @param nperm Number of permutations, if \code{0} then asymptotic distribution is used. May not be valid when shrinked estimator is used. 
#' @param alphaV Numeric, the threshold for variance test. The calculation of mean test depends on the result of variance test.
#' @param both.directions,maxNodes,minEdges,commonTh,filterSPIA,convertTo,convertBy Arguments for the \code{preparePathways()}
#'
#' @return A list:
#' \item{res}{A list. First slot is a data frame containing p-values and q-values of mean and variance tests on pathways. The second slot is a list containing data.frames of the most affected paths in each pathway. The columns of the data frames contain: 1 - Index of the starting clique 2 - Index of the ending clique 3 - Index of the clique where the maximum value is reached 4 - length of the path 5 - maximum score of the path 6 - average score along the path 7 - percentage of path activation 8 - impact of the path on the entire pathway 9 - clique involved and significant 10 - clique forming the path 11 - genes forming the significant cliques 12 - genes forming the path}
#' \item{topo.sig}{if \code{testCliques=TRUE}, a list where each slot contains the pvalues and a list of cliques in one pathway. \code{NULL} otherwise}
#' \item{degtest}{A data.frame of gene-level differential expression statistics}
#'
#' @references Martini P, Sales G, Massa MS, Chiogna M, Romualdi C. Along signal paths: an empirical gene set approach exploiting pathway topology. Nucleic Acids Res. 2013 Jan 7;41(1):e19. doi: 10.1093/nar/gks866. Epub 2012 Sep 21. PubMed PMID: 23002139; PubMed Central PMCID: PMC3592432.
#'
#' @author Ivana Ihnatova
#' @keywords htest
#'
#' @examples
#' \dontrun{
#' if (require(breastCancerVDX)) {
#' data("vdx")
#' pathways<-pathways("hsapiens","biocarta")[1:3]
#' MAdata<-Biobase::exprs(vdx)[,1:10]
#' rownames(MAdata)<-Biobase::fData(vdx)[,"Gene.symbol"]
#' MAdata<-MAdata[!duplicated(rownames(MAdata)),]
#'
#' clipper(MAdata, Biobase::pData(vdx)[,"er"][1:10], pathways, type="MA", convertTo="SYMBOL", nperm=10)
#' }
#' }
#' @export
clipper<-function(x, group, pathways, type,which="proteins", edgeType=NULL, preparePaths=TRUE, norm.method=NULL, test.method=NULL,
                  method="mean", testCliques=FALSE, nperm=1000, alphaV=0.05,
                  both.directions=TRUE, maxNodes=150, minEdges=0, commonTh=2, filterSPIA=FALSE, convertTo="none", convertBy=NULL){
  gedm<-.prepareData(x, group, type, method="clipper", norm.method)
  if (preparePaths) paths<-.preparePathways(pathways, method="clipper", both.directions, rownames(gedm[[1]]), which=which, edgeType=edgeType,maxNodes, minEdges, commonTh, filterSPIA, convertTo, convertBy ) else paths<-pathways
  res<-.CLIPPER(paths, gedm[[1]], gedm[[2]], method, testCliques, nperm, alphaV)
  if (type=="MA") deg.table<-.testMA(gedm[[1]], gedm[[2]])
  if (is.null(test.method)) test.method<-"voomlimma"
  if (type=="RNASeq") deg.table<-.testRNAseq(x, group, test.method)

  out<-list(res=res$results[1:2], topo.sig=res$results[[3]], degtable=deg.table)
  class(out)<-c(class(out), "topResultC","topResult")
  return(out)
}




.CLIPPER<-function (pathways, expr, classes, method, testCliques=FALSE, nperms, alphaV=0.05){
    #if (!require(clipper))        stop("library clipper is missing")

    classes<-factor(as.numeric(factor(classes)))
message("Analysing pathway:\n")
      out<-.catchErr(pathways, function(p) {
      cat(p[[2]],"\n")
        expr<-expr[rownames(expr) %in% nodes(p[[1]]),]  
       if (method=="mean") {
         pvar<-.pathwayVarianceTest(t(expr), classes, p[[1]], nperms, "decide")
         .pathwayMeanTest(t(expr), classes, p[[1]], nperms, "decide", pvar$pval<=alphaV, paired=FALSE)
       } else if (method=="var")
         .pathwayVarianceTest(t(expr), classes, p[[1]], nperms, "decide")
         
})
cliq.test<-list()
if (testCliques) {
message("Testing cliques\n")
 cliq.test<-lapply(pathways, function(p) {
 cat(p[[2]],"\n")
 if (method=="mean") cliq<-.cliqueMeanTest(t(expr), classes, p[[1]], nperms, "decide", TRUE)
 if (method=="var") cliq<-.cliqueVarianceTest(t(expr), classes, p[[1]], nperms, "decide")

 return(cliq)
 })

  }
if (length(out[[1]])>0) {
res<-data.frame(t(sapply(out[[1]],function(x) x[[1]])))
res$q.value<-p.adjust(res$pval,"fdr")

paths<-lapply(out[[1]],function(x) x[[2]])

out[[1]]<-list(results=res, paths=paths, cliq.test)
} else
out[[1]]<-list(results=out[1:2], paths=list(), cliq.test)
return(out)
}



