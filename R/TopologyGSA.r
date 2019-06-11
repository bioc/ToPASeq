#' Gene set analysis exploiting the topology of a pathway (TopologyGSA)
#'
#' TopologyGSA method uses graphical models to test the differential expression of a pathway. It also highlights pathway components involved in the deregulation.
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
#' @param method Either \code{"var"} and \code{"mean"}. Determine the type of test used by topologyGSA.
#' @param nperm Numeric, number of permutations.
#' @param alpha Numeric, threshold for statistical significance of variance test. It influences the method for the mean test
#' @param testCliques Logical, if \code{TRUE}, then the test is also performed on individual cliques. It can be very computationally complex.
#' @param both.directions,maxNodes,minEdges,commonTh,filterSPIA,convertTo,convertBy Arguments for the \code{preparePathways()}
#'
#' @details
#' The method requires a Directed Acyclic Graph (DAG). Therefore if a pathway contain also undirected or bidirected edges and error is thrown.
#'
#' The user can further specify for the mean test:
#' \enumerate{
#'  \item \strong{perms} number of permutations of the test,
#'  \item \strong{paired} logical, if TRUE Hotelling test for paired samples is calculated and the test on the variances is not performed
#' }
#'
#' Or for the variance test:
#' \enumerate{
#' \item \strong{variance} logical, if TRUE the estimates of the covariance matrices are included in the result.
#' \item \strong{s1} First group covariance matrix estimation.
#' \item \strong{s2} Second group covariance matrix estimation.
#' }
#'
#' @return A list
#' \item{res}{a list with one entry for each successfully analyzed pathway }
#' \item{topo.sig}{if \code{testCliques=TRUE}, a list where each slot contains the pvalues and a list of cliques in one pathway. \code{NULL} otherwise}
#' \item{degtest}{A numeric vector of gene-level differential expression statistics}
#'
#' @references Massa MS, Chiogna M, Romualdi C. Gene set analysis exploiting the topology of a pathway. BMC System Biol. 2010 Sep 1;4:121.
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
#' TopologyGSA(MAdata, Biobase::pData(vdx)[,"er"][1:10], pathways, type="MA", convertTo="SYMBOL", nperm=10)
#' }
#' }
#' @export

TopologyGSA<-function(x, group, pathways, type, which="proteins", edgeType=NULL, preparePaths=TRUE, norm.method=NULL, test.method=NULL , method="mean", nperm=1000, alpha=0.05, testCliques=FALSE, 
                      both.directions=TRUE, maxNodes=150, minEdges=0, commonTh=2, filterSPIA=FALSE, convertTo="none", convertBy=NULL ){
  gedm<-.prepareData(x, group, type, method="TopologyGSA", norm.method)
  if (preparePaths) paths<-.preparePathways(pathways, method="TopologyGSA", both.directions,rownames(gedm[[1]]), which=which, edgeType=edgeType,maxNodes, minEdges, commonTh, filterSPIA, convertTo, convertBy ) else paths<-pathways
  res<-.topologyGSA(gedm[[1]], gedm[[2]], paths, method, alpha, nperm=nperm )
  if (testCliques) {
    message("Testing cliques:\n")
    topo.sig<-.testCliques(gedm[[1]], gedm[[2]], paths, method, nperm)
  } else topo.sig<-NULL
  if (type=="MA") deg.table<-.testMA(gedm[[1]], gedm[[2]])
  if (is.null(test.method)) test.method<-"voomlimma"
  if (type=="RNASeq") deg.table<-.testRNAseq(x, group, test.method)

  out<-list(res=res, topo.sig=topo.sig, degtable=deg.table)
  class(out)<-c(class(out), "topResultC","topResult")
  return(out)
}


.topologyGSA<-function(x, group, pathways, test="mean", alpha, both.directions, nperm ){

  perms<-.preparePerms(group=group, nperm=nperm, method="TopologyGSA")

  message("Analysing pathway:\n")
  out<-.catchErr(pathways, function(p) {
    cat(p[[2]],"\n")
    if (test=="mean") {
      pvar<-.pathwayVarianceTest(t(x), group, p[[1]], nperm, "never")
      .pathwayMeanTest(t(x), group, p[[1]], nperm, "never", pvar$pval<=alpha, paired=FALSE)
    } else if (test=="var")
      .pathwayVarianceTest(t(x), group, p[[1]], nperm, "never")
      
  })

  if (length(out[[1]])>0)  {res<-sapply(out[[1]], function(x) unlist(x[sapply(x, function(y) mode(y) %in% c("numeric", "logical"))]))
  #graphs<-sapply(out[[1]], function(x) (x[sapply(x, function(y) mode(y) %in% c("S4", "list"))]))
  out[[1]]<-data.frame(t(res))
  if (length(out[[1]]>0)) out[[1]]$q.value<-p.adjust(out[[1]][,"pval"],"fdr")
  }
  return(out)
}


.testCliques<-function(x, group, pathways, test, nperm){
  xx<-list()
  if (test=="mean") {
    xx<-lapply(pathways, function(p) {
      cat(p[[2]],"\n")
      g<-.nodes(p[[1]])
      g<-g[g %in% rownames(x)]
      res<-.cliqueMeanTest(x[g,], group, p[[1]], nperm, "never", TRUE)
      res<-list(p.value=res$p.value, cliques=res$cliques)
    })
  }
  if (test=="var") {
    xx<-lapply(pathways, function(p) {
      cat(p[[2]],"\n")
      res<-.cliqueVarianceTest(x, group, p[[1]],nperm, "never")
      res<-list(p.value=res$p.value, cliques=res$cliques)
    })


  }
  return(xx)
}
