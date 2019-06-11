#' @import graphite
#' @importFrom parallel detectCores makeCluster clusterExport parSapply stopCluster
#' @importFrom SummarizedExperiment assay colData colData<-
#' @importClassesFrom Biobase ExpressionSet
#' @importMethodsFrom Biobase exprs pData fData
#' @importFrom Rcpp sourceCpp
#' @importFrom methods as is new slot slot<-
#' @importFrom stats sd setNames model.matrix p.adjust t.test cov pchisq pf qchisq density median na.omit phyper pnorm qnorm cor rbinom rnorm
#' @importFrom utils head
#' @importFrom graphics abline legend par points text
#' @importFrom gRbase is.DAG getCliques
#' @importFrom graph subGraph validGraph degree
#' @importFrom RBGL johnson.all.pairs.sp connectedComp floyd.warshall.all.pairs.sp
#' @importFrom limma voom lmFit eBayes topTable
#' @importFrom corpcor cov.shrink
#' @useDynLib ToPASeq
NULL
