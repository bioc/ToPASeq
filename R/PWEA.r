#' PathWay Enrichment Analysis (PWEA)
#'
#' The function runs PWEA method (please see References for the details) on gene expression data matrix, vector specifing to which group a sample belongs and a list of pathway graphs. Briefly, it is a weighted GSEA-like method. The weightes are based on the distance and Pearson's correlation between genes in a pathway.

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
#' @param tif A list of Topology Influence Factor's. One slot refers to one pathway. Use \code{prepareTIF()} to create it. It is required only if \code{type=="DEtable"}
#' @param alpha Numeric, a theshold value used during TIF calculation
#' @param nperm Numeric, number of permutations. Used only if \code{x \%in\% c("MA", "RNASeq")}
#' @param ncores Numeric, number of cores. Used only if \code{x \%in\% c("MA", "RNASeq")}. The permutations are calculated in parallel way
#' @param both.directions,maxNodes,minEdges,commonTh,filterSPIA,convertTo,convertBy Arguments for the \code{preparePathways()}
#'
#' @return A list:
#' \item{res }{A data frame, rows refer to pathways. It contains: Enrichment score for a pathway, p-value and p-value adjusted for multiple hypothesis testing by Benjamini-Hochberg's FDR method. \code{NA}'s if less than 2 nodes are present in the data}
#' \item{topo.sig }{A list, topology influence factors for the genes in individual pathways. \code{NULL} if less than 2 nodes are present in the data}
#' \item{degtest }{A named vector of statistics from testing the differential expression}
#'
#' @references Hung, JH., Whitfield, T. W., Yang, TH., Hu, Z., Weng, Z., DeLisi, Ch. (2010) Identification of functional modules that correlate with phenotypic difference: the influence of network topology, Genome Biology, 11:R23
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
#' PWEA(MAdata, Biobase::pData(vdx)[,"er"][1:10], pathways, type="MA", convertTo="SYMBOL", nperm=10)
#' }
#' }
#' @export
PWEA<-function(x, group, pathways, type, which="proteins", edgeType=NULL, preparePaths=TRUE, norm.method=NULL, test.method=NULL, tif=NULL, alpha=0.05, nperm=1000,  ncores=1,
               both.directions=TRUE, maxNodes=150, minEdges=0, commonTh=2, filterSPIA=FALSE, convertTo="none", convertBy=NULL){
  gedm<-.prepareData(x, group, type, method="PWEA", norm.method, test.method, nperm, ncores)
  if (preparePaths) paths<-.preparePathways(pathways, method="PWEA", both.directions, rownames(gedm[[1]]), which=which, edgeType=edgeType,maxNodes, minEdges, commonTh, filterSPIA, convertTo, convertBy ) else paths<-pathways
  if (is.null(tif) & type=="DEtable") stop("Argument 'tif' is missing. Calculate TIF with prepareTIF() first.")
  if (is.null(tif)) tif<-.prepareTIF(paths, x, alpha)
  res<-.pwea(gedm[[1]], gedm[[2]], tif, paths, alpha)
  deg.table<-gedm[[1]]
  topo.sig<-lapply(paths, function(p) {
    g<-nodes(p[[1]])
    g<-g[g %in% gedm[[1]]$ID]
    tif[[p[[2]]]][g]

    # exprs.valid<-extractsubset(x, p[[1]])$x
    # return(TIF(x, exprs.valid))
  }
  )
  out<-list(res=res, topo.sig, degtable=deg.table)
  class(out)<-c(class(out), "topResultW", "topResult")
  return(out)
}

#TIF
#path = graphNEL draha
#x - expresia genov z drahy, rows=genes
.TIF<-function(path, x, alpha){
d<- RBGL::johnson.all.pairs.sp(path)
pcc<-stats::cor(t(x), method="pearson", use="pairwise.complete.obs")
d<-d[rownames(pcc),rownames(pcc)]
f<-d/abs(pcc)
diag(f)<-NA
valid<- f <= -log(alpha)
f<-as.data.frame(t(f))
valid<-as.data.frame(t(valid))
tif<-unlist(Map(function(x,y) if (sum(y, na.rm=TRUE)>0) mean(x[y], na.rm=TRUE) else NA, f,valid))
tif<-exp(-tif)
tif[is.na(tif)]<-0

return(tif)
}


#tif list tif vsetkych genov zo vsetkych drah
# alpha - alpha z povodnej metody
# n - pocet genov mimo drahy
.notPTIF<-function(tif, n){
tif.all<-unlist(tif)
prop.pass<-sum(tif.all>1)/length(tif)

tif.pass<-tif.all[tif.all > 0]

input<-rbinom(n,size=1,prob=prop.pass)
input[input==1]<-rnorm(sum(input), mean=mean(tif.pass), sd=sd(tif.pass))
return(input)
}

#geneList gene level statistics (t) ORDERED
# tif of all genes (exponent)
# g genes in pathway (geneSet)
#  zbalika DOSE, fortify = vsetky hodnoty (pre graf)
.gseaScores<-function (geneList, geneSet, exponent = 1, fortify = FALSE)
{
    geneSet <- intersect(geneSet, names(geneList))
    N <- length(geneList)
    Nh <- length(geneSet)
    Phit <- Pmiss <- numeric(N)
    hits <- names(geneList) %in% geneSet
 if (length(exponent)>1) Phit[hits] <- abs(geneList[hits])^exponent[names(geneList[hits])] else
     Phit[hits] <- abs(geneList[hits])^exponent

    NR <- sum(Phit)
    Phit <- cumsum(Phit/NR)
    Pmiss[!hits] <- 1/(N - Nh)
    Pmiss <- cumsum(Pmiss)
    runningES <- Phit - Pmiss
    max.ES <- max(runningES)
    min.ES <- min(runningES)
    if (abs(max.ES) > abs(min.ES)) {
        ES <- max.ES
    }    else {
        ES <- min.ES
    }
    if (fortify == TRUE) {
        df <- data.frame(x = seq_along(runningES), runningScore = runningES,
            position = as.integer(hits))
        return(df)
    }
    return(ES)
}

.gseaScoresCols<-function (geneList, geneSet, exponent = 1) {
    geneSet <- intersect(geneSet, rownames(geneList))
    N <- nrow(geneList)
    Nh <- length(geneSet)

    Phit <- Pmiss <- matrix(0, nrow=N, ncol=ncol(geneList))

    hits <- rownames(geneList) %in% geneSet

     if (length(exponent)>1) Phit[hits,] <- abs(geneList[hits,])^exponent[names(geneList[hits,])] else
     Phit[hits,] <- abs(geneList[hits,])^exponent

    NR <- colSums(Phit)
    NR <- matrix(rep(NR, each=N), nrow=N)
    Phit <- apply(Phit/NR, 2, cumsum)
    Pmiss[!hits,] <- 1/(N - Nh)
    Pmiss <- apply(Pmiss, 2, cumsum)
    runningES <- Phit - Pmiss

    max.ES <- colMax(runningES)
    min.ES <- colMin(runningES)
    ES<-ifelse(abs(max.ES) > abs(min.ES), max.ES, min.ES)

    return(ES)
}

.PWEAscore<-function(tif, tifout, geneList, geneSet){

genPow<-function(x, pow){ return(sign(x)*abs(x)^(pow)) }

if (length(geneList) != length(tif)+length(tifout)) stop("Number of Topology Impact Factors differs from the number of gene-level statistics")

tif.complete<-vector("numeric", length(geneList))
names(tif.complete)<-names(geneList)

tif.complete[names(tif)]<-tif
tif.complete[names(tifout)]<-tifout

r<-genPow(abs(geneList),(1+tif.complete))
r<-sort(r, decreasing=TRUE)

sc<-.gseaScores(r, geneSet, 1+tif)
return(sc)
}

.preparePermsPWEA<-function(x, gr, nperm, test){
out<-replicate(nperm, test(x, sample(gr))$stat)
return(out)
}

.prepareTIF<-function(pathways, exprs, alpha){

all.genes<-rownames(exprs)

inP<-lapply(pathways, function(x) {
 exprs.valid<-.extractsubset(exprs, x[[1]])$x
 return(.TIF(x[[1]], exprs.valid, alpha))
 }
 )

outP<-lapply(pathways, function(p){
g<-sum(! all.genes %in% nodes(p[[1]]))
return(.notPTIF(inP, g))
})

return(Map(c, inP, outP))
}

.PWEASingle<-function(p, geneList, tif, perms){
obs<-.gseaScores(geneList, nodes(p), tif, fortify=FALSE)
rnd<-.gseaScoresCols(perms, nodes(p))
p<-sum(obs>= rnd)/ length(rnd)
return(c(ES=obs, p.value=p))
}

.parsePerms<-function(xx){
ids<-xx[,1]$ID
out<-sapply(seq_len(ncol(xx)), function(i) {
m<-match(ids, as.character(xx[,i]$ID))
xx[,i]$t[m]
})
rownames(out)<-ids
return(out)
}

.pwea<-function(obs, perms, tif, pathways, alpha){


geneList<-obs$t
names(geneList)<-obs$ID
perms<-.parsePerms(perms)
out<-.catchErr(pathways, function(p) .PWEASingle(p[[1]], geneList, tif[[p[[2]]]], perms))
if(length(out[[1]])>0){
out[[1]]<-data.frame(t(vapply(out[[1]], function(x) x, numeric(2))))
out[[1]]$q.value<-p.adjust(out[[1]]$p.value,"fdr")
}
return(out)
}
