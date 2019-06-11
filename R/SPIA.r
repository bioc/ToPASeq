#' Signaling Pathway Impact Analysis (SPIA)
#'
#' The function runs SPIA method on microarray or RNA-Seq data. The implementatio includes the identification of differentially expressed genes and transformation of pathways' topologies to an appropriate form. The SPIA method combines two independent p-values. One p-value comes from overrepresentation analysis and the other is so called pertubation factor.
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
#' @param combine Character, the method to combine p-values. Defaults to \code{"fisher"} for Fisher's method. The other possible value is \code{"norminv"} for the normal inversion method.
#' @param both.directions,maxNodes,minEdges,commonTh,filterSPIA,convertTo,convertBy Arguments for the \code{preparePathways()}
#'
#' @return A list:
#' \item{res}{A matrix with columns as descibed below:
#'    pSize - Pathway size, number of genes,
#'    NDE - Number of differentially expressed genes,
#'    pNDE - P-value of the overrepresentation part of the method,
#'    tA - The observed total preturbation accumulation in the pathway,
#'    pPERT - P-value of the pertubation part of the method,
#'    p - Combined p-value (overrepresentation and pertubation),
#'    pFdr - False discovery rate adjusted \code{p},
#'    pFWER - FWER adjusted \code{p},
#'    Status - If a pathway was identified as Acivated or Inhibited
#'    }
#' \item{topo.sig}{A list of accumulated pertubation factors and log fold-changes for genes in individual pathways}
#' \item{degtest}{A numeric vector of gene-level differential expression statistics of all genes in the dataset}
#'
#' @references Tarca AL, Draghici S, Khatri P, Hassan SS, Mittal P, Kim JS, Kim CJ, Kusanovic JP, Romero R. A novel signaling pathway impact analysis. Bioinformatics. 2009 Jan 1;25(1):75-82.
#'
#' Adi L. Tarca, Sorin Draghici, Purvesh Khatri, et. al, A Signaling Pathway Impact Analysis for Microarray Experiments, 2008, Bioinformatics, 2009, 25(1):75-82.
#'
#' Draghici, S., Khatri, P., Tarca, A.L., Amin, K., Done, A., Voichita, C., Georgescu, C., Romero, R.: A systems biology approach for pathway level analysis. Genome Research, 17, 2007. Massa MS, Chiogna M, Romualdi C. Gene set analysis exploiting the topology of a pathway. BMC System Biol. 2010 Sep 1;4:121.
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
#' SPIA(MAdata, Biobase::pData(vdx)[,"er"][1:10], pathways, type="MA", convertTo="SYMBOL", logFC.th=-1)
#' }
#' @export
SPIA<-function(x, group, pathways, type,  which="proteins", edgeType=NULL,preparePaths=TRUE, norm.method=NULL, test.method=NULL, p.th=0.05, logFC.th=2, nperm=1000, combine="fisher",
               both.directions=TRUE, maxNodes=150, minEdges=0, commonTh=2, filterSPIA=FALSE, convertTo="none", convertBy=NULL){
  degs<-.prepareData(x, group, type, method="SPIA", norm.method, test.method, p.th=p.th, logFC.th=logFC.th)

  if (preparePaths) paths<-.preparePathways(pathways, method="SPIA", both.directions, degs[[2]],  which=which, edgeType=edgeType,maxNodes, minEdges, commonTh, filterSPIA, convertTo, convertBy ) else paths<-pathways
  res<-.spia(degs[[1]], degs[[2]], paths, nperm, combine)

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
  topo.sig<-.collectWeightsSPIA(degs[[1]], degs[[2]], paths)
  out<-list(res=res, topo.sig=topo.sig, degtable=deg.table)
  class(out)<-c(class(out), "topResultW","topResult")
  return(out)
}


.combfunc<-function (p1 = NULL, p2 = NULL, combine = "fisher") {
    tm = stats::na.omit(c(p1, p2))
    if (!all(tm >= 0 & tm <= 1)) {
        stop("values of p1 and p2 have to be >=0 and <=1 or NAs")
    }
    if (combine == "fisher") {
        k = p1 * p2
        comb = k - k * log(k)
        comb[is.na(p1)] <- p2[is.na(p1)]
        comb[is.na(p2)] <- p1[is.na(p2)]
        return(comb)
    }
    if (combine == "norminv") {
        comb = stats::pnorm((stats::qnorm(p1) + stats::qnorm(p2))/sqrt(2))
        comb[is.na(p1)] <- p2[is.na(p1)]
        comb[is.na(p2)] <- p1[is.na(p2)]
        return(comb)
    }
}

.SPIASingle<-function(de, all, M, nB, combine="fisher"){

 ph <- pb <- pcomb <- nGP <- pSize <- smPFS <- tA <- tAraw <- NULL

diag(M) <- diag(M) - 1
X <- de[rownames(M)]
noMy <- sum(!is.na(X))
okg <- intersect(rownames(M), all)
ok <- rownames(M) %in% all

if (!((noMy) > 0 & (abs(det(M)) > 1e-07))) {
  pb <- ph <- smPFS <- pcomb <- tAraw <- tA <- ob <- NA
  pfs<- rep(NA, length(X)) } else
{
X[is.na(X)] <- 0
pfs <- solve(M, -X)

pfstmp<-replicate(nB, {
x <- rep(0, length(X))
names(x) <- rownames(M)
x[ok][sample(1:sum(ok), noMy)] <- as.vector(sample(de,noMy))
tt <- solve(M, -x)
sum(tt-x)
})

mnn <- stats::median(pfstmp)
pfstmp <- pfstmp - mnn
ob <- sum(pfs - X) - mnn

if (ob > 0) pb <- sum(pfstmp >= ob)/length(pfstmp) * 2
if (ob < 0) pb <- sum(pfstmp <= ob)/length(pfstmp) * 2
if (ob == 0) if (all(pfstmp == 0)) pb <- NA else pb <- 1

if (!is.na(pb)) {
if (pb <= 0)  pb <- 1/nB/100
if (pb > 1)   pb <- 1
                }

}

nGP   <- noMy
pSize <- length(okg)
smPFS <- round(sum(pfs - X),3)
tAraw <- round(smPFS,3)
ph    <- round(stats::phyper(q = noMy - 1, m = length(okg), n = length(all) - length(okg), k = length(de), lower.tail = FALSE),3)
tA    <- round(ob,3)
pcomb <- round(.combfunc(pb, ph, combine),3)

return(c(NDE=nGP, pSize=pSize, smPFS=smPFS, tAraw=tAraw, pNDE=ph, tA=tA, pG=pcomb, pPERT=pb))
}

.spia<-function(de, all, pathways, perm, combine){

.checkDEandAll(de, all)

out<-.catchErr(pathways, function(p) .SPIASingle(de, all, p[[1]], perm, combine))
out[[1]]<-data.frame(t(sapply(out[[1]], function(x) x)))

if (length(out[[1]])>0) {
tmp<-out[[1]]

    tmp$pGFdr = p.adjust(tmp$pG, "fdr")
    tmp$pGFWER = p.adjust(tmp$pG, "bonferroni")
    tmp$Status = ifelse(tmp$tA > 0,"Activated", "Inhibited"  )

tmp<-tmp[,c("pSize", "NDE", "pNDE", "tA", "pPERT", "pG", "pGFdr", "pGFWER", "Status")]
out[[1]]<-tmp
}
return(out)
}

.SPIAweights<-function(de, all, M){

diag(M) <- diag(M) - 1
X <- de[rownames(M)]
noMy <- sum(!is.na(X))
okg <- intersect(rownames(M), all)
ok <- rownames(M) %in% all

if (!((noMy) > 0 & (abs(det(M)) > 1e-07))) {
 pfs<-setNames(rep(NA, length(X)), rownames(M))

 } else {
X[is.na(X)] <- 0
pfs <- solve(M, -X)
}
return(pfs-X)
}

.collectWeightsSPIA<-function(de, all, pathways){
out<-.catchErr(pathways, function(p) .SPIAweights(de, all, p[[1]]))
return(out[[1]])
}

.SPIA4plot<-function(de,all,M, nB, name, combine="fisher"){


 ph <- pb <- tA <-  NULL

diag(M) <- diag(M) - 1
X <- de[rownames(M)]
noMy <- sum(!is.na(X))
okg <- intersect(rownames(M), all)
ok <- rownames(M) %in% all

if (!(noMy) > 0 & (abs(det(M)) > 1e-07)) {
  pb <- ph  <- tA <- NA } else
{
X[is.na(X)] <- 0
pfs <- solve(M, -X)

pfstmp<-replicate(nB, {
x <- rep(0, length(X))
names(x) <- rownames(M)
x[ok][sample(1:sum(ok), noMy)] <- as.vector(sample(de,noMy))
tt <- solve(M, -x)
sum(tt-x)
})


mnn <- stats::median(pfstmp)
pfstmp <- pfstmp - mnn
ob <- sum(pfs - X) - mnn


if (ob > 0) pb <- sum(pfstmp >= ob)/length(pfstmp) * 2
if (ob < 0) pb <- sum(pfstmp <= ob)/length(pfstmp) * 2

if (pb <= 0)  pb <- 1/nB/100
if (pb > 1)   pb <- 1

if (ob == 0 & all(pfstmp == 0)) pb <- NA  else pb <- 1

}

tA    <- ob
return(list(log2FCH=X, Pertubation=pfs, NetPertubation=pfs-X, RandomPertubation=pfstmp, pPERT=pb, ObservedTotalPERT=tA))}

.drawSPIA<-function(x, name){
X<-x$log2FCH
pfs<-x$Pertubation
pfstmp<-x$RandomPertubation
pb<-x$pPERT
tA<-x$ObservedTotalPERT
xr<-range(X)
yr<-range(pfs-X)
par(mfrow = c(1, 2))
plot(X, pfs - X, xlim=xr, ylim=yr, main = paste("pathway ID=", name, sep = ""), xlab = "Log2 FC", ylab = "Perturbation accumulation (Acc)", cex.main = 0.8, cex.lab = 1.2)
abline(h = 0, lwd = 2, col = "darkgrey")
abline(v = 0, lwd = 2, col = "darkgrey")
points(X[abs(X) > 0 & X == pfs], pfs[abs(X) > 0 & X == pfs] - X[abs(X) > 0 & X == pfs], col = "blue", pch = 19, cex = 1.4)
points(X[abs(X) > 0 & X != pfs], pfs[abs(X) > 0 & X != pfs] - X[abs(X) > 0 & X != pfs], col = "red", pch = 19, cex = 1.4)
points(X[abs(X) == 0 & X == pfs], pfs[abs(X) == 0 & X == pfs] - X[abs(X) == 0 & X == pfs], col = "black", pch = 19, cex = 1.4)
points(X[abs(X) == 0 & X != pfs], pfs[abs(X) == 0 & X != pfs] - X[abs(X) == 0 & X != pfs], col = "green", pch = 19, cex = 1.4)
legend(x="topright", legend=c("DE, not pertubed","DE, pertubed","nonDE, not pertubed","nonDE, pertubed"), col=c("blue", "red", "black","green"), pch=19)
if (requireNamespace("plotrix")) plotrix::thigmophobe.labels(X, pfs-X, labels=names(X)) else text(X, pfs-X, labels=names(X), pos=3)
plot(density(pfstmp, bw = sd(pfstmp)/4), cex.lab = 1.2, col = "black", lwd = 2, main = paste("pathway ID=", name, "  P PERT=",
   round(pb, 5), sep = ""), xlim = c(min(c(tA - 0.5, pfstmp)), max(c(tA + 0.5, pfstmp))), cex.main = 0.8,
   xlab = "Total Perturbation Accumulation (TA)")
abline(v = 0, col = "grey", lwd = 2)
abline(v = tA, col = "red", lwd = 3)

}




