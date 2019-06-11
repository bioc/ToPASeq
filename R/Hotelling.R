# variance estimation
.estimateCov<-function(x1, x2=NULL, type, equal) { #clipper estimateCov
  if (type=="scale") { #toplogyGSA .estimateCov()
    s1<-.estimateScaledCov(x1)
    if (is.null(x2)) s2<-NULL else s2<-.estimateScaledCov(x2)
  } else if (type=="shrink") { #clipper estimateExprCov(x, TRUE)
    s1<-.estimateShrinkedCov(x1)
    if (is.null(x2)) s2<-NULL else s2<-.estimateShrinkedCov(x2)
  } else if (type=="cov") {#clipper estimateExprCov(x, FALSE)
    s1<-cov(x1)
    if (is.null(x2)) s2<-NULL else s2<-cov(x2)
  }
  cov <- list(s1 = s1, s2 = s2)
  if (is.null(x2)) cov$s<-s1 else {
    n1 <- NROW(x1)
    n2 <- NROW(x2)
    if (equal) cov$s <- (s1 * (n1 - 1) + s2 * (n2 - 1))/(n1 + n2 - 2) else  cov$s <- (s1 / n1) + (s2 / n2) 
  }
  return(cov) 
}

.estimateScaledCov<-function(x){
  x.mean <- colMeans(x)
  names(x.mean) <- NULL
  x.scal <- x - matrix(rep(x.mean, NROW(x) * NCOL(x)), NROW(x), NCOL(x), byrow = TRUE)
  s <- cov(x.scal)
  return(s)
}

.estimateShrinkedCov<-function(x){ #estimateExprCov(x, TRUE)
  s<-unclass(corpcor::cov.shrink(x, verbose = FALSE))
  return(s)
}

.graphCov<-function(covList, cliques){
  covList <- lapply(covList, function(x) if (is.null(x)) NULL else qpIPF(x, cliques))
  return(covList)
}

.subsetCov<-function(covList, cli) {
  lapply(covList, function(x) x[cli,cli])
}

# mean test
# x  expression data matrix, rows samples, columns genes
# group numeric or factor with two values or levels. Assigns samples to two groups
# nperm number of permutations, if 0 then chi-square distribution is used
# 
.hotellingIPF<-function(x, group, cliques=NULL, nperm, type, equal, paired) {
  group<-as.numeric(factor(group))
  x1<-x[group==1,, drop=FALSE]
  x2<-x[group==2,, drop=FALSE]
  if (paired) {
    cov<-.estimateCov(x1-x2, type=type, equal=equal)
  } else {
    cov<-.estimateCov(x1,x2, type, equal)
  }
  if (!is.null(cliques)) cov<-.graphCov(cov,cliques)

  if (nperm == 0) {
    if (paired) {
      t.obs<-.hote_mean_paired(x1,x2, cov$s)
      pval <- 1 - pf(t.obs$t.obs, ncol(x), nrow(x1) - ncol(x))
    } else {
      t.obs<- .hote_mean(x1,x2, cov$s)
      pval <- 1 - pf(t.obs$t.obs, ncol(x), nrow(x1) + nrow(x2) - ncol(x) - 1)
    }
  } else {
    if (paired) {
      t.obs<-.hote_mean_paired(x1,x2, cov$s)
      stat.perm<-vapply(seq_len(nperm), function(i) {
        x2perm<-x2[sample(NROW(x2)),, drop=FALSE]
        cov<-.estimateCov(x1-x2perm, type=type, equal=equal)
        cov<-.graphCov(cov,cliques)
        .hote_mean_paired(x1,x2perm, cov$s)
      }, numeric(1))
    } else {
      t.obs<- .hote_mean(x1,x2, cov$s)
      stat.perm<-vapply(seq_len(nperm), function(i) {
        rgroup <- sample(group)
        x1perm <- x[rgroup==1, , drop=FALSE]
        x2perm <- x[rgroup==2, , drop=FALSE]
        cov<-.estimateCov(x1perm, x2perm, type, equal)
        cov<-.graphCov(cov,cliques)
        .hote_mean(x1perm, x2perm, cov$s )  
      }, numeric(1))
    }
    pval <- sum(stat.perm >= t.obs)/nperm
  }
  out<-list(pval = pval, t.obs = t.obs)
}

.hote_mean<-function(x1,x2, s) {
  n1 <- nrow(x1)
  n2 <- nrow(x2)
  ngene <- ncol(x1)
  m1 <- colMeans(x1)
  m2 <- colMeans(x2)
  mdiff <- m1 - m2
  
  t2 <- tryCatch(((n1 * n2)/(n1 + n2)) * (mdiff %*% solve(s) %*% mdiff), error = function(e) return(NA)) #hottelling's two sample t-squared statistic
  k <- n1 + n2 - ngene - 1  
  t.obs <- as.vector(t2 * k/(ngene * (n1 + n2 - 2))) #correction for approximation to F-distribution
  pval <- 1 - pf(t.obs, ngene, n1 + n2 - ngene - 1)
  return(list(t.obs=t.obs, pval=pval))
}

.hote_mean_paired<-function(x1,x2, s) {
  n1 <- nrow(x1)
  xdiff <- x1 - x2
  mx <- colMeans(xdiff)
  
  t2 <- as.numeric(n1 * (t(mx) %*% solve(s) %*% mx))
  ngene <- ncol(x1)
  k <- n1 - ngene
  t2<- t2 * k/(ngene * (n1 - 1))
  pval <- 1 - pf(t2, ngene, n1 - ngene)
  return(list(t.obs=t2, pval=pval))
}



# variance

.varianceTest<-function(covList, n1, n2, ngene, nedge){
  
  if (is.null(covList$s2)) stop("One-sample or paired test not supported.") 
  
  lambda.value<-tryCatch({
    k1.hat <- solve(covList$s1)
    k2.hat <- solve(covList$s2)
    k.hat <- solve(covList$s)
    k1.det <- det(k1.hat)
    k2.det <- det(k2.hat)
    k.det <- det(k.hat)
    n1 * log(k1.det/k.det) + n2 * log(k1.det/k.det)  
  }, error = function(e) return(NA))
  
  df <- nedge + ngene
  #qchisq.value <- qchisq(1 - alpha, df)
  pval <- 1 - pchisq(lambda.value, df)
  
  out<-list(lambda=lambda.value, pval=pval, df=df)
  return(out)
}

#check data
.finalizeData <-function(x, group, graph, paired, shrink){
  group<-as.numeric(factor(group))
  x1<-x[group==1,, drop=FALSE]
  x2<-x[group==2,, drop=FALSE]
  
  if (paired) {
    if (nrow(x1) != nrow(x2)) {
      stop("Your are working woth paired mode. The number of samples per class must be equal (and paired).")
    }
  }
  common <- intersect(colnames(x1), nodes(graph))
  #cat(paste0("Expression of ", length(common), " genes (out of ",length(nodes(graph)),") is available\n"))
  x1 <- x1[, common, drop = FALSE]
  x2 <- x2[, common, drop = FALSE]
  
  n1 <- nrow(x1)
  n2 <- nrow(x2)
  ngene <- ncol(x1)
  
  dag <- graph::subGraph(common, graph)
  
  graph <- .cli(dag)
  
  maxCliqueSize <- max(sapply(graph$cliques, length))
  
  type<-.chooseVariance(shrink, n1,n2, maxCliqueSize)
  list(x1 = x1, x2 = x2, n1=n1, n2=n2, ngene=ngene, type=type, graph = graph)
    
}

.cli<-function (dag){
  if (sum(diag(as(dag, "matrix"))) != 0) {
    dag <- .removeSelfLoops(dag)
  }
  moral <- .moralizeMG(dag)
  cli<-gRbase::getCliques(moral)
  cli<-lapply(cli, function(x) sort(match(x, nodes(dag))))
  list( cliques = cli,  graph = moral)
}

.moralizeMG<-function (graph){
  m <- gRbase::graphNEL2M(graph)
  m <- gRbase::moralizeMAT(m)
  gRbase::coerceGraph(m, "graphNEL")
}

.removeSelfLoops<-function (graph){
  edgeL <- graph@edgeL
  for (i in 1:length(edgeL)) {
    pos <- match(i, edgeL[[i]]$edges)
    if (!(is.na(pos))) 
      edgeL[[i]]$edges <- edgeL[[i]]$edges[-pos]
  }
  graph@edgeL <- edgeL
  return(graph)
}

.chooseVariance<-function(shrink, n1, n2, maxClique){
  if (shrink=="always") {
    type <- "shrink"
  } else {
    check <- n1 <= maxClique || n2 <= maxClique 
    if (check) {
      if (shrink=="never") stop("Both groups should have more than ", maxClique, " rows (samples)") 
      if (shrink=="decide") {
        type <- "shrink"
      }
    } else {
      if (shrink=="never") {
        type <-"scale"
      }
      if (shrink=="decide") {
        type <-"cov"
      }
    }
  }
  return(type)
}

.pathwayVarianceTest<-function(x, group, graph, nperm, shrink) {
  D<-.finalizeData(x,group, graph, FALSE, shrink)
  x1<-D$x1
  x2<-D$x2
  n1<-D$n1
  n2<-D$n2
  ngene<-D$ngene
  type<-D$type
  graph<-D$graph$graph
  cliques<-D$graph$cliques
  nedge<-(sum(as(graph, "matrix"))/2) 
 
  
  maxCliqueSize <- max(sapply(cliques, length))
  type<-.chooseVariance(shrink, NROW(x1), NROW(x2), maxCliqueSize)
  cov <-.estimateCov(x1, x2, type, TRUE)
  cov<-.graphCov(cov, cliques)
  obs<-.varianceTest(cov, n1, n2, ngene, nedge)
  
  if (type=="shrink" || nperm!=0) {
    if (type=="shrink" & nperm==0) stop("Permutation-based test must be performed. Please set the number of permutations")
    #cat("Performing permutations..")
    stat.perm<-vapply(seq_len(nperm), function(i) {
      rgroup <- sample(group)
      x1perm <- x[rgroup==1, , drop=FALSE]
      x2perm <- x[rgroup==2, , drop=FALSE]
      cov<-.estimateCov(x1perm, x2perm, type, TRUE)
      cov<-.graphCov(cov, cliques)
      .varianceTest(cov, n1, n2, ngene, nedge )$lambda  
      
    }, numeric(1))
    pval <- mean(stat.perm >= obs$lambda, na.rm=TRUE) #some stat.perm can be NA
    obs$pval<-pval
    obs$df<-NULL
    
    
  }
  obs$type<-type
  return(obs)
}

# equal TRUE ak pathwayVarianceTest nevyznamny, FALSE ak vyznamny
.pathwayMeanTest<-function(x, group, graph, nperm, shrink, equal, paired){
  D<-.finalizeData(x,group, graph, FALSE)
  
  x1<-D$x1
  x2<-D$x2
  n1 <- nrow(x1)
  n2 <- nrow(x2)
  ngene <- ncol(x1)
  
  graph<-D$graph$graph
  nedge<-(sum(as(graph, "matrix"))/2) 
  cliques<-D$graph$cliques
  
  maxCliqueSize <- max(sapply(cliques, length))
  type<-.chooseVariance(shrink, n1, n2, maxCliqueSize)
  
  if (type=="shrink" & nperm==0) stop("Permutation-based test must be performed. Please set the number of permutations")
  obs<-.hotellingIPF(x, group, cliques, nperm, type, equal, paired)
  obs$type<-type
  return(obs)
  
}

.cliqueVarianceTest<-function(x, group, graph, nperm, shrink){
  D<-.finalizeData(x,group, graph, FALSE)
  
  x1<-D$x1
  x2<-D$x2
  n1 <- nrow(x1)
  n2 <- nrow(x2)
  
  graph<-D$graph$graph
  cliques<-D$graph$cliques
  
  maxCliqueSize <- max(sapply(cliques, length))
  type<-.chooseVariance(shrink, n1, n2, maxCliqueSize)
  cov <- .estimateCov(x1, x2, type, TRUE) 
  
  sapply(cliques, function(cli) {
    nedge <- length(cli)*(length(cli)-1)/2
    cov<-.subsetCov(cov, cli)
    obs<-.varianceTest(cov, n1, n2, length(cli), nedge)
    if (type=="shrink" || nperm!=0) {
      if (type=="shrink" & nperm==0) stop("Permutation-based test must be performed. Please set the number of permutations")
      #cat("Performing permutations..")
      stat.perm<-vapply(seq_len(nperm), function(i) {
        rgroup <- sample(group)
        x1perm <- x[rgroup==1, , drop=FALSE]
        x2perm <- x[rgroup==2, , drop=FALSE]
        cov<-.estimateCov(x1perm, x2perm, type, TRUE)
        cov<-.subsetCov(cov, cli)
        .varianceTest(cov, n1, n2, length(cli), nedge )$lambda  
      }, numeric(1))
      pval <- mean(stat.perm >= obs$lambda, na.rm=TRUE) #some stat.perm can be NA
      obs$pval<-pval
      obs$df<-NULL
    }
    return(obs)   
  })
  
}

.cliqueMeanTest<-function(x, group, graph, nperm, shrink, equal, paired){
  D<-.finalizeData(x,group, graph, FALSE)
  
  x1<-D$x1
  x2<-D$x2
  n1 <- nrow(x1)
  n2 <- nrow(x2)
  ngene <- ncol(x1)
  
  graph<-D$graph$graph
  nedge<-(sum(as(graph, "matrix"))/2) 
  cliques<-D$graph$cliques
  
  maxCliqueSize <- max(sapply(cliques, length))
  type<-.chooseVariance(shrink, NROW(x1), NROW(x2), maxCliqueSize)

  
  if (type=="shrink" & nperm==0) stop("Permutation-based test must be performed. Please set the number of permutations")
  obs<-.hotellingIPF(x, group, cliques, nperm, type, equal, paired)
  
  return(obs)
  
}
