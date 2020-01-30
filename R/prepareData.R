#############################
# Preparing Expression Data #
#############################

.prepareData<-function(x, group, type, method, norm.method=NULL, test.method=NULL, 
                       nperm=NULL, ncores=NULL, p.th=0.05, logFC.th=1.5){

  if (method %in% c("TAPPA", "DEGraph", "TopologyGSA", "clipper")) {
    if (type=="MA") {
      out<-.processMA(x, group)
    } else
      if (type=="RNASeq") {
        message("Selected methods use directly expression data. Normalization may not be optimal.")
        if (is.null(norm.method)) norm.method<-"vst"
        out<-.processRNAseq(x, group, norm.method)
      } else stop("This inpout type is not supported for the selected method")

  } else
    if (method=="PWEA") {
      if (type=="RNASeq"){
        if (!is.null(norm.method)) warning("Normalization method ignored")
        if (is.null(test.method)) test.method<-"voomlimma"
        out<-.testRNAseq(x, group, test.method)
        if (test.method=="DESeq2") perm<-.testRNAseq(x, group, "DESeq2perm", nperm) else {
          test<-get(paste("test", test.method, sep=""))
          perm<-.parDE(x, group,test , nperm, ncores)
        }
        out<-list(out, perm)
      }   else
        if (type=="MA"){
          out<-.processMA(x, group)
          obs<-.testMA(out[[1]], out[[2]])
          test<-.testMA
          message("Preparing permutations..\n")
          perm<-.parDE(out[[1]], out[[2]], test, nperm, ncores)
          out<-list(obs, perm)
        } else
          if (type=="DEtable") {
            if (length(x)!=2) stop("The input data must be a list of length two (observed and random statistics)")
            if (!all(colnames(x[[1]])==c("ID", "logFC","t","pval","padj"))) stop('Colnames of the table differ from c("ID", "logFC","t","pval"."padj")"')
            out<-x
          } else stop("This inpout type is not supported for the selected method")

    } else
      if (method %in% c("PRS","SPIA", "CePa")) {
        if (type=="RNASeq"){
          if (!is.null(norm.method)) warning("Normalization method ignored")
          if (is.null(test.method)) test.method<-"voomlimma"
          out<-.testRNAseq(x, group, test.method)
          out<-.applyTh(out, p.th, logFC.th)
        } else
          if (type=="MA"){
            out<-.processMA(x, group)
            out<-.testMA(out[[1]], out[[2]])
            out<-.applyTh(out, p.th, logFC.th)
            if (length(out[[1]])==0) stop("Found NO differentially expressed genes") else message(paste("Found", length(out[[1]]), "differentially expressed genes"))
          } else
            if (type=="DEtable") {
              out<-.applyTh(x, p.th, logFC.th)
              if (length(out[[1]])==0) stop("Found NO differentially expressed genes") else message(paste("Found", length(out[[1]]), "differentially expressed genes"))
            } else
              if (type=="DElist") {
                if (length(x)!=2) stop("The input data must be a list of length two (named statistics of differentially expressed genes and names of all genes)")
                .checkDEandAll(x[[1]], x[[2]])
                out<-x
              } else stop("This inpout type is not supported for the selected method")

      } else stop(paste("Method",method,"is not implemented"))
  return(out)
}

.processMA<-function(x, group){
  if (any(is(x) == "ExpressionSet")) {
    exprs<-exprs(x)
    
    if (length(group)>1) 
      stop("'group' must be of length one (number or name of the column of pData())")
    
    if (is.character(group) & ! any(group == colnames(pData(x)))) 
      stop( paste( group,"is not present in phenoData"))  else
        if (is.numeric(group)  & group > ncol(pData(x))) 
          stop(paste("phenoData contain only", ncol(pData(x)), "columns")) else
            if (is.numeric(group) | (is.character(group)) ) 
              group<-pData(x)[,group] else
                stop ("'group' must be either numeric or character")
  }

  if (!any(is(x)%in% c("data.frame", "matrix"))) 
    stop("Invalid 'x' type")

  if (length(group)!=ncol(x)) stop("length of 'group' does not match number of samples (columns)")
  if (nlevels(factor(group))!=2) 
    stop(paste("'group' must contain two distict values, found", levels(factor(group))))

  return(list(x, factor(group)))
}

.testMA<-function(x, group){
  if (any(is.na(group))) {
    message(paste("Removing", sum(is.na(group)), "samples with missign group value"))
    x<-x[,!is.na(group)]
    group<-group[!is.na(group)]
  }
  
  if (requireNamespace("limma")) {
    message("Using limma to test differential expression of genes")
    cat(levels(group)[1],"denoted as 0 \n", levels(group)[2],"denoted as 1\n", "Contrasts: ", levels(group)[2], "-", levels(group)[1],"\n"    )
    design <- model.matrix(~0+factor(group))
    colnames(design) <- c("V1","V2")
    x<-data.frame(x)
    fit <- limma::lmFit(x, design)
    contrast.matrix = limma::makeContrasts(contrasts="V2-V1", levels=design)
    fit2 = limma::contrasts.fit(fit, contrast.matrix)
    fit2 = limma::eBayes(fit2)
    results <- limma::decideTests(fit2)
    deg.table = limma::topTable(fit2, coef=1, adjust.method="BH", number=nrow(x), genelist=rownames(x), sort.by="none")
    deg.table<-data.frame(ID=rownames(deg.table), logFC=deg.table$logFC, t=deg.table$t,  pval=deg.table$P.Value, padj=deg.table$adj.P.Val)
    rownames(deg.table)<-deg.table$ID
  } else {
    message("Using t-test to test differential expression of genes")
    x<-as.matrix(x)
    deg.table<-apply(x,1, function(r) {
      out<-t.test(r~group)
      out<-c(logFC=unname(diff(out$estimate)), out$statistic, pval=out$p.value)
      return(out)
    } )
    deg.table<-data.frame(t(deg.table))
    deg.table$ID<-rownames(deg.table)
    deg.table$padj<-p.adjust(deg.table$pval,"BH")
    deg.table<-deg.table[,c("ID","logFC","t","pval","padj")]
  }
  return(deg.table)
}


.processRNAseq<-function(x, group,norm.method){
  if (! is.integer(x)) stop("Integer count data expected")
  x<-x[rowSums(x!=0)>0,]
  #normalization
  if (norm.method=="edgeR") {
    if (requireNamespace("edgeR")) {
      x <-edgeR::DGEList(counts=x)
      x <- edgeR::calcNormFactors(x, method = 'TMM')
      norm.matrix <- edgeR::cpm(x, log=TRUE)
  }} else
  if (norm.method=="vst") {
    if (requireNamespace("DESeq2")) {
      norm.matrix <- DESeq2::varianceStabilizingTransformation(x)
  }} else
  if (norm.method=="rLog") {
    if (requireNamespace("DESeq2")) {
      norm.matrix<-DESeq2::rlogTransformation(x)
      rownames(norm.matrix)<-rownames(x)
  }} else
  if (norm.method=="none") 
    norm.matrix<-x  else 
    stop("Invalid normalization method")
  return(list(x=norm.matrix, factor(group)))
}


.testRNAseq<-function(x, group, test.method, nperm=0){
  if (test.method %in% c("DESeq2", "voomlimma", "vstlimma", "edgeR")) {
    deg.table<-get(paste("test", test.method, sep=""))(x,group)
  } else
    if (test.method=="DESeq2perm") {
      deg.table<-testDESeq2perm(x, group, nperm)
    } else
      stop("Invalid method for differential expression analysis")

  return(deg.table)
}


.parDE<-function(x, group, test, nperm, ncores=NULL){
  if (is.null(ncores)) ncores<-parallel::detectCores()

  perm.test<-function(i, test, x, group){test(x, sample(group))}

  if (.Platform$OS.type=="unix") cl <- parallel::makeCluster(ncores, type="FORK") else
    cl <- parallel::makeCluster(rep("localhost",ncores))

  i<-1
  tmp<-parallel::parSapply(cl, seq_len(nperm), perm.test, test, x, group)

  parallel::stopCluster(cl)

  return(tmp)
}



.applyTh<-function(x, p.th, logFC.th){
  select<-rep(FALSE, nrow(x))
  if (!is.na(p.th)) 
    select<- x$pval < p.th 
  if (!is.na(logFC.th)) 
    select<-abs(x$logFC) > logFC.th 
  if (!is.na(p.th) & !is.na(logFC.th)) 
    select<-(x$pval < p.th) & (abs(x$logFC) > logFC.th) 
  
  select[is.na(select)]<-FALSE
  de<-setNames(x$logFC[select], x$ID[select])
  out<-list(de=de, all=as.character(x$ID))
  return(out)
}


.checkDEandAll<-function(de, all){
  IDsNotP <- setdiff(names(de),all)
  if (length(IDsNotP)/length(de) > 0.01) {
    stop("More than 1% of your differentially expressed genes have IDs are not present in the reference array!. Are you sure you use the right reference array?")
  }
  if (!length(IDsNotP) == 0) {
    cat("The following differentially expressed genes are missing from all vector...:\n")
    cat(paste(IDsNotP, collapse = ","))
  }
  if (any(is.na(names(de)))) stop("Found NA's in the differentially expressed genes")
  if (length(intersect(names(de), all)) != length(de)) {
    stop("de must be a vector of log2 fold changes. The names of de should be included in the all!")
  }
}
