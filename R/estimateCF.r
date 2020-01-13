estimateCF<-function(graph, estimateNames=FALSE){
if (!any(class(graph)=="Pathway")) stop("graph must be of class 'Pathway',found:", paste(class(graph)))
  e<-edges(graph,"proteins")
  n<-length(nodes(graph,"proteins"))
  nod<-nodes(graph, "proteins")
  nodid<-unname(sapply(nod, function(x) strsplit(x,":")[[1]][2]))
  
  
  #complexex
e.sub<-e[e[,"direction"]=="undirected" & e[,"type"]=="Binding",]
compNod<-unique(c(e.sub[,"src"],e.sub[,"dest"]))

if (length(compNod)>0) {
  comp<-gRbase::getCliques(subGraph(compNod, buildGraphNEL(e, FALSE, NULL)))# buildGraphNEL(nod, e, TRUE))) #
  comp[sapply(comp, length)==2]<-NULL
} else comp<-NULL
if (length(comp)>0) names(comp)<-paste("complex", seq_len(length(comp)),sep="")


#families

outedg <- tapply(e[,"dest"], e[,"src"], function(x) x)[nodid]; names(outedg)<-nodid
inedg <- tapply(e[, "src"], e[, "dest"], function(x) x)[nodid]; names(inedg)<-nodid

noout<-sapply(outedg, is.null)
noin<-sapply(inedg, is.null)

outMatch<-outer(outedg, outedg, function(x,y) vapply(seq_along(x), function(i) identical(x[[i]],y[[i]]), logical(1)))
outGroups<-apply(outMatch, 1, function(x) paste(which(x), collapse=","))

inMatch<-outer(inedg, inedg, function(x,y) vapply(seq_along(x), function(i) identical(x[[i]],y[[i]]), logical(1)))
inGroups<-apply(inMatch, 1, function(x) paste(which(x), collapse=","))

io<-cbind(inGroups, outGroups)

splint<-function(x){
paste(intersect(strsplit(x,",")[[1]], strsplit(x,",")[[2]]), collapse=",")
}

famMid<-unique(apply(io, 1, splint))
famMid<-lapply(famMid, function(x) nodid[as.numeric(strsplit(x,",")[[1]])])

fam<-famMid

single.fam<-sapply(fam, function(x) length(x)==1)
fam<-fam[!single.fam]

if (estimateNames) {
namesF<-sapply(fam, function(X) paste(suppressWarnings(Reduce(f<-function(x,y){x[x %in% y]},strsplit(X,""))), collapse=""))
w<-which(nchar(namesF)==0)
if (length(w)>0) warning("Some representative names could not be estimated")
if (length(namesF) != length(unique(namesF))) warning("Found duplicities in the representative names")
if (any(namesF %in% nodid)) warning("Found duplicities in the representative names and node names")

w<-which(nchar(namesF)==0 | namesF %in% nodid | duplicated(namesF) )
    
if (length(w)>0) namesF[w]<-sapply(fam[w], function(x) paste(x, collapse="\n"))
names(fam)<-namesF
} else names(fam)<-sapply(fam, function(x) paste(x, collapse="\n")) #names(fam)<-paste("family", seq_len(length(fam)), sep="")

out<-list(complexes=comp, families=fam)
return(out)
}
