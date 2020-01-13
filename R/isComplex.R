isComplex<-function (graph, x) 
{   
if (!any(class(graph)=="Pathway")) stop("graph must be of class 'Pathway',found:", paste(class(graph)))    
    e <- edges(graph,"proteins")
    sube <- e[e[, "src"] %in% x & e[, "dest"] %in% x, ]
    nrow(sube) == ncol(combn(x,2)) & all(x %in% unlist(sube[, c("src","dest")])) & 
        all(sube[, "direction"] == "undirected") & 
        all( (regexpr("Binding",sube[, "type"]) != -1) )
}

#example
#g<-convertIdentifiers(kegg[["p53 signaling pathway"]], "symbol")
#x<-c("CCNE1", "CDK2",  "CCNE2")
#isComplex(g,x)
#
#
#x<-c(x, nodes(g)[15])
#isComplex(g,x)

###########




