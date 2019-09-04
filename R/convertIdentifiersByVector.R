#' Convert pathway identifiers
#'
#' This function converts identifiers of nodes in a pathway. It uses the user specified named vector for the conversion.
#'
#' @param pathway An object of class \code{Pathway}
#' @param  conv.table A data.frame, in which the first column contains the type
#' and the identifiers present in the pathway separeated by \code{:}
#' and the second column contains the new identifiers and
#' the third columns contains the types of the new identifiers
#' @return A \code{Pathway} with new identifiers of the nodes
#' @author Ivana Ihnatova
#' @seealso \code{\link[graphite]{Pathway-class}}
#' @keywords manip
#' @examples
#' g<-pathways("hsapiens","kegg")
#' ng<-sapply(g, function(x) length(nodes(x,"mixed")))
#' g<-g[[which.min(ng)]]
#' conv<-data.frame(orig=nodes(g,"mixed"), new=LETTERS[seq_len(min(ng))],newtype=rep("LETTERS",min(ng)))
#' gc<-convertIdentifiersByVector(g, conv.table = conv)@protEdges
#'
#' @export



convertIdentifiersByVector<-function (pathway, conv.table)
{
    ns <- graphite::nodes(pathway,"mixed")
    if (!all(ns %in% as.character(conv.table[,1])))
      warning(paste("These pathway nodes are missing in the 'conv.table':",
                    paste(ns[! ns %in% as.character(conv.table[,1])], collapse=", "),
                    ". The original identifiers were kept."))


    miss <- ns[! ns %in% as.character(conv.table[,1])]
    missT <- sapply(miss, function(x) strsplit(x,":")[[1]][2:1])

    conv.table<-as.matrix(conv.table)
    if (length(miss)>0) conv.table<-rbind(conv.table, cbind(miss,t(missT)))

    conv.table<-cbind(conv.table, t(sapply(as.character(conv.table[,1]), function(x) strsplit(x,":")[[1]])))
    rownames(conv.table)<-conv.table[,5]


    for (which in c("protEdges","protPropEdges","metabolEdges","metabolPropEdges","mixedEdges")) {
      es <- slot(pathway,which)
      if (nrow(es)>0) {
        es[,"src_type"]<- factor(conv.table[es[,"src"],3])
        es[,"dest_type"]<- factor(conv.table[es[,"dest"],3])
        es[,"src"]<-conv.table[es[,"src"],2]
        es[,"dest"]<- conv.table[es[,"dest"],2]

        slot(pathway,which)<-es
      }

    }

    return(pathway)
}

