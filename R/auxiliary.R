#' Functions to extract and display results
#'
#' Functions to extract and display results
#' @aliases res
#' @aliases topo.sig
#' @aliases degtable
#'
#' 
#' @param object,x an output from one of following function \code{"SPIA","PRS", "CePA", "PWEA", "TAPPA", "TopologyGSA","DEGraph", "clipper" }
#' @param ... other arguments
#' 
#' @export
res<-function(object){
  UseMethod("res")
}
#' @export
topo.sig<-function(object){
  UseMethod("topo.sig")
}
#' @export
degtable<-function(object){
  UseMethod("degtable")
}


#' @export
res.topResult<-function(object){
  return(object$res)
}
#' @export
topo.sig.topResult<-function(object){
  return(object$topo.sig)
}
#' @export
degtable.topResult<-function(object){
  return(object$degtable)
}
#' @export
print.topResult<-function(x,...){
  cat("$res\n")
  print(head(x$res))
  cat("... print truncated\n")
  cat("$topo.sig\n")
  print(x$topo.sig[1:3])
  cat("... print truncated\n")
  cat("$degtable\n")
  print(head(x$degtable))
  cat("... print truncated\n")
}


.catchErr<-function (l, f){
  log <- lapply(l, function(x) {
    tryCatch(list("ok", f(x)),
             error = function(e) list("err", e))})
  list(results = Filter(Negate(is.null),
                        .filterByTag("ok", log)),
       errors = sapply(.filterByTag("err", log), gettext))
}


.filterByTag <-function (tag, l)
{
  isTagged <- sapply(l, function(x) x[[1]] == tag)
  lapply(l[isTagged], function(x) x[[2]])
}


