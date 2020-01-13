.buildGraphNEL<-function (edges, sym, edge.types=NULL)
{
  if (!is.null(edge.types))
    edges <- .selectEdges(edges, edge.types)
  if (nrow(edges) == 0)
    g <- new("graphNEL", character(), list(), "directed")
  else {
    prep <- .prepareEdges(edges, sym)
    nodes <- union(unique(prep$src), unique(prep$dest))
    g <- new("graphNEL", nodes, .edLi(nodes, prep), "directed")
    graph::edgeDataDefaults(g, "edgeType") <- "undefined"
    graph::edgeData(g, prep$src, prep$dest, "edgeType") <- prep$type
  }
  return(g)
}

.selectEdges<-function (m, types)
{
  missing <- setdiff(types, edgeAttrs[[1]][,1])
  if (length(missing) > 0) {
    warning("the following edge types will be ommited: ", paste(sort(missing),
                                                         collapse = ", "))
  }
  m[m$type %in% types, , drop = FALSE]
}

.prepareEdges<-function (e, sym)
{
  e[] <- lapply(e, as.character)
  if (sym) {
    e <- .symmetric(e)
  }
  ends <- .endpoints(e)
  types <- tapply(e$type, ends, function(group) {
    paste(sort(unique(group)), collapse = ";")
  })
  binder <- function(...) rbind.data.frame(..., stringsAsFactors = FALSE)
  merged <- do.call(binder, strsplit(names(types), "|", fixed = TRUE))
  colnames(merged) <- c("src", "dest")
  cbind(merged, data.frame(type = as.character(types), stringsAsFactors = FALSE))
}

.symmetric<-function (e)
{
  mask <- e$direction == "undirected" & (e$src_type != e$dest_type | e$src != e$dest)
  dird <- e[!mask, ]
  undir <- e[mask, ]
  revdir <- e[mask, ]
  revdir$src_type <- undir$dest_type
  revdir$src <- undir$dest
  revdir$dest_type <- undir$src_type
  revdir$dest <- undir$src
  rbind(dird, undir, revdir)
}


.endpoints<-function (e)
{
  paste0( e$src_type,":",e$src,"|", e$dest_type,":", e$dest)
}




.edLi <- function(n, e) {
  sapply(n, function(n) list(edges = e[e[, 1] == n, 2]),
         simplify=FALSE, USE.NAMES = TRUE)
}

#Rgraphviz
getRenderPar<-function (g, name, what = c("nodes", "edges", "graph")) 
{
  what <- match.arg(what)
  nms <- switch(what, nodes = nodes(g), edges = edgeNames(g, 
                                                          recipEdges = graphRenderInfo(g, "recipEdges")), graph = "graph")
  ans <- switch(what, nodes = nodeRenderInfo(g, name), edges = edgeRenderInfo(g, 
                                                                              name), graph = graphRenderInfo(g, name))
  if (!is.null(ans) && !any(is.na(ans))) {
    if (!is.null(names(ans))) 
      ans <- ans[nms]
  }    else {
    default <- parRenderInfo(g, what)[[name]][1]
    if (is.null(default)) 
      default <- graph.par.get(what)[[name]][1]
    if (is.null(ans)) {
      ans <- rep(default, length(nms))
    }        else {
      if (!is.null(default)) 
        ans[is.na(ans)] <- default
      ans <- ans[nms]
    }
  }
  ans
}

renderSpline<-function (spline, arrowhead = FALSE, arrowtail = FALSE, len = 1, 
                        col = "black", lwd = 1, lty = "solid", bbox, ...) 
{
  mylty <- as.numeric(lty)
  if (!is.na(mylty)) 
    lty <- mylty
  lapply(spline, lines, col = col, lwd = lwd, lty = lty, ...)
  warn <- FALSE
  xyhead <- tail(bezierPoints(spline[[length(spline)]]), 2)
  if (is.function(arrowhead[[1]])) {
    xy <- list(x = xyhead[2, 1], y = xyhead[2, 2])
    try(arrowhead[[1]](xy, col = col, lwd = lwd, lty = lty))
  }    else {
    warn <- drawHead(arrowhead, xyhead, bbox, col, lwd, lty, 
                     len, out = TRUE)
  }
  xytail <- head(bezierPoints(spline[[length(spline)]]), 2)
  if (is.function(arrowtail[[1]])) {
    xy <- list(x = xytail[1, 1], y = xytail[1, 2])
    try(arrowtail[[1]](xy, col = col, lwd = lwd, lty = lty))
  }    else {
    warn <- warn | drawHead(arrowtail, xytail[2:1, ], bbox, 
                            col, lwd, lty, len, out = FALSE)
  }
  warn
}
