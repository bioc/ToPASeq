#from qpgraph
qpIPF<-function (vv, clqlst, tol = 0.001, verbose = FALSE) 
{
    vv <- as(vv, "matrix")
    if (verbose) {
        n.var <- nrow(vv)
        if (clqlst[[1]][1] <= n.var) {
            n.var <- 0
        }
        cat(paste(paste("qpIPF: ", length(clqlst) - n.var), "cliques\n"))
    }
    V <- diag(length(vv[, 1]))
    precision <- 1
    while (precision > tol) {
        Vold <- V
        V <- .qpIPFpass(vv, V, clqlst)
        precision <- max(abs(V - Vold))
        if (verbose) 
            cat("qpIPF: precision =", precision, "\n")
    }
    return(as(V, "dspMatrix"))
}

.qpIPFstep<-function (Vf, Vn, wh, clqlst) 
{
    a <- clqlst[[wh]]
    b <- (1:length(Vf[, 1]))[-a]
    Vfaa <- Vf[a, a]
    Vni <- solve(Vn[a, a])
    Bnba <- Vn[b, a] %*% Vni
    Vnbba <- Vn[b, b] - Vn[b, a] %*% Vni %*% Vn[a, b]
    V <- Vf
    V[b, a] <- Bnba %*% Vfaa
    V[a, b] <- t(V[b, a])
    V[b, b] <- Vnbba + Bnba %*% Vfaa %*% t(Bnba)
    return(V)
}

.qpIPFpass<-function (Vf, Vn, clqlst) 
{
  n.var <- nrow(Vf)
  firstclq <- 1
  if (clqlst[[1]][1] > n.var) {
    firstclq <- n.var + 1
  }
  V <- Vn
  for (i in firstclq:length(clqlst)) {
    V <- .qpIPFstep(Vf, V, i, clqlst)
  }
  return(V)
}