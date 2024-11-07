# additive helpers ---------------------------
form_Bmat <- function(mat, basis_f, normalized = TRUE,
                      center = NULL, scale = NULL) {
  mat <- as.matrix(mat)
  if (!normalized) mat <- normalize(mat, center, scale)
  p <- dim(mat)[2]
  Bmat <- do.call(cbind,
    lapply(1:ncol(mat),
           function(j) {fda::eval.basis(mat[,j],basis_f[[j]])[,-rm_id]} )
  )
  # rm_id should be defined in the outermost function environment
  return(Bmat)
}

normalize <- function(theta_mat, center_th, scale_th) {
  theta_mat <- as.matrix(theta_mat)
  theta_mat <- sweep(theta_mat, 2, center_th, FUN = "-")
  theta_mat <- sweep(theta_mat, 2, scale_th, FUN = "/")
  return(g.fun(theta_mat))
}

## g fun: maps [-1,1] (probably wider) -> [0,1]

if (TRANSFORMER=="sigmoid") {
  g.fun <- function(x) { 1/(1+exp(-x)) }
  dg.fun <- function(x) {
    g.value <- g.fun(x)
    return(g.value * (1 - g.value))
  }
  hj.fun <- function(x, center_j, scale_j) { (x - center_j) / scale_j }
  Tj.fun <- function(x, center_j, scale_j) { g.fun(hj.fun(x, center_j, scale_j)) }
  dTj.fun <- function(T.value, scale_j) { T.value * (1 - T.value) / scale_j }
  d2Tj.fun <- function(T.value, scale_j) { (1 - 2*T.value) * T.value * (1 - T.value) / (scale_j^2) }
}

if (TRANSFORMER=="linear") {
  g.fun <- function(x) {(x+1)/2}
  dg.fun <- function(x) {x[] <- 1/2; x}
  hj.fun <- function(x, center_j, scale_j) { (x - center_j) / scale_j }
  Tj.fun <- function(x, center_j, scale_j) { g.fun(hj.fun(x, center_j, scale_j)) }
  dTj.fun <- function(T.value, scale_j) {dg.fun(T.value)/scale_j}
  d2Tj.fun <- function(T.value, scale_j) {T.value[] <- 0; T.value}
}


# h.fun <- function(X, Xcenter, Xscale) { sweep(sweep(X,2,Xcenter,"-"),2,Xscale,"/") }
# T.fun <- function(X, Xcenter, Xscale) { g.fun(h.fun(X,Xcenter,Xscale)) }
# dT.fun <- function(Tmat, Xscale) { sweep(Tmat*(1-Tmat),2,Xscale,"/") }
# d2T.fun <- function(Tmat, Xscale) { sweep((1-2*Tmat)*Tmat*(1-Tmat),2,Xscale,"/") }

