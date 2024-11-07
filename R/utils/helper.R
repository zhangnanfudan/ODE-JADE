
# computation helpers -----------------------
rowProd <- function(v, A) {
  # diag(v) %*% A
  v <- as.vector(v)
  if (!is.matrix(A)) stop("A is not matrix")
  if (nrow(A) != length(v)) stop("v and A, dims not matched")
  return(v * A)
}
# NOTE: colSums() and rowSums() in base R, their names have mislead me
# rowSums is "the sums for each row" instead of "sum all the rows together"
colProd <- function(A, v) {
  # A %*% diag(v)
  v <- as.vector(v)
  if (!is.matrix(A)) stop("A is not matrix")
  if (ncol(A) != length(v)) stop("v and A, dims not matched")
  return(sweep(A, 2, v, FUN="*"))
}

foldRep <- function(v, R, foldmethod = "sum") {
  return(apply(matrix(v, ncol = R),1,foldmethod))
}

repRbind <- function(v, R) {
  return(do.call(rbind,replicate(R,v,simplify=FALSE)))
}

repCbind <- function(v, R) {
  return(do.call(cbind,replicate(R,v,simplify=FALSE)))
}

keepDim <- function(a, dims) {
  return(apply(a,dims,as.vector))
}

newCol <- function(mat,j,x) {
  mat[,j] <- x
  return(mat)
}

selected_group <- function(v,grpsize) {
  if (length(v)%%gprsize!=0) stop("length(v) is not a multiple of grpsize.")
  nr <- length(v)/grpsize
  v <- matrix(v, nrow = grpsize)
  v <- v!=0
  return(apply(v,2,any))
}
