source("./R/utils/helper.R", local=TRUE)


# determine ydata's distribution by b()

### a transformation function
# currently use: sigmoid function



partial_Gj_cj <- function(j, gam0j, gam1j, c_mat, inputs){ # CHECKED
  basis_f <- inputs[['basis_f']]
  center_th <- inputs[['center_th']]
  scale_th <- inputs[['scale_th']]
  quadwts <- inputs[['quadwts']]
  lam_th <- inputs[['lam_th']]
  db <- inputs[['db']]
  d2b <- inputs[['d2b']]
  
  yj <- inputs$ydata[,j]
  nt <- inputs$nt
  
  # likelihood term ~ tvals
  theta_j <- as.vector(inputs$Hmat %*% c_mat[,j])
  # term_l <- (-1/Nt) * as.vector((yj - db(theta_j)) %*% inputs$Hmat)  # nbasis_th
  
  
  # ode-fidelity term ~ t4int
  theta_mat <- inputs$H4int %*% c_mat  # nt4int x p
  dtheta_j <- as.vector(inputs$dH4int %*% c_mat[,j])
  
  # normalize
  T_theta_mat <- normalize(theta_mat, center_th, scale_th)
  T_theta_j <- T_theta_mat[,j]
  B_mat <- form_Bmat(T_theta_mat, basis_f)
  dB_j <- eval.basis(T_theta_j, basis_f, 1)
  
  term_f1 <- dtheta_j - gam0j - as.vector(B_mat %*% as.vector(gam1j))  # nt4int
  term_f2 <- inputs$dH4int -
    rowProd(dB_j%*%gam1j[,j]*dTj.fun(T_theta_j,scale_th[j]), inputs$H4int)  # nt4int x nbasis_th
  term_f <- 2 * lam_th * as.vector((quadwts * term_f1) %*% term_f2)  # nbasis_th
  
  return(term_l + term_f)
}



partial2_Gj_cj <- function(j, gam0j, gam1j, c_mat, inputs) { # CHECKED
  basis_f <- inputs[['basis_f']]
  center_th <- inputs[['center_th']]
  scale_th <- inputs[['scale_th']]
  quadwts <- inputs[['quadwts']]
  lam_th <- inputs[['lam_th']]
  db <- inputs[['db']]
  d2b <- inputs[['d2b']]
  
  yj <- inputs$ydata[,j]
  nt <- inputs$nt
  
  # likelihood term ~ tvals
  theta_j <- as.vector(inputs$Hmat %*% c_mat[,j])
  term_l <- 1/nt * t(inputs$Hmat) %*% rowProd(d2b(theta_j), inputs$Hmat)
  
  # ode fidelity term ~ t4int
  theta_mat <- inputs$H4int %*% c_mat  # nt4int x p
  dtheta_j <- as.vector(inputs$dH4int %*% c_mat[,j])
  # normalize
  T_theta_mat <- normalize(theta_mat, center_th, scale_th)
  T_theta_j <- T_theta_mat[,j]
  dT_theta_j <- dTj.fun(T_theta_j, scale_th[j])
  d2T_theta_j <- d2Tj.fun(T_theta_j, scale_th[j])
  # B_mat
  B_mat <- form_Bmat(T_theta_mat, basis_f) 
  dB_j <- fda::eval.basis(T_theta_j, basis_f, 1)
  d2B_j <- fda::eval.basis(T_theta_j, basis_f, 2)
  
  A <- inputs$dH4int - 
    rowProd(dB_j %*% gam1j[,j] * dT_theta_j, inputs$H4int)  # nt4int x nbasis_f
  term_f1 <- t(A) %*% rowProd(quadwts, A)
  
  w1 <- dtheta_j - gam0j - as.vector(B_mat %*% as.vector(gam1j))  # nt4int
  tmp_mat <- rowProd((dT_theta_j)^2, d2B_j) + rowProd(d2T_theta_j, dB_j)  # nt4int x nbasis_f
  w2 <- (-1) * as.vector(tmp_mat %*% gam1j[,j])  # nt4int
  term_f2 <- t(inputs$H4int) %*% rowProd(w1*w2*quadwts, inputs$H4int)  # nbasis_f x nbasis_f
  
  return(term_l + 2 * lam_th *(term_f1 + term_f2))
}



# partial2_Gj_cj_gam0j <- function(j, gam0j, gam1j, c_mat, inputs) {
#   center_th <- inputs[['center_th']]
#   scale_th <- inputs[['scale_th']]
#   basis_f <- inputs[['basis_f']]
#   quadwts <- inputs[['quadwts']]
#   lam_th <- inputs[['lam_th']]
#   enlarge_amount <- inputs[['enlarge_amount']]
#   db <- inputs[['db']]
#   d2b <- inputs[['d2b']]
#   
#   p <- dim(ydata)[2]
#   nbasis_th <- dim(H4int)[2]
#   
#   # ode fidelity term ~ t4int
#   theta_mat <- H4int %*% c_mat  # nt4int x p
#   dtheta_j <- as.vector(dH4int %*% c_mat[,j])
#   
#   # normalize
#   T_theta_mat <- normalize(theta_mat, center_th, scale_th)
#   T_theta_j <- T_theta_mat[,j]
#   
#   B_mat <- form_Bmat(T_theta_mat, basis_f)
#   dB_j <- fda::eval.basis(T_theta_j, basis_f, 1)
#   
#   A <- dH4int - rowProd(dB_j %*% gam1j[,j] * dTj.fun(T_theta_j, scale_th[j]), H4int)  # nt4int x nbasis_th
#   
#   # partial2_Gj_cj_gam0j
#   term_f <- -2 * lam_th * as.vector(quadwts %*% A)
#   
#   return(term_f)
# }
# 
# partial2_Gj_cj_gam1jk <- function(j, k, gam0j, gam1j, c_mat, inputs) {
#   # Input:
#   #   j: current dim
#   #   gam0j: scalar;  gam1j: nbasis x p
#   #   c_mat: nbasis_th x p
#   #   inputs: general inputs for constant settings
#   ydata <- inputs[['ydata']]
#   Hmat <- inputs[['Hmat']]
#   H4int <- inputs[['H4int']]
#   dH4int <- inputs[['dH4int']]
#   center_th <- inputs[['center_th']]
#   scale_th <- inputs[['scale_th']]
#   basis_f <- inputs[['basis_f']]
#   nbasis_f <- inputs[['nbasis_f']]
#   quadwts <- inputs[['quadwts']]
#   lam_th <- inputs[['lam_th']]
#   enlarge_amount <- inputs[['enlarge_amount']]
#   db <- inputs[['db']]
#   d2b <- inputs[['d2b']]
#   
#   p <- dim(ydata)[2]
#   nbasis_th <- dim(H4int)[2]
#   
#   # ode fidelity term ~ t4int
#   theta_mat <- H4int %*% c_mat  # nt4int x p
#   dtheta_j <- as.vector(dH4int %*% c_mat[,j])
#   
#   # normalize
#   T_theta_mat <- normalize(theta_mat, center_th, scale_th)
#   T_theta_j <- T_theta_mat[,j]
#   
#   B_mat <- form_Bmat(T_theta_mat, basis_f)
#   dB_j <- fda::eval.basis(T_theta_j, basis_f, 1)
#   
#   A <- dH4int - rowProd(dB_j %*% gam1j[,j] * dTj.fun(T_theta_j, scale_th[j]), H4int)  # nt4int x nbasis_th
#   
#   # partial2_Gj_cj_gam1jk
#   B_k <- B_mat[,((k-1)*nbasis_f+1):(k*nbasis_f)]
#   term_f <- -2 * lam_th * t(A) %*% rowProd(quadwts, B_k)
#   
#   return(term_f)
# }
# 
# 
# partial2_Gj_cj_gam1jj <- function(j, gam0j, gam1j, c_mat, inputs) {
#   # Input:
#   #   j: current dim
#   #   gam0j: scalar;  gam1j: nbasis x p
#   #   c_mat: nbasis_th x p
#   #   inputs: general inputs for constant settings
#   ydata <- inputs[['ydata']]
#   Hmat <- inputs[['Hmat']]
#   H4int <- inputs[['H4int']]
#   dH4int <- inputs[['dH4int']]
#   center_th <- inputs[['center_th']]
#   scale_th <- inputs[['scale_th']]
#   basis_f <- inputs[['basis_f']]
#   nbasis_f <- inputs[['nbasis_f']]
#   quadwts <- inputs[['quadwts']]
#   lam_th <- inputs[['lam_th']]
#   enlarge_amount <- inputs[['enlarge_amount']]
#   db <- inputs[['db']]
#   d2b <- inputs[['d2b']]
#   
#   p <- dim(ydata)[2]
#   nbasis_th <- dim(H4int)[2]
#   
#   # ode fidelity term ~ t4int
#   theta_mat <- H4int %*% c_mat  # nt4int x p
#   dtheta_j <- as.vector(dH4int %*% c_mat[,j])
#   
#   # normalize
#   T_theta_mat <- normalize(theta_mat, center_th, scale_th)
#   T_theta_j <- T_theta_mat[,j]
#   dT_theta_j <- dTj.fun(T_theta_j, scale_th[j])
#   
#   B_mat <- form_Bmat(T_theta_mat, basis_f)
#   dB_j <- fda::eval.basis(T_theta_j, basis_f, 1)
#   
#   A <- dH4int - rowProd(dB_j %*% gam1j[,j] * dTj.fun(T_theta_j, scale_th[j]), H4int)  # nt4int x nbasis_th
#   
#   # partial2_Gj_cj_gam1jj
#   B_j <- B_mat[,((j-1)*nbasis_f+1):(j*nbasis_f)]
#   w <- (dtheta_j - gam0j - as.vector(B_mat %*% as.vector(gam1j))) * dT_theta_j
#   term_f <- -2 * lam_th * (
#     t(A) %*% rowProd(quadwts, B_j) + t(H4int) %*% rowProd(w * quadwts, dB_j)
#   )
#   
#   return(term_f)
# }


partial2_Gj_cj_gamj <- function(j, gam0j, gam1j, c_mat, inputs) {
  center_th <- inputs[['center_th']]
  scale_th <- inputs[['scale_th']]
  basis_f <- inputs[['basis_f']]
  nbasis_f <- inputs[['nbasis_f']]
  quadwts <- inputs[['quadwts']]
  lam_th <- inputs[['lam_th']]
  db <- inputs[['db']]
  d2b <- inputs[['d2b']]
  
  p <- inputs$p
  nbasis_th <- inputs$nbasis_th
  
  # ode fidelity term ~ t4int
  theta_mat <- inputs$H4int %*% c_mat  # nt4int x p
  dtheta_j <- as.vector(inputs$dH4int %*% c_mat[,j])
  
  # normalize
  T_theta_mat <- normalize(theta_mat, center_th, scale_th)
  T_theta_j <- T_theta_mat[,j]
  dT_theta_j <- dTj.fun(T_theta_j, scale_th[j])
  
  B_mat <- form_Bmat(T_theta_mat, basis_f)
  dB_j <- fda::eval.basis(T_theta_j, basis_f, 1)
  
  term_f1 <- dtheta_j - gam0j - as.vector(B_mat %*% as.vector(gam1j))  # nt4int
  A <- dH4int - rowProd(dB_j %*% gam1j[,j] * dTj.fun(T_theta_j, scale_th[j]),
                        inputs$H4int)  # nt4int x nbasis_th
  
  # partial2_Gj_cj_gam0j
  term_f <- as.vector(quadwts %*% A)
  # partial2_Gj_cj_gam1jk
  for (k in 1:p) {
    B_k <- B_mat[,((k-1)*nbasis_f+1):(k*nbasis_f)]
    term_f <- cbind(term_f, t(A) %*% rowProd(quadwts, B_k))
  }
  if (dim(term_f)[1]!=nbasis_th || dim(term_f)[2]!=(nbasis_f*p+1)) stop("wrong dim for partial2_Gj_cj_gamj()!")
  # partial2_Gj_cj_gam1jj has an additional term, add it on.
  w <- (dtheta_j - gam0j - as.vector(B_mat %*% as.vector(gam1j))) * dT_theta_j
  term_f[,((j-1)*nbasis_f+2):(j*nbasis_f+1)] <- term_f[,((j-1)*nbasis_f+1):(j*nbasis_f)] +
    t(inputs$H4int) %*% rowProd(w * quadwts, dB_j)
  
  term_f <- -2 * lam_th * term_f
  return(term_f)
}


# partial2_Gj_cj_gamj.num <- function(j, gam0j, gam1j, c_mat, inputs) {
#   return(
#     pracma::jacobian(f = function(gamj) {
#       gam0j <- gamj[1]
#       gam1j <- matrix(gamj[-1], nbasis_f, p)
#       partial_Gj_cj(j, gam0j, gam1j, c_mat, inputs)
#     }, x0 = c(gam0j, gam1j))
#   )
# }




#### Objective functions

Gj <- function(j, gam0j, gam1j, c_mat, inputs) {  # CHECKED
  basis_f <- inputs[['basis_f']]
  center_th <- inputs[['center_th']]
  scale_th <- inputs[['scale_th']]
  quadwts <- inputs[['quadwts']]
  lam_th <- inputs[['lam_th']]
  b <- inputs[['b']]
  
  yj <- inputs$ydata[,j]
  nt <- inputs$nt
  
  # likelihood term ~ tvals
  theta_j <- as.vector(inputs$Hmat %*% c_mat[,j])
  term_l <- (-1/nt) * sum(yj * theta_j - b(theta_j))
  
  # ode fidelity term ~ t4int
  dtheta_j <- as.vector(inputs$dH4int %*% c_mat[,j])
  theta_mat <- inputs$H4int %*% c_mat
  T_theta_mat <- normalize(theta_mat, center_th, scale_th)
  
  B_mat <- form_Bmat(T_theta_mat, basis_f) 
  
  term_f <- lam_th * sum(
    quadwts * (dtheta_j - gam0j - as.vector(B_mat %*% as.vector(gam1j)))^2
  )
  
  return(term_l + term_f)
}


# Hj <- function(j, gam0j, gam1j, c_mat, inputs, reg = "l1") {
#   ydata <- inputs[['ydata']]
#   Hmat <- inputs[['Hmat']]
#   basis_f <- inputs[['basis_f']]
#   lam_gam <- inputs[['lam_gam']]
#   b <- inputs[['b']]
#   
#   
#   theta_j <- as.vector(Hmat %*% c_mat[,j])
#   term_l <- (-1/nt) * sum(yj * theta_j - b(theta_j))
#   
#   if (reg == "l1") {
#     return(term_l + lam_gam * sum(c(gam0j, as.vector(gam1j))^2))
#   } else if (reg == "l2") {
#     return(term_l + lam_gam * sum(c(gam0j, as.vector(gam1j))^2))
#   } else if (reg == "none") {
#     return(term_l)
#   }
# }



