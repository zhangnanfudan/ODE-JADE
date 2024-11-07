source("./R/utils/helper.R", local = TRUE)
source("./R/utils/utils.R", local = TRUE)
source("./R/utils/transformer.R", local = TRUE)


### optimize functions =====================


descent.theta <- function(j, gam0, gam1, Cmat, inputs, lam_th, con) {
  cj_old <- Cmat[, j]
  pGpc <- partial_Gj_cj(j, gam0, gam1, Cmat, inputs, lam_th)
  cj <- cj_old - con$stepsize * pGpc
  return(cj)
}


estimate.theta <- function(j, gam0, gam1, Cmat, inputs, lam_th, con) {
  if (con$theta_method == "1gd") {
    cj <- Cmat[, j]
    pGpc <- partial_Gj_cj(j, gam0, gam1, Cmat, inputs, lam_th)
    d <- con$stepsize * pGpc
    while (TRUE) { ## Armijo rule
      pred <- sum(pGpc * d)
      ared <- Gj(j, gam0, gam1, Cmat, inputs, lam_th) -
        Gj(j, gam0, gam1, newCol(Cmat, j, cj - d), inputs, lam_th)
      if (pred * con$amj.alpha <= ared) {
        Cmat[, j] <- cj - d
        break
      } else {
        d <- d * (con$amj.beta)
      }
      if (max(abs(d)) < con$eps.conv) break
    }
    return(Cmat[, j])
  } else if (con$theta_method == "1dH") {
    cj <- Cmat[, j]
    pGpc <- partial_Gj_cj(j, gam0, gam1, Cmat, inputs, lam_th)
    p2Gpc2 <- partial2_Gj_cj(j, gam0, gam1, Cmat, inputs, lam_th)
    diagH <- diag(p2Gpc2)
    # if (min(diagH) < 0) browser()
    diagH <- pmin(pmax(diagH, 1e-2), 1e9)
    d <- diagH^(-1) * pGpc
    pred <- sum(pGpc * d) - 0.5 * sum(d^2 * diagH)
    while (TRUE) { ## Armijo rule
      ared <- Gj(j, gam0, gam1, Cmat, inputs, lam_th) -
        Gj(j, gam0, gam1, newCol(Cmat, j, cj - d), inputs, lam_th)
      if (pred * con$amj.alpha <= ared) {
        Cmat[, j] <- cj - d
        break
      } else {
        d <- d * (con$amj.beta)
        pred <- pred * (con$amj.beta)
      }
      if (max(abs(d)) < con$eps.conv) break
    }
  }
  return(Cmat[, j])
}


estimate.gamma <- function(j, gam0, gam1, Cmat, inputs, lam_th, lam_gam, con) {
  p <- inputs[["p"]]
  nL <- inputs[["nL"]]
  theta <- inputs$H4int %*% Cmat
  Bmat <- form_Bmat(theta, inputs$basis_f,
    center = inputs[["center_th"]], scale = inputs[["scale_th"]],
    normalized = FALSE
  )
  Z <- Bmat
  dtheta_j <- as.vector(inputs$dH4int %*% Cmat[, j])
  fit_j <- grpSparseReg(
    X = Z, y = dtheta_j,
    group = inputs$group,
    sparse_penalty = con$sparse_penalty,
    lambda = lam_gam, intercept = TRUE,
    pen.factor = inputs$w.mat[, j]
  )
  gam1j <- as.numeric(fit_j$beta)
  gam0j <- as.numeric(fit_j$b0)
  
  # TODO: modify this
  # meanfjk <- colSums(matrix(Bc*gam1j,nrow=nL))
  # gam1jshift <- rep(meanfjk,each=nL)
  # gam1j <- gam1j - gam1jshift

  # sanity check
  # TODO: whether necessary? how to do it?
  # old_obj <- sum((dtheta_j - gam0[j] - Bmat %*% c(gam1[,,j]))^2) +
  #   lam_gam * sum(apply(gam1[,,j],2,function(x) {sqrt(sum(x^2))}))
  # new_obj <- sum((dtheta_j - gam0j - Bmat %*% c(gam1j))^2) +
  #   lam_gam * sum(apply(matrix(gam1j,nL,p),2,function(x) {sqrt(sum(x^2))}))
  # if (old_obj < new_obj) {
  #   gam1j <- c(gam1[,,j])
  #   gam0j <- gam0[j]
  # }
  return(c(gam0j, gam1j))
}


### derivatives =========================

partial_Gj_cj <- function(j, gam0, gam1, Cmat, inputs, lam_th) {
  basis_f <- inputs[["basis_f"]]
  center_th <- inputs[["center_th"]]
  scale_th <- inputs[["scale_th"]]
  quadwts <- inputs[["quadwts"]]
  db <- inputs[["db"]]
  d2b <- inputs[["d2b"]]

  yj <- inputs$ydata[, j]
  Nt <- inputs$Nt
  R <- inputs$R
  p <- inputs$p
  nt4int <- inputs$nt4int

  # likelihood term ~ tvals
  theta_j <- as.vector(inputs$Hmat %*% Cmat[, j])
  term_l <- (-1 / Nt) * as.vector((foldRep(yj, R) - R * db(theta_j)) %*% inputs$Hmat)

  # ode-fidelity term ~ t4int
  theta_mat <- inputs$H4int %*% Cmat # nt4int x p
  dtheta_mat <- inputs$dH4int %*% Cmat # nt4int x p
  T_theta_mat <- normalize(theta_mat, center_th, scale_th) # nt4int x p
  Bmat <- form_Bmat(T_theta_mat, basis_f) # nt4int x (nbasis_f*p)
  dB_j <- fda::eval.basis(T_theta_mat[, j], basis_f[[j]], 1)[, -rm_id] # nt4int x nbasis_f

  w1 <- dtheta_mat - repRbind(gam0, nt4int) - Bmat %*% keepDim(gam1, 3) # nt4int x p, 2b summed by s
  w2 <- -dB_j %*% gam1[, j, ] # nt4int x p, 2b summed by s
  w3 <- dTj.fun(T_theta_mat[, j], scale_th[j]) # nt4int
  term_f <- 2 * lam_th * (rowSums(w1 * w2) * w3 * quadwts) %*% inputs$H4int
  term_f <- as.vector(term_f + 2 * lam_th * (w1[, j] * quadwts) %*% inputs$dH4int)

  return(term_l + term_f)
}



partial2_Gj_cj <- function(j, gam0, gam1, Cmat, inputs, lam_th) { # CHECKED
  basis_f <- inputs[["basis_f"]]
  center_th <- inputs[["center_th"]]
  scale_th <- inputs[["scale_th"]]
  quadwts <- inputs[["quadwts"]]
  db <- inputs[["db"]]
  d2b <- inputs[["d2b"]]
  yj <- inputs$ydata[, j]
  Nt <- inputs$Nt
  nt <- inputs$nt
  R <- inputs$R
  p <- inputs$p
  nt4int <- inputs$nt4int

  # likelihood term ~ tvals
  theta_j <- as.vector(inputs$Hmat %*% Cmat[, j])
  term_l <- 1 / nt * t(inputs$Hmat) %*% rowProd(d2b(theta_j), inputs$Hmat)
  # browser()
  # ode fidelity term ~ t4int
  theta_mat <- inputs$H4int %*% Cmat # nt4int x p
  dtheta_mat <- inputs$dH4int %*% Cmat # nt4int x p
  T_theta_mat <- normalize(theta_mat, center_th, scale_th) # nt4int x p
  Bmat <- form_Bmat(T_theta_mat, basis_f) # nt4int x (nbasis_f*p)
  dB_j <- fda::eval.basis(T_theta_mat[, j], basis_f[[j]], 1)[, -rm_id]
  d2B_j <- fda::eval.basis(T_theta_mat[, j], basis_f[[j]], 2)[, -rm_id]
  dT_theta_j <- dTj.fun(T_theta_mat[, j], scale_th[j])
  d2T_theta_j <- d2Tj.fun(T_theta_mat[, j], scale_th[j])

  # term_f1
  w1 <- -dB_j %*% gam1[, j, ] # nt4int x p, 2b summed by s
  A <- inputs$dH4int + rowProd(w1[, j] * dT_theta_j, inputs$H4int) # nt4int x nbasis_f
  term_f11 <- t(A) %*% rowProd(quadwts, A)
  term_f12 <- t(inputs$H4int) %*% rowProd(
    rowSums(w1[, -j, drop = FALSE]^2) * dT_theta_j^2 * quadwts, inputs$H4int
  )
  term_f1 <- term_f11 + term_f12
  # term_f2
  w1 <- dtheta_mat - repRbind(gam0, nt4int) - Bmat %*% keepDim(gam1, 3) # nt4int x p, 2b summed by s
  w2 <- -(rowProd((dT_theta_j)^2, d2B_j) + rowProd(d2T_theta_j, dB_j)) %*%
    gam1[, j, ] # nt4int x p, 2b summed by s
  term_f2 <- t(inputs$H4int) %*% rowProd(rowSums(w1 * w2) * quadwts, inputs$H4int)

  return(term_l + 2 * lam_th * (term_f1 + term_f2))
}



#### Objective functions

Gj <- function(j, gam0, gam1, Cmat, inputs, lam_th) { # CHECKED
  basis_f <- inputs[["basis_f"]]
  center_th <- inputs[["center_th"]]
  scale_th <- inputs[["scale_th"]]
  quadwts <- inputs[["quadwts"]]
  b <- inputs[["b"]]
  yj <- inputs$ydata[, j]
  Nt <- inputs$Nt
  R <- inputs$R
  p <- inputs$p
  nt4int <- inputs$nt4int

  # likelihood term ~ tvals
  theta_j <- as.vector(inputs$Hmat %*% Cmat[, j])
  term_l <- (-1 / Nt) * sum(foldRep(yj, R) * theta_j - R * b(theta_j))

  # ode fidelity term ~ t4int
  dtheta_mat <- inputs$dH4int %*% Cmat
  theta_mat <- inputs$H4int %*% Cmat
  T_theta_mat <- normalize(theta_mat, center_th, scale_th)
  Bmat <- form_Bmat(T_theta_mat, basis_f)
  term_f <- as.numeric(lam_th * quadwts %*%
    rowSums((dtheta_mat - repRbind(gam0, nt4int) - Bmat %*% keepDim(gam1, 3))^2))

  return(term_l + term_f)
}


G <- function(gam0, gam1, Cmat, inputs, lam_th) {
  basis_f <- inputs[["basis_f"]]
  center_th <- inputs[["center_th"]]
  scale_th <- inputs[["scale_th"]]
  quadwts <- inputs[["quadwts"]]
  b <- inputs[["b"]]
  Nt <- inputs$Nt
  nt4int <- inputs$nt4int
  R <- inputs$R
  p <- inputs$p
  # likelihood term ~ tvals
  theta_mat <- inputs$Hmat %*% Cmat
  YR <- apply(inputs$ydata, 2, foldRep, R = R)
  term_l <- (-1 / Nt) * sum(YR * theta_mat - R * b(theta_mat))
  # ode fidelity term ~ t4int
  dtheta_mat <- inputs$dH4int %*% Cmat
  theta_mat <- inputs$H4int %*% Cmat
  T_theta_mat <- normalize(theta_mat, center_th, scale_th)
  Bmat <- form_Bmat(T_theta_mat, basis_f)
  term_f <- as.numeric(lam_th * quadwts %*%
    rowSums((dtheta_mat - repRbind(gam0, nt4int) - Bmat %*% keepDim(gam1, 3))^2))
  return(term_l + term_f)
}



GH <- function(gam0, gam1, Cmat, inputs, lam_th, lam_gam) {
  basis_f <- inputs[["basis_f"]]
  center_th <- inputs[["center_th"]]
  scale_th <- inputs[["scale_th"]]
  quadwts <- inputs[["quadwts"]]
  b <- inputs[["b"]]
  Nt <- inputs$Nt
  nt4int <- inputs$nt4int
  R <- inputs$R
  p <- inputs$p
  # likelihood term ~ tvals
  theta_mat <- inputs$Hmat %*% Cmat
  YR <- apply(inputs$ydata, 2, foldRep, R = R)
  term_l <- (-1 / Nt) * sum(YR * theta_mat - R * b(theta_mat))
  # ode fidelity term ~ t4int
  dtheta_mat <- inputs$dH4int %*% Cmat
  theta_mat <- inputs$H4int %*% Cmat
  T_theta_mat <- normalize(theta_mat, center_th, scale_th)
  Bmat <- form_Bmat(T_theta_mat, basis_f)
  term_f <- as.numeric(lam_th * quadwts %*%
    rowSums((dtheta_mat - repRbind(gam0, nt4int) - Bmat %*% keepDim(gam1, 3))^2))
  term_p <- lam_th * sum((inputs[["w.mat"]] * apply(gam1, c(2, 3), function(x) sqrt(sum(x^2)))) %*% lam_gam)
  return(term_l + term_f + term_p)
}