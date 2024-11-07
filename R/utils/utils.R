source("./R/utils/transformer.R", local = TRUE)



# computation tools ------------------
get.quadwts <- function(x) {
  quadts <- x
  nquad <- length(quadts)
  quadwts <- rep(1, nquad)
  even.ind <- seq(2, (nquad - 1), by = 2)
  odd.ind <- seq(3, (nquad - 2), by = 2)
  quadwts[even.ind] <- 4
  quadwts[odd.ind] <- 2
  h <- quadts[2] - quadts[1]
  quadwts <- quadwts * (h / 3)
  return(quadwts)
}



# metric helpers -------------------------

neg_log_likelihood <- function(ydata, theta_est, fam) {
  b <- fam$b
  # negative log-likelihood for exp.family
  theta_est <- as.matrix(theta_est)
  ydata <- as.matrix(ydata)
  nt <- dim(theta_est)[1]
  p <- dim(theta_est)[2]
  if (dim(ydata)[1] != nt || dim(ydata)[2] != p) stop("preds, y, df not matched")
  btheta_est <- b(theta_est)
  return(-1 / (nt * p) * sum(ydata * theta_est - btheta_est)) # or -1/(nt)
}


get.fidelity <- function(Cmat, gam0, gam1, inputs,
                         type = c("all", "data", "ode")) {
  type <- match.arg(type)

  ydata <- inputs$ydata
  Hmat <- inputs$Hmat
  H4int <- inputs$H4int
  dH4int <- inputs$dH4int
  p <- dim(ydata)[2]
  nt <- dim(ydata)[1]
  R <- inputs$R

  theta_bspl <- Hmat %*% Cmat
  theta_bspl_4int <- H4int %*% Cmat
  dtheta_bspl_4int <- dH4int %*% Cmat
  Bmat <- form_Bmat(theta_bspl_4int, inputs$basis_f,
    normalized = FALSE,
    center = inputs$center_th, scale = inputs$scale_th
  )
  dtheta_ode_4int <-
    do.call(
      cbind,
      lapply(1:p, function(j) {
        gam0[j] + (Bmat %*% as.vector(gam1[, , j]))
      })
    )
  out <- list()
  if (type %in% c("all", "data")) {
    fidelity.data <-
      neg_log_likelihood(
        ydata,
        do.call(rbind, replicate(R, theta_bspl, FALSE)), inputs$fam
      )
    out <- append(out, list(fidelity.data = fidelity.data))
  }
  if (type %in% c("all", "ode")) {
    fidelity.ode <- Metrics::mse(dtheta_ode_4int, dtheta_bspl_4int)
    out <- append(out, list(fidelity.ode = fidelity.ode))
  }
  return(out)
}


get.fidelity.vec <- function(Cmat, gam0, gam1, inputs,
                             type = c("all", "data", "ode")) {
  type <- match.arg(type)

  ydata <- inputs$ydata
  Hmat <- inputs$Hmat
  H4int <- inputs$H4int
  dH4int <- inputs$dH4int
  p <- dim(ydata)[2]
  nt <- dim(ydata)[1]
  R <- inputs$R

  theta_bspl <- Hmat %*% Cmat
  theta_bspl_4int <- H4int %*% Cmat
  dtheta_bspl_4int <- dH4int %*% Cmat
  Bmat <- form_Bmat(theta_bspl_4int, inputs$basis_f,
    normalized = FALSE,
    center = inputs$center_th, scale = inputs$scale_th
  )
  dtheta_ode_4int <-
    do.call(
      cbind,
      lapply(1:p, function(j) {
        gam0[j] + (Bmat %*% as.vector(gam1[, , j]))
      })
    )
  out <- list()
  if (type %in% c("all", "data")) {
    fidelity.data <- sapply(
      seq_len(p),
      function(j) {
        neg_log_likelihood(ydata[,j], theta_bspl[,j], inputs$fam)
      }
    )
  }
  if (type %in% c("all", "ode")) {
    fidelity.ode <- apply((dtheta_ode_4int - dtheta_bspl_4int)^2, 2, mean)
    out <- append(out, list(fidelity.ode = fidelity.ode))
  }
  return(out)
}


## sparse regression helper ------------
grpSparseReg <- function(X, y, group,
                         family = c("gaussian", "poisson", "binomial"),
                         sparse_penalty = c("grLasso", "grSCAD", "grMCP"),
                         lambda = NULL, intercept = TRUE,
                         pen.factor = NULL) {
  # group: integer sequence; 0 for un-penalized group
  # must be strictly increasing integers
  sparse_penalty <- match.arg(sparse_penalty)
  family <- match.arg(family)
  is_y_const <- (diff(range(y)) < 1e-8)
  if (!is_y_const) {
    if (is.null(lambda)) {
      fit <- grpreg::grpreg(
        X = X, y = y,
        group = group, penalty = sparse_penalty,
        family = family,
        group.multiplier = pen.factor
      )
    } else {
      fit <- grpreg::grpreg(
        X = X, y = y,
        group = group, penalty = sparse_penalty,
        family = family,
        lambda = lambda,
        group.multiplier = pen.factor
      )
    }
    lambdas <- fit$lambda
    pred <- predict(fit, X)
    b0 <- as.matrix(fit$beta)[1, ]
    beta <- as.matrix(fit$beta)[-1, , drop = FALSE]
  } else {
    lambdas <- lambda
    b0 <- rep(mean(y), length(lambdas))
    beta <- matrix(0, nrow = ncol(X), ncol = length(lambdas))
    pred <- matrix(b0, nrow = length(y), ncol = length(lambdas))
  }
  df <- apply(beta != 0, 2, sum)
  return(
    list(
      lambda = lambdas,
      pred = as.matrix(pred),
      b0 = b0,
      beta = beta,
      df = df
    )
  )
}

adjustGroup <- function(group) {
  # adjust the group argument to grpreg format
  idx <- sort(unique(group[group != 0]))
  group_new <- numeric(length(group))
  for (i in seq_along(idx)) {
    group_new[group == idx[i]] <- i
  }
  return(group_new)
}

# fix identifiability ----

set.F.meanzero <- function(gam, B, nbasis, nL, rm_id) {
  p <- ncol(B) / nL
  const <- numeric(p)
  for (k in 1:p) {
    subk <- ((k - 1) * nL + 1):(k * nL) # (nb,)
    mean.f_k <- colMeans(B[, subk] %*% gam[-rm_id, k, ]) # (p,)
    gam[, k, ] <- gam[, k, ] - t(replicate(nbasis, mean.f_k)) # nb x p
    const <- const + mean.f_k
  }
  return(list(gam = gam, const = const))
}


get.penalty.weights <- function(basis_f, gam1, Cmat, inputs) {
  p <- inputs$p
  rm_id <- inputs$rm_id
  w.mat <- matrix(NA, p, p)
  theta_mat <- inputs$H4int %*% Cmat
  for (j in 1:p) {
    norm_j <- numeric(p)
    for (k in 1:p) {
      if (length(rm_id) != 0) {
        gamjk <- numeric(length(inputs$rm_id) + length(gam1[, k, j]))
        gamjk[-rm_id] <- gam1[, k, j]
      } else { gamjk <- gam1[, k, j] }
      sigma_theta <- transformer(scale(theta_mat[, k], inputs$center_th[k], inputs$scale_th[k]))
      y <- as.vector(fda::eval.fd(sigma_theta, fda::fd(gamjk, basis_f[[k]])))
      norm_j[k] <- sqrt(mean((y - mean(y))^2))
    }
    w.mat[, j] <- norm_j^(-1)
  }
  return(w.mat)
}


get.ode_BIC <- function(
  nlambda, inputs, Cmat.all,
  gam1.all, gam0.all,
  L_lower_bound_factor = NULL
) {
  p <- inputs[["p"]]
  BICs <- sapply(
    1:nlambda,
    function(i) {
      theta_mat <- inputs[["H4int"]] %*% Cmat.all[, , i]
      dtheta_mat <- inputs[["dH4int"]] %*% Cmat.all[, , i]
      if (!is.null(L_lower_bound_factor)) {
        dtheta_dtrended <- pracma::detrend(dtheta_mat)
        L_lower_bound <- apply(dtheta_dtrended, 2, sd) * L_lower_bound_factor
      } else {
        L_lower_bound <- rep(0, p)
      }
      Bmat <- form_Bmat(
        theta_mat, inputs$basis_f,
        normalized = FALSE,
        center = inputs$center_th, scale = inputs$scale_th
      )
      dtheta_ode_4int <-
        do.call(
          cbind,
          lapply(1:p, function(j) {
            gam0.all[j, i] + (Bmat %*% as.vector(gam1.all[, , j, i]))
          })
        )
      L <- colMeans((dtheta_mat - dtheta_ode_4int)^2) |>
        pmax(L_lower_bound) |>
        log()
      bic <- L + apply(gam1.all[, , , i] != 0, 3, sum) * log(Nt) / Nt
      return(bic)
    }
  )
  return(BICs)
}

get.tol_idx <- function(idx, n_tol, idx_max) {
  out <- seq(
    idx - floor((n_tol - 1) / 2),
    idx + ceiling((n_tol - 1) / 2)
  )
  out[out < 1] <- 1
  out[out > idx_max] <- idx_max
  # out <- seq(idx, idx + n_tol - 1)
  # out[out < 1] <- 1
  # out[out > idx_max] <- idx_max
  return(out)
}

get.j_list <- function(p, random_order) {
  if (random_order) {
    j_list <- sample(1:p)
  } else {
    j_list <- 1:p
  }
  return(j_list)
}