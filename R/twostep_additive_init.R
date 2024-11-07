twostep.additive.4jade <- function(
  ydata, times,
  family = c("gaussian", "poisson", "binomial"),
  methods = c("GRADE", "SAODE", "baseline"),
  dist_args = NULL,
  lam_gam_list = 10^seq(-4, 2, 0.1),
  lam_gam_list2 = 10^seq(-4, 4, 0.2),
  dt.dense = 0.01,
  dt.int = 0.002,
  scaling_factor = 1,
  norder_f = 4, nknots_f = 4,
  nbasis_reduce = 3,
  knots_position = c("equalspace", "quantile"),
  smooth_package = c("stats", "fda", "gss"),
  loglam_smooth = seq(-8, 6, 0.2),
  ada.gam = 1,
  verbose = TRUE,
  plot_smooth = FALSE
) {
  if (verbose) cat("Start two-step collocation:\n")
  ## 0. Preparation ----------------------
  ## arguments ##
  glmnet.family <- family <- match.arg(family)
  methods <- match.arg(methods)
  knots_position <- match.arg(knots_position)
  if (length(methods) == 0) stop("methods not valid")
  smooth_package <- match.arg(smooth_package)
  lam_gam_list <- sort(lam_gam_list, decreasing = TRUE)
  lam_gam_list2 <- sort(lam_gam_list2, decreasing = TRUE)

  ## functions ##
  source("R/utils/helper.R", local = TRUE)
  source("R/utils/exp_family.R", local = TRUE)
  source("R/utils/transformer.R", local = TRUE)
  source("R/utils/utils.R", local = TRUE)
  source("R/smoothing.R", local = TRUE)

  # other helper functions
  get.BIC <- function(nlambda, N, response, preds, dfs, fam = NULL) {
    if (is.null(fam) || fam$familyname == "gaussian") {
      BICs <- sapply(
        seq_len(nlambda),
        function(i) {
          MSE <- mean((response - preds[, i])^2)
          bic <- N * log(MSE) + dfs[i] * log(N)
          return(bic)
        }
      )
    } else if (fam$familyname %in% c("poisson", "binomial")) {
      BICs <- sapply(
        seq_len(nlambda),
        function(i) {
          L <- -mean(response * preds[, i] - fam$b(preds[, i]))
          bic <- 2 * L + dfs[i] * log(N) / N
          return(bic)
        }
      )
    } else {
      stop("NOT IMPLEMENTED")
    }
    return(BICs)
  }

  ## inputs ##
  R <- length(ydata)
  p <- dim(ydata[[1]])[2]
  nt <- dim(ydata[[1]])[1]
  ydata <- do.call(rbind, ydata)
  Nt <- nt * R
  # IMPORTANT: tvals is assumed to be equally spaced on time range
  if (nt != length(times)) stop("tvals and ydata not match")
  tvals <- rep(times, R) # tvals is a replicated version of times

  td <- seq(range(times)[1], range(times)[2], by = dt.dense)
  ntd <- length(td) # d for dense
  t4int <- seq(range(times)[1], range(times)[2], by = dt.int)
  nt4int <- length(t4int) # int for integral

  # family function
  fam <- get.family(family)

  ## create input object
  inputs <- list(
    p = p, nt = nt, Nt = Nt, ntd = ntd, nt4int = nt4int,
    R = R, ydata = ydata,
    times = times, tvals = tvals, td = td, t4int = t4int,
    dist_args = dist_args,
    b = fam$b, db = fam$db, d2b = fam$d2b,
    fam = fam
  )


  # %%%%%%%%%%%%%%%%%%%%%%%%%%%
  # 1. Smoothing ====

  smth_obj <- smoothing(family, smooth_package,
    inputs, loglam_smooth,
    by_basis = FALSE,
    eval_arg = c(times, td, t4int)
  )
  theta <- smth_obj$theta[1:nt, ]
  dtheta <- smth_obj$dtheta[1:nt, ]
  theta.d <- smth_obj$theta[(nt + 1):(nt + ntd), ]
  dtheta.d <- smth_obj$dtheta[(nt + 1):(nt + ntd), ]
  theta4int <- smth_obj$theta[(nt + ntd + 1):(nt + ntd + nt4int), ]
  dtheta4int <- smth_obj$dtheta[(nt + ntd + 1):(nt + ntd + nt4int), ]
  theta.center <- 0.5 * (apply(theta4int, 2, max) + apply(theta4int, 2, min))
  theta.scale <- scaling_factor * (apply(theta4int, 2, max) - apply(theta4int, 2, min))
  # TODO: tune the scaling parameter
  T_theta_hat <- normalize(theta4int, theta.center, theta.scale)

  const_id <- which(abs(theta.scale) < 1e-8)
  theta.scale[const_id] <- 1

  ## create basis ##
  nbasis_f <- nknots_f + norder_f - 2
  nL <- nbasis_f - nbasis_reduce
  if (nbasis_reduce == 0) {
    rm_id <- NULL
  } else {
    if (nbasis_reduce %% 2 == 0) {
      rm_id <- c(
        1:(nbasis_reduce / 2),
        (nbasis_f - nbasis_reduce / 2 + 1):nbasis_f
      )
    } else {
      rm_id <- c(
        1:((nbasis_reduce + 1) / 2),
        (nbasis_f - (nbasis_reduce - 1) / 2) + seq_len((nbasis_reduce - 1) / 2)
      )
    }
  }
  basis_f <- list()
  Bmat <- c()
  for (j in 1:p) {
    T_theta_j <- T_theta_hat[, j]
    if (knots_position == "quantile") {
      percents <- seq(0, 1, length.out = nknots_f - 2)
      knots_j <- quantile(T_theta_j, probs = percents)
      knots_j <- c(0, knots_j, 1)
      basis_f <- append(
        basis_f, list(
          fda::create.bspline.basis(
            rangeval = c(0, 1), breaks = knots_j,
            nbasis = nbasis_f, norder = norder_f
          )
        )
      )
    } else if (knots_position == "equalspace") {
      knots_j <- seq(0, 1, length.out = nknots_f)
      basis_f <- append(
        basis_f, list(
          fda::create.bspline.basis(
            rangeval = c(0, 1), breaks = knots_j,
            nbasis = nbasis_f, norder = norder_f
          )
        )
      )
    }
    Bmat <- cbind(Bmat, fda::eval.basis(T_theta_j, basis_f[[j]])[, -rm_id])
  }


  if (plot_smooth) {
    if (family == "gaussian") {
      dat <- cbind(
        data.frame(cls = rep("est", ntd), times = td),
        data.frame(theta.d)
      )
      dat <- tidyr::pivot_longer(dat,
        cols = paste0("X", 1:p),
        names_to = "j",
        values_to = "val"
      ) %>%
        mutate(j = factor(j, levels = paste0("X", 1:p)))
      fig <- ggplot(dat) +
        geom_line(aes(x = times, y = val, col = cls, linetype = cls)) +
        facet_wrap(vars(j), scales = "free") +
        theme_light()
      print(fig)
      # TODO: plot the observations
    }
  }

  out <- list(
    times = times,
    theta = theta, dtheta = dtheta,
    td = td,
    theta.d = theta.d, dtheta.d = dtheta.d,
    theta.center = theta.center,
    theta.scale = theta.scale,
    basis_f = basis_f,
    nbasis_f = nbasis_f,
    nL = nL, rm_id = rm_id
  )

  # FIXME: BIC not given for GRADE/SAODE intialization
  # %%%%%%%%%%%%%%%%%%%%%%%%%%%
  # 2. Fitting ODE components ====

  # 2.1 GRADE ------
  # if ("GRADE" %in% methods) {
  #   IntBmat <- apply(Bmat * diff(range(t4int)) / nt4int, 2, cumsum)
  #   IntBc <- apply(IntBmat, 2, mean)

  #   eid <- sapply(times, function(x) {
  #     ind <- which(abs(t4int - x) < 1e-6)
  #     ind <- ifelse(length(ind) == 0, NA, ind)
  #   })
  #   Z <- cbind(
  #     scale(tvals, TRUE, FALSE),
  #     repRbind(scale(IntBmat[eid, ], IntBc, FALSE), R)
  #   )
  #   group <- c(0, rep(1:p, each = nL))

  #   ### group Lasso
  #   gam0_est <- numeric(p)
  #   gam1_est <- array(dim = c(nL, p, p))
  #   gam0.all <- matrix(NA, p, length(lam_gam_list))
  #   gam1.all <- array(dim = c(nL, p, p, length(lam_gam_list)))
  #   best.lam <- numeric(p)
  #   best.lam_id <- numeric(p)
  #   for (j in 1:p) {
  #     if (glmnet.family == "gaussian") {
  #       fit <- gglasso::gglasso(
  #         Z, ydata[, j],
  #         group = group + 1,
  #         lambda = lam_gam_list, intercept = TRUE,
  #         pf = c(0, rep(sqrt(nL), p))
  #       )
  #       lambdas <- fit$lambda
  #       nlambda <- length(lambdas)
  #       preds <- predict(fit, Z)
  #       dfs <- fit$df
  #       gam0.all[j, 1:nlambda] <- fit$beta[1, ]
  #       gam1.all[, , j, 1:nlambda] <- fit$beta[-1, ]
  #     } else {
  #       fit <- grpSparseReg(
  #         Z, ydata[, j],
  #         group = group, family = glmnet.family,
  #         sparse_penalty = "grLasso",
  #         lambda = lam_gam_list, intercept = TRUE,
  #         pen.factor = rep(sqrt(nL), p)
  #       )
  #       lambdas <- fit$lambda
  #       nlambda <- length(lambdas)
  #       preds <- fit$pred
  #       dfs <- fit$df
  #       gam0.all[j, 1:nlambda] <- fit$beta[1, ]
  #       gam1.all[, , j, 1:nlambda] <- fit$beta[-1, ]
  #     }
  #     BICs <- get.BIC(nlambda, Nt, ydata[, j], preds, dfs, fam)
  #     minbic.id <- which.min(BICs)
  #     gam0_est[j] <- gam0.all[j, minbic.id]
  #     gam1_est[, , j] <- gam1.all[, , j, minbic.id]
  #     best.lam[j] <- lambdas[minbic.id]
  #     best.lam_id[j] <- minbic.id
  #   }

  #   out[['GRADE']] <- list(
  #       gam0 = gam0_est,
  #       gam1 = gam1_est,
  #       graph = apply(gam1_est != 0, c(2, 3), any),
  #       best.lam = best.lam,
  #       gam0.all = gam0.all,
  #       gam1.all = gam1.all,
  #       graph.all = apply(gam1.all[, , , ] != 0, c(2, 3, 4), any),
  #       lambdas = lam_gam_list,
  #       best.lam_id = best.lam_id,
  #       Bmat = Bmat,
  #       IntBc = IntBc
  #     )
  # }

  # # 2.1 SAODE ------
  # if ("SAODE" %in% methods) {

  #   Bc <- apply(Bmat, 2, mean)

  #   eid <- sapply(td, function(x) {
  #     ind <- which(abs(t4int - x) < 1e-6)
  #     ind <- ifelse(length(ind) == 0, NA, ind)
  #   })
  #   Z <- scale(Bmat[eid, ], Bc, FALSE)
  #   # currently we don't subtract Bc from Bmat
  #   # as it requires to shift gam1 after it is estimated
  #   # which may make the code a bit messy
  #   group <- rep(1:p, each = nL)

  #   # 2.2.i SAODE group Lasso ----
  #   gam0_grp <- numeric(p)
  #   gam1_grp <- array(dim = c(nL, p, p))
  #   gam0_grp.all <- matrix(NA, p, length(lam_gam_list))
  #   gam1_grp.all <- array(dim = c(nL, p, p, length(lam_gam_list)))
  #   best.lam_grp <- numeric(p)
  #   bicvalues_grp <- numeric(p)

  #   for (j in 1:p) {
  #     fit <- grpSparseReg(Z, dtheta.d[, j],
  #       group = group,
  #       lambda = lam_gam_list, intercept = TRUE,
  #       pen.factor = rep(sqrt(nL), p)
  #     )
  #     lambdas <- fit$lambda
  #     nlambda <- length(lambdas)
  #     preds <- fit$pred
  #     dfs <- fit$df
  #     gam0_grp.all[j, 1:nlambda] <- fit$b0
  #     gam1_grp.all[, , j, 1:nlambda] <- fit$beta
  #     BICs <- get.BIC(nlambda, Nt, dtheta.d[, j], preds, dfs)
  #     minbic.id <- which.min(BICs)
  #     gam0_grp[j] <- gam0_grp.all[j, minbic.id]
  #     gam1_grp[, , j] <- gam1_grp.all[, , j, minbic.id]
  #     best.lam_grp[j] <- lambdas[minbic.id]
  #     bicvalues_grp[j] <- BICs[minbic.id]
  #   }

  #   # 2.2.ii SAODE adaptive group Lasso ----

  #   # adaptive weights
  #   # FIXME: w.mat should be modified
  #   # as the Bmat here has not been standardized
  #   w.mat <- sqrt(apply(gam1_grp^2, 2:3, sum))^(-ada.gam)

  #   gam0_adagrp <- numeric(p)
  #   gam1_adagrp <- array(dim = c(nL, p, p))
  #   gam0_adagrp.all <- matrix(NA, p, length(lam_gam_list2))
  #   gam1_adagrp.all <- array(dim = c(nL, p, p, length(lam_gam_list2)))
  #   best.lam <- numeric(p)
  #   bicvalues <- numeric(p)
  #   for (j in 1:p) {
  #     fit <- grpSparseReg(Z, dtheta.d[, j],
  #       group = group,
  #       lambda = lam_gam_list2, intercept = TRUE,
  #       pen.factor = w.mat[, j]
  #     )
  #     lambdas <- fit$lambda
  #     nlambda <- length(lambdas)
  #     preds <- fit$pred
  #     dfs <- fit$df
  #     gam0_adagrp.all[j, 1:nlambda] <- fit$b0
  #     gam1_adagrp.all[, , j, 1:nlambda] <- fit$beta
  #     BICs <- get.BIC(nlambda, Nt, dtheta.d[, j], preds, dfs)
  #     minbic.id <- which.min(BICs)
  #     gam0_adagrp[j] <- gam0_adagrp.all[j, minbic.id]
  #     gam1_adagrp[, , j] <- gam1_adagrp.all[, , j, minbic.id]
  #     best.lam[j] <- lambdas[minbic.id]
  #     bicvalues[j] <- BICs[minbic.id]
  #   }

  #   # 2.2.iii SAODE individual adaptive lasso ----

  #   # adaptive weights
  #   w.mat <- abs(structure(gam1_adagrp, dim = c(nL * p, p)))^(-ada.gam)

  #   gam0_indv <- numeric(p)
  #   gam1_indv <- array(dim = c(nL, p, p))
  #   gam0_indv.all <- matrix(NA, p, length(lam_gam_list2))
  #   gam1_indv.all <- array(dim = c(nL, p, p, length(lam_gam_list2)))
  #   best.lam <- numeric(p)
  #   bicvalues <- numeric(p)

  #   for (j in 1:p) {
  #     idx <- !is.infinite(w.mat[, j])
  #     if (sum(idx) > 0) {
  #       fit <- glmnet::glmnet(
  #         Z[, idx], dtheta.d[, j],
  #         lambda = lam_gam_list2,
  #         penalty.factor = c(w.mat[idx, j])
  #       )
  #       lambdas <- fit$lambda
  #       nlambda <- length(lambdas)
  #       preds <- predict(fit, Z[, idx])
  #       dfs <- fit$df
  #       gam0_indv.all[j, 1:nlambda] <- fit$a0
  #       gam1_indv.all[, , j, 1:nlambda][idx] <- as.matrix(fit$beta)
  #       gam1_indv.all[, , j, 1:nlambda][!idx] <- 0
  #       BICs <- get.BIC(nlambda, Nt, dtheta.d[, j], preds, dfs)
  #       minbic.id <- which.min(BICs)
  #       gam0_indv[j] <- gam0_indv.all[j, minbic.id]
  #       gam1_indv[, , j] <- gam1_indv.all[, , j, minbic.id]
  #       best.lam[j] <- lambdas[minbic.id]
  #       bicvalues[j] <- BICs[minbic.id]
  #     } else {
  #       nlambda <- length(lam_gam_list2)
  #       gam0_indv.all[j, 1:nlambda] <- gam0_indv[j] <- mean(dtheta.d[, j])
  #       gam1_indv.all[, , j, 1:nlambda] <- gam1_indv[, , j] <- 0
  #       best.lam[j] <- Inf
  #       bicvalues[j] <- NA
  #     }
  #   }

  #   out[["SAODE"]] <- list(
  #     gam0 = gam0_indv,
  #     gam1 = gam1_indv,
  #     graph = apply(gam1_indv != 0, c(2, 3), any),
  #     best.lam = best.lam,
  #     bicvalues = bicvalues,
  #     gam0.all = gam0_indv.all,
  #     gam1.all = gam1_indv.all,
  #     lambdas = lam_gam_list2,
  #     Bmat = Bmat,
  #     Bc = Bc
  #   )
  # }


  # 2.3 baseline ------
  if ("baseline" %in% methods) {

    eid <- sapply(td, function(x) {
      ind <- which(abs(t4int - x) < 1e-6)
      ind <- ifelse(length(ind) == 0, NA, ind)
    })
    Z <- Bmat[eid, ]
    group <- rep(1:p, each = nL)

    gam0_est <- numeric(p)
    gam1_est <- array(dim = c(nL, p, p))
    gam0.all <- matrix(NA, p, length(lam_gam_list))
    gam1.all <- array(dim = c(nL, p, p, length(lam_gam_list)))
    best.lam <- numeric(p)
    best.lam_id <- numeric(p)
    bicvalues <- numeric(p)
    for (j in 1:p) {
      fit <- grpSparseReg(
        Z, dtheta.d[, j],
        group = group, sparse_penalty = "grLasso",
        lambda = lam_gam_list, intercept = TRUE,
        pen.factor = rep(sqrt(nL), p)
      )
      lambdas <- fit$lambda
      nlambda <- length(lambdas)
      preds <- fit$pred
      dfs <- fit$df
      gam0.all[j, 1:nlambda] <- fit$b0
      gam1.all[, , j, 1:nlambda] <- fit$beta
      BICs <- get.BIC(nlambda, Nt, dtheta.d[, j], preds, dfs)
      minbic.id <- which.min(BICs)
      gam0_est[j] <- gam0.all[j, minbic.id]
      gam1_est[, , j] <- gam1.all[, , j, minbic.id]
      best.lam[j] <- lambdas[minbic.id]
      best.lam_id[j] <- minbic.id
      bicvalues[j] <- BICs[minbic.id]
    }

    out$baseline <- append(
      out$baseline,
      list(
        gam0 = gam0_est,
        gam1 = gam1_est,
        graph = apply(gam1_est != 0, c(2, 3), any),
        best.lam_id = best.lam_id,
        best.lam = best.lam,
        BICvalue = bicvalues,
        gam0.all = gam0.all,
        gam1.all = gam1.all,
        graph.all = apply(gam1.all[, , , ] != 0, c(2, 3, 4), any),
        lambda = lam_gam_list,
        Bmat = Bmat
      )
    )
  }

  return(out)
}