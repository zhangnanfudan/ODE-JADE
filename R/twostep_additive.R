# a grid search routine

twostep.additive <- function(
  ydata, times,
  family = c("gaussian", "poisson", "binomial"),
  methods = c("GRADE", "SAODE"),
  dist_args = NULL,
  dt.dense = 0.01,
  dt.int = 0.002,
  specs = list(
    GRADE = list(
      lam_gam_list = 10^seq(-6, 0, 0.2),
      lam_gam_list2 = 10^seq(-6, 4, 0.2),
      norder_f = 4, nknots_f = 4
    ),
    SAODE = list(
      lam_gam_list = 10^seq(-6, 0, 0.2),
      lam_gam_list2 = 10^seq(-6, 4, 0.2),
      norder_f = 4, nknots_f = 4
    )
  ),
  smooth_package = c("stats", "fda", "gss"),
  all_grpreg = FALSE,
  loglam_smooth = seq(-8, 6, 0.2),
  ada.gam = 1,
  verbose = TRUE,
  plot_smooth = FALSE,
  mode.test = FALSE
) {
  if (verbose) cat("Start two-step collocation:\n")
  ## 0. Preparation ----------------------
  ## arguments ##
  glmnet.family <- family <- match.arg(family)
  methods <- base::intersect(methods, c("GRADE", "SAODE"))
  if (length(methods) == 0) stop("methods not valid")
  smooth_package <- match.arg(smooth_package)
  # TODO: get default settings for specs

  ## functions ##
  source("R/utils/helper.R", local = TRUE)
  source("R/utils/exp_family.R", local = TRUE)
  source("R/smoothing.R", local = TRUE)

  # other helper functions
  neg_log_likelihood <- function(ydata, theta_est, fam) {
    b <- fam$b
    # negative log-likelihood for exp.family
    theta_est <- as.matrix(theta_est)
    ydata <- as.matrix(ydata)
    nt <- dim(theta_est)[1]
    p <- dim(theta_est)[2]
    if (dim(ydata)[1] != nt || dim(ydata)[2] != p) {
      stop("preds, y, df not matched")
    }
    btheta_est <- b(theta_est)
    return(-1 / (nt * p) * sum(ydata * theta_est - btheta_est)) # or -1/(nt)
  }

  grpSparseReg <- function(
    X, y, group,
    family = c("gaussian", "poisson", "binomial"),
    lambda = NULL, intercept = TRUE,
    pen.factor = NULL
  ) {
    # group: integer sequence; 0 for un-penalized group
    # must be strictly increasing integers
    family <- match.arg(family)
    grpsize <- table(group) # may go wrong for bad "group"
    if (is.null(pen.factor)) pen.factor <- sqrt(grpsize)
    is_y_const <- (diff(range(y)) < 1e-8)
    if (!is_y_const) {
      if (is.null(lambda)) {
        fit <- grpreg::grpreg(
          X = X, y = y,
          group = group, penalty = "grLasso",
          family = family,
          group.multiplier = pen.factor
        )
      } else {
        fit <- grpreg::grpreg(
          X = X, y = y,
          group = group, penalty = "grLasso",
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

  adjustGroup <- function(group) {
    # adjust the group argument to grpreg format
    idx <- sort(unique(group[group != 0]))
    group_new <- numeric(length(group))
    for (i in seq_along(idx)) {
      group_new[group == idx[i]] <- i
    }
    return(group_new)
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
  theta.scale <- 0.501 * (apply(theta4int, 2, max) - apply(theta4int, 2, min))
  theta_hat <- scale(theta4int, theta.center, theta.scale) # scaled processes
  # larger scale than 0.5 to prevent overflow

  const_id <- which(abs(theta.scale) < 1e-8)
  theta.scale[const_id] <- 1

  if (plot_smooth) {
    if (family == "gaussian") {
      dat <- cbind(
        data.frame(cls = rep("est", ntd), times = td),
        data.frame(theta.d)
      ) %>%
        tidyr::pivot_longer(
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
    theta.scale = theta.scale
  )


  # %%%%%%%%%%%%%%%%%%%%%%%%%%%
  # 2. Fitting ODE components ====

  # 2.1 GRADE ------
  if ("GRADE" %in% methods) {

    ## create basis ##
    nbasis_f <- specs$GRADE$nknots_f + specs$GRADE$norder_f - 2
    nL <- nbasis_f - 1 # number of used basis
    basis_f <- list()
    Bmat <- c()
    for (j in 1:p) {
      theta_j <- theta_hat[, j]
      percents <- seq(0, 1, length.out = specs$GRADE$nknots_f)
      knots_j <- quantile(theta_j, probs = percents)
      knots_j <- c(-1, knots_j[-c(1, length(knots_j))], 1)
      basis_f <- append(
        basis_f, list(
          fda::create.bspline.basis(
            rangeval = c(-1, 1), breaks = knots_j,
            nbasis = nbasis_f, norder = specs$GRADE$norder_f
          )
        )
      )
      Bmat <- cbind(Bmat, fda::eval.basis(theta_j, basis_f[[j]])[, -1])
    }
    Bc <- apply(Bmat, 2, mean)
    IntBmat <- apply(Bmat * diff(range(t4int)) / nt4int, 2, cumsum)
    IntBc <- apply(IntBmat, 2, mean)

    if (verbose) {
      cat("=========== GRADE ===========", "\n")
    }
    # re-scale theta to the support of b-spline basis
    # IMPORTANT: evaluate b-spline on a dense sequence
    eid <- sapply(times, function(x) {
      ind <- which(abs(t4int - x) < 1e-6)
      ind <- ifelse(length(ind) == 0, NA, ind)
    })
    Z <- cbind(
      scale(tvals, TRUE, FALSE),
      repRbind(scale(IntBmat[eid, ], IntBc, FALSE), R)
    )
    group <- c(0, rep(1:p, each = nL))

    # group Lasso
    lam_gam_list <- sort(specs$GRADE$lam_gam_list, decreasing = TRUE)
    gam0_est <- numeric(p)
    gam1_est <- array(dim = c(nL, p, p))
    gam0.all <- matrix(NA, p, length(lam_gam_list))
    gam1.all <- array(dim = c(nL, p, p, length(lam_gam_list)))
    best.lam <- numeric(p)
    bicvalues <- numeric(p)
    for (j in 1:p) {
      if (glmnet.family == "gaussian" && !all_grpreg) {
        fit <- gglasso::gglasso(
          Z, ydata[, j],
          group = group + 1,
          lambda = lam_gam_list, intercept = TRUE,
          pf = c(0, rep(sqrt(nL), p))
        )
        lambdas <- fit$lambda
        nlambda <- length(lambdas)
        preds <- predict(fit, Z)
        dfs <- fit$df
        gam0.all[j, 1:nlambda] <- fit$beta[1, ]
        gam1.all[, , j, 1:nlambda] <- fit$beta[-1, ]
      } else {
        fit <- grpSparseReg(
          Z, ydata[, j],
          group = group, family = glmnet.family,
          lambda = lam_gam_list, intercept = TRUE,
          pen.factor = rep(sqrt(nL), p)
        )
        lambdas <- fit$lambda
        nlambda <- length(lambdas)
        preds <- fit$pred
        dfs <- fit$df
        gam0.all[j, 1:nlambda] <- fit$beta[1, ]
        gam1.all[, , j, 1:nlambda] <- fit$beta[-1, ]
      }
      BICs <- get.BIC(nlambda, Nt, ydata[, j], preds, dfs, fam)
      minbic.id <- which.min(BICs)
      gam0_est[j] <- gam0.all[j, minbic.id]
      gam1_est[, , j] <- gam1.all[, , j, minbic.id]
      best.lam[j] <- lambdas[minbic.id]
      bicvalues[j] <- BICs[minbic.id]
    }

    # evaluate fidelity
    dtheta_ode_4int <- sapply(1:p, function(j) {
      gam0_est[j] + (Bmat %*% as.vector(gam1_est[, , j]))
    })
    fidelity.ode <- mean((dtheta_ode_4int - dtheta4int)^2)
    fidelity.data <- neg_log_likelihood(ydata, repRbind(theta, R), fam)

    out[["GRADE"]] <- list(
      gam0 = gam0_est,
      gam1 = gam1_est,
      graph = apply(gam1_est != 0, c(2, 3), any),
      best.lam = best.lam,
      bicvalues = bicvalues,
      fidelity.ode = fidelity.ode,
      fidelity.data = fidelity.data,
      gam0.all = gam0.all,
      gam1.all = gam1.all,
      graph.all = apply(gam1.all[, , , ] != 0, c(2, 3, 4), any),
      lambdas = lam_gam_list,
      basis_f = basis_f,
      Bmat = Bmat,
      IntBc = IntBc
    )

    if (verbose) {
      cat(paste("Data fidelity - likelihood:", fidelity.data), "\n")
      cat(paste("ODE fidelity", fidelity.ode), "\n")
    }
  }

  # 2.2 SAODE ------
  if ("SAODE" %in% methods) {
    if (verbose) {
      cat("=========== SAODE ===========", "\n")
    } 
    ## create basis ##
    nbasis_f <- specs$SAODE$nknots_f + specs$SAODE$norder_f - 2
    nL <- nbasis_f - 1
    basis_f <- list()
    Bmat <- c()
    for (j in 1:p) {
      theta_j <- theta_hat[, j]
      percents <- seq(0, 1, length.out = specs$SAODE$nknots_f)
      knots_j <- quantile(theta_j, probs = percents)
      knots_j <- c(-1, knots_j[-c(1, length(knots_j))], 1)
      basis_f <- append(
        basis_f, list(
          fda::create.bspline.basis(
            rangeval = c(-1, 1), breaks = knots_j,
            nbasis = nbasis_f, norder = specs$SAODE$norder_f
          )
        )
      )
      Bmat <- cbind(Bmat, fda::eval.basis(theta_j, basis_f[[j]])[, -1])
    }
    Bc <- apply(Bmat, 2, mean)

    eid <- sapply(td, function(x) {
      ind <- which(abs(t4int - x) < 1e-6)
      ind <- ifelse(length(ind) == 0, NA, ind)
    })
    Z <- scale(Bmat[eid, ], Bc, FALSE)
    group <- rep(1:p, each = nL)

    lam_gam_list <- sort(specs$SAODE$lam_gam_list, decreasing = TRUE)
    lam_gam_list2 <- sort(specs$SAODE$lam_gam_list2, decreasing = TRUE)

    # 2.2.i SAODE group Lasso ----
    gam0_grp <- numeric(p)
    gam1_grp <- array(dim = c(nL, p, p))
    gam0_grp.all <- matrix(NA, p, length(lam_gam_list))
    gam1_grp.all <- array(dim = c(nL, p, p, length(lam_gam_list)))
    best.lam_grp <- numeric(p)
    bicvalues_grp <- numeric(p)

    for (j in 1:p) {
      fit <- grpSparseReg(Z, dtheta.d[, j],
        group = group,
        lambda = lam_gam_list, intercept = TRUE,
        pen.factor = rep(sqrt(nL), p)
      )
      lambdas <- fit$lambda
      nlambda <- length(lambdas)
      preds <- fit$pred
      dfs <- fit$df
      gam0_grp.all[j, 1:nlambda] <- fit$b0
      gam1_grp.all[, , j, 1:nlambda] <- fit$beta
      BICs <- get.BIC(nlambda, Nt, dtheta.d[, j], preds, dfs)
      minbic.id <- which.min(BICs)
      gam0_grp[j] <- gam0_grp.all[j, minbic.id]
      gam1_grp[, , j] <- gam1_grp.all[, , j, minbic.id]
      best.lam_grp[j] <- lambdas[minbic.id]
      bicvalues_grp[j] <- BICs[minbic.id]
    }

    # 2.2.ii SAODE adaptive group Lasso ----

    # adaptive weights
    w.mat <- sqrt(apply(gam1_grp^2, 2:3, sum))^(-ada.gam)

    gam0_adagrp <- numeric(p)
    gam1_adagrp <- array(dim = c(nL, p, p))
    gam0_adagrp.all <- matrix(NA, p, length(lam_gam_list2))
    gam1_adagrp.all <- array(dim = c(nL, p, p, length(lam_gam_list2)))
    best.lam <- numeric(p)
    bicvalues <- numeric(p)
    for (j in 1:p) {
      fit <- grpSparseReg(Z, dtheta.d[, j],
        group = group,
        lambda = lam_gam_list2, intercept = TRUE,
        pen.factor = w.mat[, j]
      )
      lambdas <- fit$lambda
      nlambda <- length(lambdas)
      preds <- fit$pred
      dfs <- fit$df
      gam0_adagrp.all[j, 1:nlambda] <- fit$b0
      gam1_adagrp.all[, , j, 1:nlambda] <- fit$beta
      BICs <- get.BIC(nlambda, Nt, dtheta.d[, j], preds, dfs)
      minbic.id <- which.min(BICs)
      gam0_adagrp[j] <- gam0_adagrp.all[j, minbic.id]
      gam1_adagrp[, , j] <- gam1_adagrp.all[, , j, minbic.id]
      best.lam[j] <- lambdas[minbic.id]
      bicvalues[j] <- BICs[minbic.id]
    }

    # 2.2.iii SAODE individual adaptive lasso ----

    # adaptive weights
    w.mat <- abs(structure(gam1_adagrp, dim = c(nL * p, p)))^(-ada.gam)

    gam0_indv <- numeric(p)
    gam1_indv <- array(dim = c(nL, p, p))
    gam0_indv.all <- matrix(NA, p, length(lam_gam_list2))
    gam1_indv.all <- array(dim = c(nL, p, p, length(lam_gam_list2)))
    best.lam <- numeric(p)
    bicvalues <- numeric(p)

    for (j in 1:p) {
      idx <- !is.infinite(w.mat[, j])
      if (sum(idx) > 0) {
        fit <- glmnet::glmnet(
          Z[, idx], dtheta.d[, j],
          lambda = lam_gam_list2,
          penalty.factor = c(w.mat[idx, j])
        )
        lambdas <- fit$lambda
        nlambda <- length(lambdas)
        preds <- predict(fit, Z[, idx])
        dfs <- fit$df
        gam0_indv.all[j, 1:nlambda] <- fit$a0
        gam1_indv.all[, , j, 1:nlambda][idx] <- as.matrix(fit$beta)
        gam1_indv.all[, , j, 1:nlambda][!idx] <- 0
        BICs <- get.BIC(nlambda, Nt, dtheta.d[, j], preds, dfs)
        minbic.id <- which.min(BICs)
        gam0_indv[j] <- gam0_indv.all[j, minbic.id]
        gam1_indv[, , j] <- gam1_indv.all[, , j, minbic.id]
        best.lam[j] <- lambdas[minbic.id]
        bicvalues[j] <- BICs[minbic.id]
      } else {
        nlambda <- length(lam_gam_list2)
        gam0_indv.all[j, 1:nlambda] <- gam0_indv[j] <- mean(dtheta.d[, j])
        gam1_indv.all[, , j, 1:nlambda] <- gam1_indv[, , j] <- 0
        best.lam[j] <- Inf
        bicvalues[j] <- get.BIC(
          1, Nt, dtheta.d[, j],
          matrix(gam0_indv[j], nrow = ntd, ncol = 1), 0
        )
      }
    }

    # evaluate fidelity
    dtheta_ode_4int <- sapply(1:p, function(j) {
      gam0_indv[j] + (scale(Bmat, Bc, F) %*% as.vector(gam1_indv[, , j]))
    })
    fidelity.ode <- mean((dtheta_ode_4int - dtheta4int)^2)
    fidelity.data <- neg_log_likelihood(ydata, repRbind(theta, R), fam)

    out[['SAODE']] <- list(
        gam0 = gam0_indv,
        gam1 = gam1_indv,
        graph = apply(gam1_indv != 0, c(2, 3), any),
        best.lam = best.lam,
        bicvalues = bicvalues,
        fidelity.ode = fidelity.ode,
        fidelity.data = fidelity.data,
        gam0.all = gam0_indv.all,
        gam1.all = gam1_indv.all,
        graph.all = apply(gam1_indv.all != 0, 2:4, any),

        gam0_grp = gam0_grp,
        gam1_grp = gam1_grp,
        best.lam_grp = best.lam_grp,
        bicvalues_grp = bicvalues_grp,
        graph_grp = apply(gam1_grp != 0, c(2, 3), any),
        gam0_grp.all = gam0_grp.all,
        gam1_grp.all = gam1_grp.all,
        graph_grp.all = apply(gam1_grp.all != 0, 2:4, any),

        gam0_adagrp = gam0_adagrp,
        gam1_adagrp = gam1_adagrp,
        graph_adagrp = apply(gam1_adagrp != 0, c(2, 3), any),
        gam0_adagrp.all = gam0_adagrp.all,
        gam1_adagrp.all = gam1_adagrp.all,
        graph_adagrp.all = apply(gam1_adagrp.all != 0, 2:4, any),

        lambdas = lam_gam_list2,
        lambdas_grp = lam_gam_list,

        basis_f = basis_f,
        Bmat = Bmat,
        Bc = Bc
      )

    if (verbose) {
      cat(paste("Data fidelity - likelihood:", fidelity.data), "\n")
      cat(paste("ODE fidelity", fidelity.ode), "\n")
    }
  }

  return(out)
}