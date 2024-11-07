## implemented for a tolerant range in tuning

bcd.ode.additive <- function(
    ydata, times,
    family = c("gaussian", "poisson", "binomial"),
    dist_args = NULL,
    lam_th_min = 10^(-3), lam_th_max = 10^(2),
    delta_lam_th_init = 0.5,
    delta_update_method = c("const", "decrease"),
    lam_gam_list = 10^(seq(-3, 2, 0.5)),
    lam_gam_ratio = c("const", "auto"),
    delta_lam_gam = NULL,
    norder_th = 6, nknots_th = NULL,
    knots_position = c("equalspace", "quantile"),
    knots.theta.dt = 0.1,
    dt.int = 0.01,
    norder_f = 4, nknots_f = 4,
    nbasis_reduce = 3,
    scaling_factor = 1,
    init_method = c("baseline", "GRADE", "SAODE"),
    lam_gam_init = 10^seq(-2, 2, 0.1),
    lam_gam_init2 = 10^seq(-1, 4, 0.2),
    theta_method = c("1gd", "1dH"),
    penalty = c("grpLasso", "grpSCAD", "grpMCP"),
    penalty.weighted = TRUE,
    initialize_penalty_weights = FALSE,
    initw_method = c("twostage", "jade"), # deprecated
    smooth_package = c("stats", "fda", "gss"),
    loglam_smooth = seq(-8, 6, 0.2),
    MAXITER = 10,
    NUM_REWEIGHT = 1,
    final_result_only = TRUE,
    random_order = TRUE,
    tuning_tolerant = 3,
    BICobj = c("data", "ode"),
    stepsize = 1, amj.alpha = 0.1, amj.beta = 0.5,
    eps.conv = 1e-4,
    eps.th1 = 1e-3, eps.th0 = 0.8,
    ada.gam = 1,
    L_lower_bound_factor = NULL,
    verbose = TRUE,
    plot_smooth = FALSE,
    lam.message = TRUE,
    conv.message = FALSE,
    know_truth = FALSE,
    fit.trajectory = FALSE,
    mode.test = FALSE) {
  # %%%%%%%%%%%%%%%%%%%%%%%%%%%
  # 0 Preparation ====
  ## arguments ##
  if (verbose) cat("Start block coordinate descent:\n")
  delta_update_method <- match.arg(delta_update_method)
  penalty <- match.arg(penalty)
  init_method <- match.arg(init_method)
  knots_position <- match.arg(knots_position)
  sparse_penalty <- grpreg.penalty <-
    c("grLasso", "grSCAD", "grMCP", "grBridge")[
      match(penalty, c("grpLasso", "grpSCAD", "grpMCP", "grpBridge"))
    ]
  smooth_package <- match.arg(smooth_package)
  theta_method <- match.arg(theta_method)
  family <- match.arg(family)
  if (family != "gaussian") smooth_package <- "gss"
  lam_gam_list <- sort(lam_gam_list, decreasing = TRUE)
  lam_gam_init <- sort(lam_gam_init, decreasing = TRUE)
  lam_gam_init2 <- sort(lam_gam_init2, decreasing = TRUE)
  lam_gam_ratio <- match.arg(lam_gam_ratio)
  if (is.null(delta_lam_gam)) {
    delta_lam_gam <- min(10^abs(diff(log10(lam_gam_list) / 2)))
  }
  if (delta_lam_gam <= 0) stop("wrong delta_lam_gam")
  BICobj <- match.arg(BICobj)
  initw_method <- match.arg(initw_method)
  if (!penalty.weighted) {
    initialize_penalty_weights <- FALSE
    NUM_REWEIGHT <- 1
  }

  ## functions ##
  source("R/twostep_additive_init.R", local = TRUE)
  source("R/utils/utils.R", local = TRUE)
  source("R/utils/transformer.R", local = TRUE)
  source("R/utils/exp_family.R", local = TRUE)
  source("R/optimize/optimize_additive_bcd.R", local = TRUE)
  source("R/smoothing.R", local = TRUE)
  if (know_truth) source("R/evaluate.R", local = TRUE)

  # %%%%%%%%%%%%%%%%%%%%%%%%%%%
  # 1 Smoothing ====
  res_init <- twostep.additive.4jade(
    ydata, times,
    family = family,
    dt.dense = dt.int,
    methods = init_method,
    lam_gam_list = lam_gam_init,
    lam_gam_list2 = lam_gam_init2,
    scaling_factor = scaling_factor,
    norder_f = norder_f, nknots_f = nknots_f,
    nbasis_reduce = nbasis_reduce,
    knots_position = knots_position,
    smooth_package = smooth_package,
    loglam_smooth = loglam_smooth,
    ada.gam = ada.gam,
    verbose = FALSE,
    plot_smooth = FALSE
  )

  # %%%%%%%%%%%%%%%%%%%%%%%%%%%
  # 2 Initialize parameters and constants ====

  ## inputs ##
  R <- length(ydata)
  p <- dim(ydata[[1]])[2]
  nt <- dim(ydata[[1]])[1]
  ydata <- do.call(rbind, ydata)
  Nt <- nt * R
  # IMPORTANT: tvals is assumed to be equally spaced on time range
  # otherwise there will be problems on bspline basis' knots
  if (nt != length(times)) stop("tvals and ydata not match")
  tvals <- rep(times, R) # tvals is a replicated version of times

  t4int <- seq(range(times)[1], range(times)[2], by = dt.int)
  nt4int <- length(t4int)
  theta4int <- res_init$theta.d
  if (nrow(theta4int) != nt4int) stop("nt4int unmatched")

  # quadrature weights
  quadwts <- get.quadwts(t4int)
  # family function
  fam <- get.family(family)

  ## create basis ##
  trange <- range(times)
  knots <- seq(trange[1], trange[2], by = knots.theta.dt)
  basis_th <- fda::create.bspline.basis(rangeval = trange, breaks = knots, norder = norder_th)
  Hmat <- fda::eval.basis(times, basis_th)
  dHmat <- fda::eval.basis(times, basis_th, Lfdobj = 1)
  H4int <- fda::eval.basis(t4int, basis_th)
  dH4int <- fda::eval.basis(t4int, basis_th, Lfdobj = 1)

  basis_f <- res_init$basis_f
  nbasis_f <- res_init$nbasis_f
  nL <- res_init$nL
  rm_id <- res_init$rm_id

  ## create input object and control object
  inputs <- list(
    p = p, nt = nt, Nt = Nt, nt4int = nt4int,
    R = R, ydata = ydata,
    times = times, tvals = tvals, t4int = t4int,
    quadwts = quadwts,
    dist_args = dist_args,
    basis_th = basis_th, nbasis_th = basis_th$nbasis,
    basis_f = basis_f, nbasis_f = nbasis_f, nL = nL, rm_id = rm_id,
    Hmat = Hmat, dHmat = dHmat, H4int = H4int, dH4int = dH4int,
    b = fam$b, db = fam$db, d2b = fam$d2b,
    fam = fam,
    group = rep(1:p, each = nL)
  )

  con <- list(
    MAXITER = MAXITER,
    verbose = verbose,
    WARNING = FALSE,
    conv.message = conv.message,
    theta_method = theta_method,
    sparse_penalty = sparse_penalty,
    eps.conv = eps.conv,
    stepsize = stepsize,
    amj.alpha = amj.alpha,
    amj.beta = amj.beta,
    ada.gam = ada.gam
  )

  inputs[["center_th"]] <- center_th <- res_init$theta.center
  inputs[["scale_th"]] <- scale_th <- res_init$theta.scale

  # turn smoothing curve into basis expansion
  Cmat_init <- fda::smooth.basis(t4int, theta4int, fdParobj = basis_th)$fd$coefs

  gam1_array <- res_init[[init_method]]$gam1
  gam0_vec <- res_init[[init_method]]$gam0

  lam_id <- sapply( # n_tol x p
    1:p,
    function(j) {
      get.tol_idx(
        idx = res_init[[init_method]]$best.lam_id[j],
        n_tol = tuning_tolerant,
        idx_max = sum(apply(!is.na(res_init[[init_method]]$gam1.all[, , j, ]), 3, all))
      )
    }
  ) |> (\(x) matrix(x, ncol = p))()

  # consider a range of initializations
  param_init <- with(
    res_init[[init_method]],
    lapply(
      seq_len(tuning_tolerant),
      function(i) {
        indexes <- cbind(1:p, lam_id[i, ])
        out <- list(
          gam1 = sapply(1:p, \(ii) gam1.all[, , indexes[ii, 1], indexes[ii, 2]]) |>
            structure(dim = c(nL, p, p)),
          gam0 = gam0.all[indexes],
          best.lam_gam_init = lambda[lam_id[i, ]]
        )
        # the ratio of lam_gam across different processes
        out[["lam_gam_ratio"]] <- switch(lam_gam_ratio,
          const = rep(1, p),
          auto = out$best.lam_gam_init / out$best.lam_gam_init[1]
        )
        return(out)
      }
    )
  )

  # %%%%%%%%%%%%%%%%%%%%%%%%%%%
  # 3 Fitting ====
  # if (fit.trajectory) {
  #   nrecord <- MAXITER * NUM_REWEIGHT + 10
  #   GHvalues <- matrix(NA, nrow = nrecord, ncol = length(lam_gam_list))
  #   lam_th_values <- matrix(NA, nrow = nrecord, ncol = length(lam_gam_list))
  # }
  # if keep track on past estimations
  if (!final_result_only) past_estimations <- list()

  # repeated BCD iterations
  for (i_rw in seq_len(NUM_REWEIGHT)) {
    if (verbose) cat("Reweighting:", i_rw, "of", NUM_REWEIGHT, "\n")
    estimation <- list()
    for (i_init in seq_len(tuning_tolerant)) {
      if (verbose) cat("Reinitializing:", i_init, "of", tuning_tolerant, "\n")
      # empty estimations
      nlambda <- length(lam_gam_list)
      Cmat.all <- array(dim = c(basis_th$nbasis, p, nlambda))
      gam0.all <- array(dim = c(p, nlambda))
      gam1.all <- array(dim = c(nL, p, p, nlambda))
      lam_gam_ratio <- param_init[[i_init]]$lam_gam_ratio
      # stop decreasing lam_gam if err_flag==TRUE
      err_flag <- FALSE
      # block coordinate descent
      for (i_gam in seq_along(lam_gam_list)) {
        lam_th <- lam_th_min
        lam_gam <- lam_gam_list[i_gam] * lam_gam_ratio
        if (lam.message) cat("::lambda_gamma=", lam_gam[1], "\n")

        # initializations
        if (i_rw > 1) {
          Cmat_est <- param_init[[i_init]]$Cmat
        } else {
          Cmat_est <- Cmat_init
        }
        gam0_est <- param_init[[i_init]]$gam0
        gam1_est <- param_init[[i_init]]$gam1

        # initial weight matrix
        if (penalty.weighted && (initialize_penalty_weights || i_rw > 1)) {
          inputs[["w.mat"]] <- w.mat <-
            get.penalty.weights(basis_f, gam1_est, Cmat_est, inputs)^ada.gam
        } else {
          w.mat <- matrix(1, nrow = p, ncol = p)
        }
        inputs[["w.mat"]] <- w.mat

        for (it in seq_len(con$MAXITER)) {
          j_list <- get.j_list(p, random_order)
          # update theta_j ----
          for (j in j_list) {
            Cmat_est[, j] <- estimate.theta(
              j, gam0_est, gam1_est, Cmat_est, inputs, lam_th, con
            )
          }

          j_list <- get.j_list(p, random_order)
          # update gamma_j ----
          for (j in j_list) {
            tryCatch(
              {
                gamj <- estimate.gamma(
                  j, gam0_est, gam1_est, Cmat_est, inputs,
                  lam_th, lam_gam[j], con
                )
              },
              error = function(e) {
                err_flag <<- TRUE
              }
            )
            if (err_flag) break

            gam0_est[j] <- gamj[1]
            gam1_est[, , j] <- gamj[-1]
          }
        }
        if (err_flag) break
        if (conv.message && it == con$MAXITER) {
          warning("Maximum number of iterations has been reached.")
        }
        Cmat.all[, , i_gam] <- Cmat_est
        gam0.all[, i_gam] <- gam0_est
        gam1.all[, , , i_gam] <- gam1_est
      }

      nlambda <- ifelse(err_flag, i_gam - 1, i_gam)
      lambdas <- lam_gam_list[1:nlambda]

      # store the temporal results
      estimation[[i_init]] <- list(
        nlambda = nlambda,
        lambdas = lambdas,
        Cmat.all = Cmat.all,
        gam0.all = gam0.all,
        gam1.all = gam1.all,
        w.mat = w.mat,
        lam_gam_ratio = lam_gam_ratio # TODO: adapt lam_gam_ratio
      )
    }
    ## BIC select of results ----
    BICs <- lapply(
      seq_len(tuning_tolerant),
      \(i) with(
        estimation[[i]],
        get.ode_BIC(
          nlambda, inputs, Cmat.all,
          gam1.all, gam0.all,
          L_lower_bound_factor
        ) |> colSums()
      )
    )
    max_len <- max(sapply(BICs, length))
    BICs <- do.call(
      cbind,
      lapply(BICs, \(x) c(x, rep(Inf, max_len - length(x))))
    )
    min_BIC <- min(BICs)
    index <- which(BICs == min_BIC, arr.ind = TRUE)
    if (nrow(index) > 1) index <- index[1, , drop = FALSE]
    # select for the next round of iterations
    lam_id <- get.tol_idx(
      index[1], tuning_tolerant, estimation[[index[2]]]$nlambda
    )
    # adapt lan_gam_ratio ----
    lam_gam_ratio <- estimation[[index[2]]]$lam_gam_ratio
    BIC_by_dim <- with(
      estimation[[index[2]]],
      get.ode_BIC(
        nlambda, inputs, Cmat.all,
        gam1.all, gam0.all,
        L_lower_bound_factor
      )
    )
    index_by_dim <- apply(BIC_by_dim, 1, which.min)
    # make some lam_gam_ratio larger
    lam_gam_ratio[index_by_dim < index[1]] <-
      lam_gam_ratio[index_by_dim < index[1]] * delta_lam_gam # TODO: set an upper/lower limit
    # make some lam_gam_ratio smaller
    lam_gam_ratio[index_by_dim > index[1]] <-
      lam_gam_ratio[index_by_dim > index[1]] / delta_lam_gam # TODO: consider when index hits boundary
    # TODO: consider to standardize w.mat and lam_gam_list
    # TODO: such that the best lam_gam_list is in the middle of the lam_gam_list
    # pass on the results ----
    if (i_rw < NUM_REWEIGHT) {
      for (ii in seq_along(lam_id)) {
        param_init[[ii]]$Cmat <-
          estimation[[index[2]]]$Cmat.all[, , lam_id[ii]]
        param_init[[ii]]$gam0 <-
          estimation[[index[2]]]$gam0.all[, lam_id[ii]]
        param_init[[ii]]$gam1 <-
          estimation[[index[2]]]$gam1.all[, , , lam_id[ii]]
        param_init[[ii]]$lam_gam_ratio <- lam_gam_ratio
      }
      if (!final_result_only) {
        past_estimations[[i_rw]] <- list(
          Cmat = estimation[[index[2]]]$Cmat.all[, , index[1]],
          gam0 = estimation[[index[2]]]$gam0.all[, index[1]],
          gam1 = estimation[[index[2]]]$gam1.all[, , , index[1]],
          BICvalue = BICs[index],
          best.lam = estimation[[index[2]]]$lambdas[index[1]],
          lam_gam_ratio = lam_gam_ratio
        )
      }
    } else {
      Cmat <- estimation[[index[2]]]$Cmat.all[, , index[1]]
      gam0 <- estimation[[index[2]]]$gam0.all[, index[1]]
      gam1 <- estimation[[index[2]]]$gam1.all[, , , index[1]]
      nlambda <- estimation[[index[2]]]$nlambda
      lambdas <- estimation[[index[2]]]$lambdas
      best.lam <- lambdas[index[1]]
      w.mat <- estimation[[index[2]]]$w.mat
      Cmat.all <- estimation[[index[2]]]$Cmat.all
      gam0.all <- estimation[[index[2]]]$gam0.all
      gam1.all <- estimation[[index[2]]]$gam1.all
      BICvalue <- BICs[index]
    }
  }

  ## evaluate the final result ----
  fidelity <- get.fidelity(Cmat, gam0, gam1, inputs, type = "all")

  out <- list(
    Cmat = Cmat, gam0 = gam0, gam1 = gam1,
    graph = apply(gam1 != 0, 2:3, any),
    nlambda = nlambda, lambdas = lambdas,
    best.lam = best.lam,
    lam_gam_ratio = lam_gam_ratio,
    inits = list(Cmat = Cmat_init, gam0 = gam0_vec, gam1 = gam1_array),
    inputs = inputs,
    w.mat = w.mat, best.lam = best.lam,
    Cmat.all = Cmat.all,
    gam0.all = gam0.all, gam1.all = gam1.all,
    graph.all = apply(gam1.all != 0, 2:4, any),
    fidelity = fidelity,
    BICvalue = BICvalue
  )

  if (!final_result_only) {
    out <- append(out, list(past_estimations = past_estimations))
  }

  return(out)
}
