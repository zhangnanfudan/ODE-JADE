### JADE rw=1 -------------------------
res1_list <- list()
for (i in 1:num_rep) {
  res1_list[[i]] <- bcd.ode.additive(
    ydata, times,
    family = family,
    dist_args = dist_args,
    lam_th_min = lam_th_min, lam_th_max = lam_th_max,
    delta_lam_th_init = 1,
    delta_update_method = "decrease",
    lam_gam_list = lam_gam_list,
    lam_gam_ratio = "auto",
    delta_lam_gam = 10^0.2,
    knots.theta.dt = 0.1,
    scaling_factor = 1,
    norder_f = 4, nknots_f = 6,
    knots_position = "quantile",
    init_method = init_method,
    lam_gam_init = lam_gam_init,
    lam_gam_init2 = lam_gam_init2,
    penalty = "grpLasso",
    penalty.weighted = TRUE,
    smooth_package = "stats",
    theta_method = "1dH",
    initialize_penalty_weights = initw,
    MAXITER = MAXITER,
    NUM_REWEIGHT = NUM_REWEIGHT,
    final_result_only = FALSE,
    random_order = TRUE,
    tuning_tolerant = n_tol,
    BICobj = "ode",
    eps.conv = 1e-4,
    ada.gam = 1,
    L_lower_bound_factor = 1e-8,
    know_truth = FALSE,
    fit.trajectory = FALSE,
    mode.test = FALSE,
    verbose = FALSE,
    lam.message = FALSE
  )
}

BICs <- sapply(res1_list, \(x) x$BICvalue)
minbic.id <- which.min(BICs)
res1 <- res1_list[[minbic.id]]


for (ii in seq_len(NUM_REWEIGHT - 1)) {
  write(
    c(
      paste0("jointest_rw", ii),
      with(
        res1$past_estimations[[ii]],
        perf.calculator(
          Cmat, gam0, gam1, res1$inputs, truth, TRUE,
          transformer = transformer, rm_id = res1$inputs$rm_id
        )
      ) %>% unlist()
    ) %>%
      paste(collapse = ","),
    filename,
    append = TRUE
  )
}

write(
  c(
    paste0("jointest_rw", NUM_REWEIGHT),
    with(
      res1,
      perf.calculator(
        Cmat, gam0, gam1, inputs, truth, TRUE,
        transformer = transformer, rm_id = inputs$rm_id
      )
    ) %>% unlist()
  ) %>%
    paste(collapse = ","),
  filename,
  append = TRUE
)


### SA-ODE/GRADE -------------------------

res2_list <- list()
for (i in 1:num_rep) {
  res2_list[[i]] <- twostep.additive(
    ydata, times,
    family = family,
    methods = c("GRADE", "SAODE"), # "GRADE", "SAODE"
    specs = list(
      GRADE = list(
        lam_gam_list = lam_gam_list.grade,
        lam_gam_list2 = lam_gam_list.grade,
        norder_f = 4, nknots_f = 4
      ),
      SAODE = list(
        lam_gam_list = lam_gam_list.saode,
        lam_gam_list2 = lam_gam_list2.saode,
        norder_f = 4, nknots_f = 4
      )
    ),
    plot_smooth = FALSE,
    all_grpreg = all_grpreg,
    ada.gam = 2,
    verbose = FALSE
  )
}

# GRADE
BICs <- sapply(res2_list, \(x) sum(x$GRADE$bicvalues))
minbic.id <- which.min(BICs)
res2.GRADE <- res2_list[[minbic.id]]
# SAODE
BICs <- sapply(res2_list, \(x) sum(x$SAODE$bicvalues))
minbic.id <- which.min(BICs)
res2.SAODE <- res2_list[[minbic.id]]


perf2 <- with(
  res2.GRADE,
  perf.calculator.2stage(
    GRADE$basis_f, GRADE$gam1, GRADE$gam0, GRADE$graph,
    theta.d, dtheta.d,
    theta.center, theta.scale,
    GRADE$fidelity.data, truth,
    subtract_mean = TRUE,
    rm_id = 1
  )
)
write(
  paste(c("GRADE", unlist(perf2)), collapse = ","),
  filename, append = TRUE
)

perf2 <- with(
  res2.SAODE,
  perf.calculator.2stage(
    SAODE$basis_f, SAODE$gam1, SAODE$gam0, SAODE$graph,
    theta.d, dtheta.d,
    theta.center, theta.scale,
    SAODE$fidelity.data, truth,
    subtract_mean = TRUE,
    Bc = SAODE$Bc,
    rm_id = 1
  )
)
write(
  paste(c("SAODE", unlist(perf2)), collapse = ","),
  filename, append = TRUE
)