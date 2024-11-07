res1 <-  bcd.ode.additive(ydata, times, family = family,
                          dist_args = dist_args,
                          lam_th_min = lam_th_min, lam_th_max = lam_th_max,
                          delta_lam_th_init = 1,
                          delta_update_method = "decrease",
                          lam_gam_list = lam_gam_list,
                          lam_gam_list0 = lam_gam_list0,
                          lam_gam_ratio = 1,
                          knots.theta.dt = 0.1,
                          scaling = 0.5,
                          norder_f = 4, nknots_f = 6,
                          knots_position = "quantile",
                          init_method = init_method,
                          lam_gam_init = lam_gam_init,
                          lam_gam_init2 = lam_gam_init2,
                          penalty = "grpLasso",
                          penalty.weighted = TRUE,
                          smooth_package = "stats",
                          theta_method = "1dH",
                          initw_method = initw_method,
                          MAXITER = MAXITER,
                          MAXITER0 = MAXITER0,
                          random_order = FALSE,
                          BICobj = "ode",
                          eps.conv = 1e-4,
                          ada.gam = 2,
                          know_truth = FALSE,
                          fit.trajectory = FALSE,
                          mode.test = FALSE,
                          verbose = FALSE,
                          lam.message = FALSE)

write(
  c("jointest",
    with(res1,
         perf.calculator(Cmat,gam0,gam1,inputs,truth,TRUE,
                         transformer = transformer, rm_id = inputs$rm_id)
    ) %>% unlist()) %>%
    paste(collapse = ","),
  filename, append = TRUE
)

# two-step
res2 <- twostep.additive(ydata, times, family = family,
                         methods = c("GRADE", "SAODE"), # "GRADE", "SAODE"
                         specs = list(
                           GRADE = list(
                             lam_gam_list = lam_gam_list.grade,
                             norder_f = 4, nknots_f = 4
                           ),
                           SAODE =list(
                             lam_gam_list = lam_gam_list.saode,
                             lam_gam_list2 = lam_gam_list2.saode,
                             norder_f = 4, nknots_f = 4
                           )
                         ),
                         plot_smooth = FALSE,
                         ada.gam = 2,
                         verbose = FALSE)


perf2 <- with(
  res2,
  perf.calculator.2stage(GRADE$basis_f, GRADE$gam1, GRADE$graph,
                         theta.d, dtheta.d,
                         theta.center, theta.scale,
                         GRADE$fidelity.data, truth,
                         subtract_mean = TRUE,
                         rm_id = 1)
)
write(
  paste(c("GRADE", unlist(perf2)), collapse = ","),
  filename, append = TRUE
)

perf2 <- with(
  res2,
  perf.calculator.2stage(SAODE$basis_f, SAODE$gam1, SAODE$graph,
                         theta.d, dtheta.d,
                         theta.center, theta.scale,
                         SAODE$fidelity.data, truth,
                         subtract_mean = TRUE,
                         rm_id = 1)
)
write(
  paste(c("SAODE", unlist(perf2)), collapse = ","),
  filename, append = TRUE
)