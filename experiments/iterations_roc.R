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

lambdas <- res1$lambdas
for (i in seq_along(lambdas)) {
  write(
    c("jointest",
      with(res1,
           perf.calculator(Cmat.all[,,i],gam0.all[i],gam1.all[,,,i],inputs,truth,TRUE,
                           transformer = transformer, rm_id = inputs$rm_id)
      ) %>% unlist(),
      lambdas[i]) %>%
      paste(collapse = ","),
    filename, append = TRUE
  )
}



# two-step
res2 <- twostep.additive(ydata, times, family = family,
                         methods = c("GRADE", "baseline"), # "GRADE", "SAODE"
                         specs = list(
                           GRADE = list(
                             lam_gam_list = lam_gam_list.grade,
                             norder_f = 4, nknots_f = 4
                           ),
                           baseline =list(
                             lam_gam_list = lam_gam_list.saode,
                             norder_f = 4, nknots_f = 4
                           )
                         ),
                         plot_smooth = FALSE,
                         ada.gam = 2,
                         verbose = FALSE)

nlambda <- max(which(apply(!is.na(res2$GRADE$graph.all), 3, all)))
lambdas <- lam_gam_list.grade[1:nlambda]
for (i in seq_along(lambdas)) {
  perf2 <- with(
    res2,
    perf.calculator.2stage(GRADE$basis_f, GRADE$gam1.all[,,,i], GRADE$graph.all[,,i],
                           theta.d, dtheta.d,
                           theta.center, theta.scale,
                           GRADE$fidelity.data, truth,
                           subtract_mean = TRUE,
                           rm_id = 1)
  )
  write(
    paste(c("GRADE", unlist(perf2), lambdas[i]), collapse = ","),
    filename, append = TRUE
  )
}

nlambda <- max(which(apply(!is.na(res2$baseline$graph.all), 3, all)))
lambdas <- lam_gam_list.saode[1:nlambda]
for (i in seq_along(lambdas)) {
  perf2 <- with(
    res2,
    perf.calculator.2stage(baseline$basis_f, baseline$gam1.all[,,,i], baseline$graph.all[,,i],
                           theta.d, dtheta.d,
                           theta.center, theta.scale,
                           baseline$fidelity.data, truth,
                           subtract_mean = TRUE,
                           rm_id = 1)
  )
  write(
    paste(c("baseline", unlist(perf2), lambdas[i]), collapse = ","),
    filename, append = TRUE
  )
}

