out <- generate.ydata(p = p, by = d_t, R = R, family = family,
                      dist_args = dist_args,
                      plot.ode = FALSE, plot.legend = FALSE,
                      plot_ydata = FALSE)

truth <- out$truth
ydata <- out$ydata
times <- out$times