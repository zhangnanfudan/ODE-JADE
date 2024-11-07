out <- generate.ydata(p = p, by = d_t, R = R, family = family,
                      dist_args = NULL,
                      plot.ode = FALSE, plot.legend = FALSE,
                      plot_ydata = FALSE)

shifts <- apply(out$truth$eta,2,function(x) (max(x)+min(x))/2)
scales <- apply(out$truth$eta,2,function(x) diff(range(x))) * 0.2

state_init <- (out$gen$eta[1,] - shifts)/scales

out <- generate.ydata(p = p, by = d_t, R = R, family = family,
                      dist_args = NULL, state_init = state_init,
                      a0p = out$gen$a0[7:10],
                      shifts = shifts, scales = scales,
                      plot.ode = FALSE, plot.legend = FALSE,
                      plot_ydata = FALSE)

truth <- out$truth
ydata <- out$ydata
times <- out$times
