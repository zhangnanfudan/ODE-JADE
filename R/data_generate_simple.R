library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)


generate.ode <- function(p, a0, a1, state = NULL,
                         t1=0, t2=3, by=0.025,
                         state_sigma = 1,
                         plot.ode = TRUE, plot.legend = TRUE,
                         seed = NULL) {
  
  ode.model <- function(t, state, parameters) {
    a0 <- parameters[[1]]
    a1 <- parameters[[2]]
    dstate <- a0 + as.vector(state %*% a1)
    return(list(dstate))
  }
  if (!is.null(seed)) set.seed(seed)
  if (is.null(state)) {
    y0 = rnorm(p/2, mean=0, sd=state_sigma)
    state <- numeric(p)
    for (k in 1:(p/2)) {
      state[2*k-1] <- sin(y0[k])
      state[2*k] <- cos(y0[k])
    }
  }
  set.seed(NULL)
  times <- seq(t1,t2,by)
  nt <- length(times)
  
  eta <- deSolve::ode(state, times, ode.model, list(a0,a1))[,-1]
  d.eta <- sweep(eta %*% a1, 2, a0, '+')
  
  if (plot.ode==TRUE) {
    plot.eta <- data.frame(times=times, eta) %>%
      pivot_longer(cols = paste0("X",1:p),
                   names_to = "process",
                   values_to = "value") %>%
      ggplot(aes(x = times, y = value, color = process)) +
      geom_line() + ggtitle("state value") +
      labs(x = "", y = "") +
      theme_bw()
    plot.deta <- data.frame(times=times, d.eta) %>%
      pivot_longer(cols = paste0("X",1:p),
                   names_to = "process",
                   values_to = "value") %>%
      ggplot(aes(x = times, y = value, color = process)) +
      geom_line() + ggtitle("derivative value") +
      labs(x = "", y = "") +
      theme_bw()
    if (plot.legend==FALSE) {
      plot.eta <- plot.eta + theme(legend.position = "none")
      plot.deta <- plot.deta + theme(legend.position = "none")
    }
    grid.arrange(arrangeGrob(plot.eta, plot.deta, ncol = 1))
  }

  return(list(
    nt = nt,
    times = times,
    eta = eta,
    d.eta = d.eta
  ))
}


generate.ydata.ode <- function(times, eta, R = 1,
                               family = c("gaussian", "poisson",
                                          "binomial", "nbinomial",
                                          "negbinomial"),
                               dist_args = NULL,
                               plot.ydata = FALSE, plot.legend = TRUE,
                               seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  family <- match.arg(family)
  nt <- nrow(eta)
  p <- ncol(eta)
  if (family=="gaussian") {
    SNR <- dist_args[['SNR']]
    if (is.infinite(SNR)) {
      ydata <- replicate(R, eta,
                         simplify = FALSE)
    } else {
      sdy <- apply(eta, 2, sd)/SNR
      ydata <- replicate(R,
                         eta + t(replicate(nt, rnorm(p, mean = 0, sd = sdy))),
                         simplify = FALSE)
    }
  } else if (family=="poisson") {
    lam <- exp(eta)
    ydata <- replicate(R,
                       matrix(
                         rpois(length(lam), lambda = lam),
                         nrow = nt, ncol = p
                       ),
                       simplify = FALSE
    )
  } else if (family=="binomial") {
    prob <- 1/(1+exp(-eta))
    ydata <- replicate(R,
                       matrix(
                         rbinom(length(prob), size = 1, prob = prob),
                         nrow = nt, ncol = p
                       ),
                       simplify = FALSE
    )
  } else if (family=="nbinomial") {
    n <- dist_args[['n']]
    if (!is.numeric(n)) stop('n not given')
    prob <- 1/(1+exp(-eta))
    ydata <- replicate(R,
                       matrix(
                         rbinom(length(prob), size = n, prob = prob),
                         nrow = nt, ncol = p
                       ),
                       simplify = FALSE
    )
  } else if (family=="negbinomial") {
    nu <- dist_args[['nu']]
    if (!is.numeric(nu)) stop('n not given')
    mu <- exp(eta)
    prob <- nu/(mu+nu)
    ydata <- replicate(R,
                       matrix(
                         rnbinom(length(mu), nu, prob),
                         nrow = nt, ncol = p
                       ),
                       simplify = FALSE
    )
  }
  if (plot.ydata) {
    fig <- data.frame(times = times, ydata[[1]]) %>%
      pivot_longer(cols = paste0("X",1:p), names_to = "X", values_to = "value") %>%
      ggplot(aes(x = times, y = value, color=X)) +
      geom_point() + ggtitle("observations") + labs(x = "") +
      theme_bw()
    if (plot.legend==FALSE) {
      fig <- fig + theme(legend.position = "none")
    }
    print(fig)
  }
  set.seed(NULL)
  
  return(ydata)
}


generate.ydata <- function(p, state = NULL,
                           t1=0, t2=3, by=0.025, R = 1,
                           family = c("gaussian", "poisson",
                                      "binomial", "nbinomial",
                                      "negbinomial"),
                           plot.ode = TRUE, plot.ydata = FALSE, plot.legend = FALSE,
                           dist_args = NULL,
                           ode.seed = NULL, data.seed = NULL,
                           ...) {
  
  if (p %% 2 != 0) stop("This experiment needs an even p")
  
  a0 <- rep(0,p)
  a1 <- matrix(0,p,p)
  for (k in 1:(p/2)) {
    a1[2*k-1,2*k] <- 2*k*pi
    a1[2*k,2*k-1] <- -2*k*pi
  }
  
  gen <- generate.ode(p, a0, a1, state = state,
                      t1=t1, t2=t2, by=by,
                      plot.ode = plot.ode, plot.legend = plot.legend,
                      seed = ode.seed)
  
  truth <- list(
    a0 = a0,
    a1 = a1,
    graph = a1!=0,
    eta = gen$eta,
    d.eta = gen$d.eta
  )
  
  true_fun <- rep(list(rep(list(function(x) {0*x}), p)), p)
  for (k in 1:(p/2)) {
    true_fun[[2*k-1]][[2*k]] <- local({
      a12 <- a1[2*k,2*k-1]
      function(x) {a12*x}
    })
    true_fun[[2*k]][[2*k-1]] <- local({
      a21 <- a1[2*k-1,2*k]
      function(x) {a21*x}
    })
  }
  truth$true_fun <- true_fun

  ydata <- generate.ydata.ode(gen$times, gen$eta, family=family, R=R,
                              plot.ydata = plot.ydata,
                              plot.legend = plot.legend,
                              dist_args = dist_args,
                              seed = data.seed)
  
  return(list(
    gen = gen,
    truth = truth,
    ydata = ydata,
    times = gen$times
  ))
}