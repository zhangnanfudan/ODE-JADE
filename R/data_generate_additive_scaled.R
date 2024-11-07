library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)

polyd <- function(x,d=3) return(sapply(1:d,function(i) x^i))
polydc <- function(x,d=3) return(c(1,polyd(x,d)))

coeftransform.j1k <- function(ajk, shift, scale) {
  if (length(ajk)!=3) stop("invalid ajk")
  newcoef <- c(
    sum(ajk * c(shift,shift^2,shift^3)),
    sum(scale * ajk * (1:3) * c(1,shift,shift^2)),
    sum(scale^2 * ajk[2:3] * c(1,3) * c(1,shift)),
    scale^3 * ajk[3]
  )
  return(newcoef)
}

coeftransform.j1 <- function(aj, aj0, j, shifts, scales) {
  p <- length(aj)/3
  if (length(shifts)!=p || length(scales)!=p) {
    stop("length unmatched")
  }
  bj <- numeric(length(aj)); bj0 <- aj0 
  for (k in 1:p) {
    if (all(aj[(3*k-2):(3*k)]==0)) next
    newcoef <- coeftransform.j1k(aj[(3*k-2):(3*k)], shifts[k], scales[k])
    bj0 <- bj0 + newcoef[1]
    bj[(3*k-2):(3*k)] <- newcoef[2:4]
  }
  return(list(bj = bj/scales[j], bj0 = bj0/scales[j]))
}

generate.ode <- function(p=10, t1=0, t2=20, by=0.5,
                         a0p = NULL, state_init = NULL,
                         plot.ode = TRUE, plot.legend = TRUE,
                         shifts = rep(0,p), scales = rep(1,p),
                         seed = NULL) {
  set.seed(seed)
  p0 <- 10; p1 <- 6; p2 <- p0-p1
  
  a10 = 0; a11 = c(1.2, 0.3, -0.6); a12 = c(0.1, 0.2, 0.2)
  a20 = 0.4; a21 = c(-2, 0, 0.4); a22 = c(0.5, 0.2, -0.3)
  
  a30 = -0.2; a33 = c(0, 0, 0); a34 = c(0.3, 0.4, 0.1)
  a40 = -0.2; a43 = c(0.2,-0.1,-0.2); a44 = c(0, 0, 0)
  
  a50 = 0.05; a55 = c(0, 0, 0); a56 = c(0.1, 0, -0.8)
  a60 = -0.05; a65 = c(0, 0, 0.5); a66 = c(0, 0, 0)
  
  a0 <- c(a10,a20,a30,a40,a50,a60)
  A1 <- rbind(c(a11,a12),c(a21,a22))
  A2 <- rbind(c(a33,a34),c(a43,a44))
  A3 <- rbind(c(a55,a56),c(a65,a66))
  a1 <- t(as.matrix(Matrix::bdiag(list(A1,A2,A3))))
  
  if (is.null(a0p)) {
    a0p <- runif(p2, min=-0.2, max=0.2)
  }
  
  ## transformation
  b1 <- a1 * 0; b0 <- a0 * 0
  for (j in 1:p1) {
    trans <- coeftransform.j1(a1[,j],a0[j],j,shifts[1:p1], scales[1:p1])
    b1[,j] <- trans$bj; b0[j] <- trans$bj0
  }
  b0p <- a0p / scales[(p1+1):p0]
  
  ode.model <- function(t, state, parameters) {
    b0 <- parameters[[1]]
    b1 <- parameters[[2]]
    state1 <- c(sapply(state[1:p1], polyd))
    dstate <- c(b0 + as.vector(state1 %*% b1), b0p)
    return(list(dstate))
  }
  if (!is.null(seed)) set.seed(seed)
  if (is.null(state_init)) {
    state1 <- c(-2,2,2,-2,-1.5,1.5)
    state2 <- rnorm(p2)
    state  <- (c(state1, state2) - shifts)/scales
  } else {
    if (length(state_init)!=p0) stop("Length not matched.")
    state <- state_init
  }
  set.seed(NULL)
  times <- seq(t1,t2,by)
  nt <- length(times)
  td <- seq(t1, t2, 0.01)
  
  get.d.eta <- function(eta) {
    eta1 <- eta[,1:p1]
    eta2 <- eta[,(p1+1):p0]
    poly3eta1 <- t(apply(eta1,1,function(x) {
      c(sapply(x, polyd))
    }))
    d.eta1 <- sweep(poly3eta1 %*% b1, 2, b0, '+')
    d.eta2 <- t(replicate(nrow(eta), b0p))
    d.eta <- cbind(d.eta1, d.eta2)
    return(d.eta)
  }
  
  eta <- deSolve::ode(state, times, ode.model, list(b0,b1))[,-1]
  d.eta <- get.d.eta(eta)
  etadense <- deSolve::ode(state, td, ode.model, list(b0,b1))[,-1]
  d.etadense <- get.d.eta(etadense)
  
  colnames(eta) <- colnames(d.eta) <- paste0("X",1:p0)
  
  a0 <- c(a0, a0p)
  b0 <- c(b0, b0p)
  
  # subseting
  eta <- eta[,1:p]
  d.eta <- d.eta[,1:p]
  if (p < 6) {
    a0 <- a0[1:p]
    a1 <- a1[1:(3*p),1:p]
  }
  
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
    d.eta = d.eta,
    a0 = b0,
    a1 = b1,
    a0star = a0,
    a1star = a1,
    td = td,
    etadense = etadense,
    d.etadense = d.etadense
  ))
}



generate.ydata.ode <- function(times, eta, R = 1,
                               family = c("gaussian", "poisson",
                                          "binomial", "nbinomial",
                                          "negbinomial"),
                               dist_args = NULL,
                               plot.ydata = TRUE, plot.legend = TRUE,
                               seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  family <- match.arg(family)
  nt <- nrow(eta)
  p <- ncol(eta)
  if (family=="gaussian") {
    SNR <- dist_args[['SNR']]
    if (!is.numeric(SNR)) stop("SNR not given")
    sdy <- apply(eta, 2, sd)/SNR
    ydata <- replicate(R,
                       eta + t(replicate(nt, rnorm(p, mean = 0, sd = sdy))),
                       simplify = FALSE)
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



generate.ydata <- function(p=10, t1=0, t2=20, by=0.5, R = 1,
                           family = c("gaussian", "poisson",
                                      "binomial", "nbinomial",
                                      "negbinomial"),
                           shifts = 0, scales = 1, state_init = NULL,
                           plot.ode = TRUE, plot.ydata = FALSE, plot.legend = FALSE,
                           dist_args = NULL, a0p = NULL,
                           ode.seed = NULL, data.seed = NULL,
                           ...) {
  p0 <- 10; p1 <- 6; p2 <- p0-p1
  if (length(shifts)==1) shifts <- rep(0,p)
  if (length(scales)==1) scales <- rep(1,p)
  
  gen <- generate.ode(p=p, t1=t1, t2=t2, by=by,
                      a0p = a0p, state_init = state_init,
                      plot.ode = plot.ode, plot.legend = plot.legend,
                      shifts = shifts, scales = scales,
                      seed = ode.seed)
  a0 <- gen$a0
  a1 <- gen$a1
  
  G1 <- rbind(c(1,1),c(1,1))
  G2 <- rbind(c(0,1),c(1,0))
  G3 <- rbind(c(0,1),c(1,0))
  G4 <- matrix(0,4,4)
  graph <- as.matrix(Matrix::bdiag(list(G1,G2,G3,G4)))
  graph <- matrix(as.logical(graph),nrow(graph),ncol(graph))
  
  graph <- graph[1:p,1:p]
  
  truth <- list(
    eta = gen$eta,
    d.eta = gen$d.eta,
    a0 = gen$a0,
    a1 = gen$a1,
    a0star = gen$a0star,
    a1star = gen$a1star,
    td = gen$td,
    etadense = gen$etadense,
    d.etadense = gen$d.etadense,
    graph = graph
  )
  
  true_fun <- rep(list(rep(list(function(x) {0*x}), p0)), p0)
  for (k in 1:(p1/2)) {
    true_fun[[2*k-1]][[2*k-1]] <- local({
      a11 <- a1[(6*k-5):(6*k-3),2*k-1]
      function(x) {sapply(x,function(xx){sum(a11*polyd(xx))})}
    })
    true_fun[[2*k-1]][[2*k]] <- local({
      a21 <- a1[(6*k-2):(6*k),2*k-1]
      function(x) {sapply(x,function(xx){sum(a21*polyd(xx))})}
    })
    true_fun[[2*k]][[2*k-1]] <- local({
      a12 <- a1[(6*k-5):(6*k-3),2*k]
      function(x) {sapply(x,function(xx){sum(a12*polyd(xx))})}
    })
    true_fun[[2*k]][[2*k]] <- local({
      a22 <- a1[(6*k-2):(6*k),2*k]
      function(x) {sapply(x,function(xx){sum(a22*polyd(xx))})}
    })
  }
  if (p < p0) true_fun <- lapply(true_fun[1:p], function(listobj) {listobj[1:p]})
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