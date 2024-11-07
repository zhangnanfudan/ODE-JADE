library(dplyr)
library(tidyr)
library(ggplot2)
library(grid)
library(gridExtra)


# plot.fun_obj <- function(basis, gam1,
#                          transformer = function(x) {x}, 
#                          true_fun = NULL,
#                          range1=0, range2=1, by=0.005,
#                          normalize = FALSE,
#                          centers = 0, scales = 1,
#                          subtract_mean = FALSE) {
#   # transformer: the transformer function used before
#   # provided to the basis function
#   # centers, scales: length=p, scaling the theta_k, k=1,...,p
#   p <- dim(gam1)[2]
#   # old_par <- par(mfrow=c(p,p))
#   if (length(range1)==1) range1 <- rep(range1, p)
#   if (length(range2)==1) range2 <- rep(range2, p)
#   if (length(centers)==1) centers <- rep(centers, p)
#   if (length(scales)==1) scales <- rep(scales, p)
#   if (is.function(transformer)) {
#     transformer <- rep(list(transformer),p)
#   } else if (!is.list(transformer) || !all(sapply(transformer,is.function)) || length(transformer)!=p) {
#     stop("transformer should be 1 or list of p functions.")
#   }
#   plots <- list()
#   for (j in 1:p) {
#     for (k in 1:p) {
#       gamjk <- gam1[,k,j]
#       xx <- seq(range1[k],range2[k],by=by)
#       if (normalize) zz <- scale(xx,center=centers[k],scale=scales[k])
#       zz <- transformer[[k]](zz)
#       y_est <- as.vector(fda::eval.fd(zz,fda::fd(gamjk,basis)))
#       if (!is.null(true_fun)) y_true <- true_fun[[j]][[k]](xx)
#       if (subtract_mean) {
#         y_est <- y_est - mean(y_est)
#         y_true <- y_true - mean(y_true)
#       }
#       dat <- data.frame(x=xx, y_est=y_est, y_true=y_true)
#       dat <- tidyr::pivot_longer(dat, cols = c("y_est","y_true"),
#                                  names_to = "type", values_to = "val")
#       fig <- ggplot(dat, aes(x=x, y=val)) +
#         geom_line(aes(color=type, linetype=type), alpha=0.8, size=0.8) +
#         scale_color_manual(values=c("#095C6D","#ECAB15")) +
#         labs(x="",y="",title=paste0("f",j,k)) +
#         theme_light() + theme(legend.position = "none",
#                               title = element_text(size=10, face='bold'))
#       
#       plots <- append(plots, list(fig))
#     }
#   }
#   gridExtra::grid.arrange(grobs=plots)
# }



plot.fun_obj.jk <- function(j, k, basis, gam1,
                            transformer = function(x) {x}, 
                            true_fun = NULL,
                            range1=0, range2=1, by=0.01,
                            normalize = FALSE,
                            centers = 0, scales = 1,
                            subtract_mean = FALSE,
                            rm_id = c()) {
  p <- dim(gam1)[2]
  if (length(range1)==1) range1 <- rep(range1, p)
  if (length(range2)==1) range2 <- rep(range2, p)
  if (length(centers)==1) centers <- rep(centers, p)
  if (length(scales)==1) scales <- rep(scales, p)
  if (is.function(transformer)) {
    transformer <- rep(list(transformer),p)
  } else if (!is.list(transformer) || !all(sapply(transformer,is.function)) || length(transformer)!=p) {
    stop("transformer should be 1 or list of p functions.")
  }
  gamjk <- gam1[,k,j]
  if (length(gamjk)==basis[[k]]$nbasis-length(rm_id)) {
    tmp <- numeric(length(rm_id)+length(gamjk))
    if (length(rm_id)!=0) {tmp[-rm_id] <- gamjk} else {
      tmp <- gamjk
    }
    gamjk <- tmp
  } else {
    stop("Length error")
  }
  xx <- seq(range1[k],range2[k],by=by)
  if (normalize) zz <- scale(xx,center=centers[k],scale=scales[k])
  zz <- transformer[[k]](zz)
  y_est <- as.vector(fda::eval.fd(zz,fda::fd(gamjk,basis)))
  if (!is.null(true_fun)) y_true <- true_fun[[j]][[k]](xx)
  if (subtract_mean) {
    y_est <- y_est - mean(y_est)
    y_true <- y_true - mean(y_true)
  }
  dat <- data.frame(x=xx, y_est=y_est, y_true=y_true)
  dat <- tidyr::pivot_longer(dat, cols = c("y_est","y_true"),
                             names_to = "type", values_to = "val")
  fig <- ggplot(dat, aes(x=x, y=val)) +
    geom_line(aes(color=type, linetype=type), alpha=0.8, size=0.8) +
    scale_color_manual(values=c("#095C6D","#ECAB15")) +
    labs(x="",y="",title=paste0("f",j,k)) +
    theme_light() + theme(legend.position = "none",
                          title = element_text(size=10, face='bold'))
  print(fig)
}


plot.est.ode <- function(t1, theta = NULL, dtheta = NULL,
                         t2 = NULL, dtheta_ode = NULL, 
                         truth = NULL,
                         type = c("bspl", "ode")) {
  type <- match.arg(type)
  p <- ifelse(is.null(theta), ncol(dtheta), ncol(theta))
  n1 <- length(t1); if (!is.null(t2)) n2 <- length(t2)
  if (type=="bspl") {
    if (is.null(truth)) {
      dat <- cbind(data.frame(cls=rep("est",n1), times=t1), data.frame(theta))
    } else {
      dat <- rbind(
        cbind(data.frame(cls=rep("est",n1), times=t1), data.frame(theta)),
        cbind(data.frame(cls=rep("true",n1), times=t1), data.frame(truth$eta))
      )
    }
  } else if (type=="ode") {
    dat <- rbind(
      cbind(data.frame(cls = rep("bspl",n1), times=t1), data.frame(dtheta)),
      cbind(data.frame(cls = rep("ode",n2), times=t2), data.frame(dtheta_ode))
    )
  }
  dat <- tidyr::pivot_longer(dat, cols = paste0("X",1:p),
                             names_to = "j",
                             values_to = "val") %>% 
    mutate(j = factor(j, levels = paste0("X",1:p)))
  fig <- ggplot(dat) +
    geom_line(aes(x = times, y = val, col = cls, linetype = cls)) + 
    facet_wrap(vars(j), scales = "free") +
    theme_light()
  if (!is.null(truth)) {fig <- fig +
    scale_color_manual(values = c("chocolate3","darkgreen")) }
  print(fig)
}



plot.bcd.ode <- function(res, truth) {
  p <- res$inputs$p
  ## 1. bspl estimate
  t1 <- res1$inputs$t4int
  theta1 <- res1$inputs$H4int %*% res1$Cmat
  
  n1 <- length(t1); if (!is.null(t2)) n2 <- length(t2)
  if (type=="bspl") {
    if (is.null(truth)) {
      dat <- cbind(data.frame(cls=rep("est",n1), times=t1), data.frame(theta))
    } else {
      dat <- rbind(
        cbind(data.frame(cls=rep("est",n1), times=t1), data.frame(theta)),
        cbind(data.frame(cls=rep("true",n1), times=t1), data.frame(truth$eta))
      )
    }
  } else if (type=="ode") {
    dat <- rbind(
      cbind(data.frame(cls = rep("bspl",n1), times=t1), data.frame(dtheta)),
      cbind(data.frame(cls = rep("ode",n2), times=t2), data.frame(dtheta_ode))
    )
  }
  dat <- tidyr::pivot_longer(dat, cols = paste0("X",1:p),
                             names_to = "j",
                             values_to = "val") %>% 
    mutate(j = factor(j, levels = paste0("X",1:p)))
  fig <- ggplot(dat) +
    geom_line(aes(x = times, y = val, col = cls, linetype = cls)) + 
    facet_wrap(vars(j), scales = "free") +
    theme_light()
  if (!is.null(truth)) {fig <- fig +
    scale_color_manual(values = c("chocolate3","darkgreen")) }
  print(fig)
}



compare.curves <- function(t1, curves1, cap1 = "A",
                           t2, curves2, cap2 = "B",
                           j = NULL, title = "") {
  if (is.null(j)) {
    p <- ncol(curves1)
    if (ncol(curves2)!=p) stop("dim not matched")
    curvedat <- data.frame(rbind(curves1,curves2))
  } else {
    p <- 1
    curves1 <- curves1[,j]
    curves2 <- curves2[,j]
    curvedat <- data.frame(X = c(curves1, curves2))
  }
  dat <- cbind(
    data.frame(times = c(t1, t2),
               type = c(rep(cap1,length(t1)),rep(cap2,length(t2)))),
    curvedat
  )
  dat <- dat %>% 
    mutate(type = factor(type, levels = c(cap1, cap2)))
  if (is.null(j)) {
    dat <- dat %>% 
      pivot_longer(cols = paste0("X",1:p), names_to = "j", values_to = "val") %>% 
      mutate(j = factor(j, levels = paste0("X",1:p)))
    fig <- ggplot(dat) +
      geom_line(aes(x = times, y = val, col = type, linetype = type)) + 
      facet_wrap(vars(j), scales = "free") +
      theme_light() +
      scale_color_manual(values = c("chocolate3","darkgreen"))
    print(fig)
  } else {
    fig <- ggplot(dat) +
      geom_line(aes(x = times, y = X, col = type, linetype = type)) + 
      theme_light() +
      scale_color_manual(values = c("chocolate3","darkgreen"))
    print(fig)
  }
}


plot_fj.2stage <- function(basis, gam1, gam0, theta_est,
                           true_fun, a0, theta_true,
                           time_grid,
                           centers, scales,
                           Bc = 0) {
  source("./R/evaluate_helper.R", local = TRUE)
  p <- ncol(theta_est)
  nL <- dim(gam1)[1]
  ntheta <- nrow(theta_est)
  fmat_est <- matrix(0, nrow = ntheta, ncol = p)
  fmat_est <- sweep(fmat_est, 2, gam0, "+")
  fmat_true <- matrix(0, nrow = ntheta, ncol = p)
  fmat_true <- sweep(fmat_true, 2, a0, "+")
  if (length(Bc)==1) Bc <- rep(Bc, nL*p)
  
  for (j in 1:p) {
    for (k in 1:p) {
      fjk_est <-
        get.fjk_theta.2stage(basis[[k]], gam1[,k,j],
                             theta_est[,k], centers[k], scales[k],
                             subtract_mean = FALSE,
                             Bc_k = Bc[((k-1)*nL+1):(k*nL)],
                             rm_id = 1)
      fjk_true <-
        get.fjk_theta.truth(true_fun[[j]][[k]],
                            theta_true[,k],
                            subtract_mean = FALSE)
      fmat_est[,j] <- fmat_est[,j] + fjk_est$y
      fmat_true[,j] <- fmat_true[,j] + fjk_true$y
    }
  }
  
  dat <- cbind(
    type = rep(c("Est", "Truth"), each=ntheta),
    linewidth = rep(c("thick", "medium"), each=ntheta),
    t = rep(time_grid, 2),
    data.frame(
      rbind(fmat_est, fmat_true)
    )
  )
  
  dat <- dat %>% 
    pivot_longer(cols = paste0("X",1:p),
                 names_to = "j",
                 values_to = "f_theta_j") %>% 
    mutate(j = str_remove(j, "X")) %>% 
    mutate(j = factor(j, levels = as.character(1:p))) %>% 
    arrange(type, j, t)
  
  fig <- ggplot(dat) +
    geom_line(aes(x = t, y = f_theta_j, col = type, linetype = type,
                  size = linewidth)) + 
    ggpubr::theme_pubr() +
    labs(x = expression(t), y = "") + 
    scale_color_manual(values = c("#1BB438", "#625559")) +
    scale_linetype_manual(values=c("solid", "dashed")) +
    scale_size_manual(values=c(thick = 1, medium = 0.5)) + 
    theme(legend.title = element_blank(),
          legend.position = "none") + 
    guides(size = "none")+
    ggtitle(expression(f[j](theta)))+
    facet_wrap(vars(j), ncol = 5, nrow = 2, scales = "free")
  
  return(fig)
}

plot_fj.jade <- function(basis, gam1, gam0, theta_est,
                         true_fun, a0, theta_true,
                         time_grid,
                         centers, scales,
                         transformer = pracma::sigmoid,
                         rm_id = c()) {
  source("./R/evaluate_helper.R", local = TRUE)
  p <- ncol(theta_est)
  ntheta <- nrow(theta_est)
  fmat_est <- matrix(0, nrow = ntheta, ncol = p)
  fmat_est <- sweep(fmat_est, 2, gam0, "+")
  fmat_true <- matrix(0, nrow = ntheta, ncol = p)
  fmat_true <- sweep(fmat_true, 2, a0, "+")
  
  for (j in 1:p) {
    for (k in 1:p) {
      fjk_est <-
        get.fjk_theta.jade(basis[[k]], gam1[,k,j],
                           theta_est[,k], centers[k], scales[k],
                           subtract_mean = FALSE,
                           transformer = transformer,
                           rm_id = rm_id)
      fjk_true <-
        get.fjk_theta.truth(true_fun[[j]][[k]],
                            theta_true[,k],
                            subtract_mean = FALSE)
      fmat_est[,j] <- fmat_est[,j] + fjk_est$y
      fmat_true[,j] <- fmat_true[,j] + fjk_true$y
    }
  }
  
  dat <- cbind(
    type = rep(c("Est", "Truth"), each=ntheta),
    linewidth = rep(c("thick", "medium"), each=ntheta),
    t = rep(time_grid, 2),
    data.frame(
      rbind(fmat_est, fmat_true)
    )
  )
  
  dat <- dat %>% 
    pivot_longer(cols = paste0("X",1:p),
                 names_to = "j",
                 values_to = "f_theta_j") %>% 
    mutate(j = str_remove(j, "X")) %>% 
    mutate(j = factor(j, levels = as.character(1:p))) %>% 
    arrange(type, j, t)
  
  fig <- ggplot(dat) +
    geom_line(aes(x = t, y = f_theta_j, col = type, linetype = type,
                  size = linewidth)) + 
    ggpubr::theme_pubr() +
    labs(x = expression(t), y = "") + 
    scale_color_manual(values = c("#1BB438", "#625559")) +
    scale_linetype_manual(values=c("solid", "dashed")) +
    scale_size_manual(values=c(thick = 1, medium = 0.5)) + 
    theme(legend.title = element_blank(),
          legend.position = "none") + 
    guides(size = "none")+
    ggtitle(expression(f(theta[j])))+
    facet_wrap(vars(j), ncol = 5, nrow = 2, scales = "free")
  
  return(fig)
}
