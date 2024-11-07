library(dplyr)
library(tidyr)
library(ggplot2)
library(grid)
library(gridExtra)


### est of additive components -----------

f.relerr <- function(basis, gam1,
                     transformer = function(x) {x}, 
                     true_fun = NULL,
                     range1=0, range2=1, by=0.001,
                     normalize = FALSE,
                     centers = 0, scales = 1,
                     subtract_mean = FALSE,
                     rm_id = c()) {
  # transformer: the transformer function used before
  # provided to the basis function
  # centers, scales: length=p, scaling the theta_k, k=1,...,p
  p <- dim(gam1)[2]
  # old_par <- par(mfrow=c(p,p))
  if (length(range1)==1) range1 <- rep(range1, p)
  if (length(range2)==1) range2 <- rep(range2, p)
  if (length(centers)==1) centers <- rep(centers, p)
  if (length(scales)==1) scales <- rep(scales, p)
  if (is.function(transformer)) {
    transformer <- rep(list(transformer),p)
  } else if (!is.list(transformer) || !all(sapply(transformer,is.function)) || length(transformer)!=p) {
    stop("transformer should be 1 or list of p functions.")
  }
  err_mat <- matrix(NA, nrow = p, ncol = p)
  for (j in 1:p) {
    for (k in 1:p) {
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
      y_est <- as.vector(fda::eval.fd(zz,fda::fd(gamjk,basis[[k]])))
      if (!is.null(true_fun)) y_true <- true_fun[[j]][[k]](xx)
      if (subtract_mean) {
        y_est <- y_est - mean(y_est)
        y_true <- y_true - mean(y_true)
      }
      err_mat[j,k] <- mean((y_est-y_true)^2)
    }
  }
  return(err_mat)
}


ftheta.relerr <- function(basis, gam1,
                          transformer = function(x) {x}, 
                          true_fun = NULL,
                          xmat1 = NULL, xmat2 = NULL,
                          normalize = FALSE,
                          centers = 0, scales = 1,
                          subtract_mean = FALSE,
                          rm_id = c()) {
  ## calculate the error matrix along a given sequence
  if (is.null(xmat1) || is.null(xmat2)) stop("xmat not given")
  
  # transformer: the transformer function used before
  # provided to the basis function
  # centers, scales: length=p, scaling the theta_k, k=1,...,p
  p <- dim(gam1)[2]
  if (!is.matrix(xmat1) || ncol(xmat1)!=p) stop("Invalid xmat")
  if (!is.matrix(xmat2) || ncol(xmat2)!=p) stop("Invalid xmat")
  if (nrow(xmat1)!=nrow(xmat2)) stop("Invalid xmat")
  # old_par <- par(mfrow=c(p,p))
  if (length(centers)==1) centers <- rep(centers, p)
  if (length(scales)==1) scales <- rep(scales, p)
  if (is.function(transformer)) {
    transformer <- rep(list(transformer),p)
  } else if (!is.list(transformer) || !all(sapply(transformer,is.function)) || length(transformer)!=p) {
    stop("transformer should be 1 or list of p functions.")
  }
  err_mat <- matrix(NA, nrow = p, ncol = p)
  for (j in 1:p) {
    for (k in 1:p) {
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
      xx1 <- xmat1[,k]
      xx2 <- xmat2[,k]
      if (normalize) {
        zz1 <- scale(xx1,center=centers[k],scale=scales[k])
      }
      zz1 <- transformer[[k]](zz1)
      y_est <- as.vector(fda::eval.fd(zz1,fda::fd(gamjk,basis[[k]])))
      if (!is.null(true_fun)) y_true <- true_fun[[j]][[k]](xx2)
      if (subtract_mean) {
        y_est <- y_est - mean(y_est)
        y_true <- y_true - mean(y_true)
      }
      err_mat[j,k] <- mean((y_est-y_true)^2)
    }
  }
  return(err_mat)
}



ftheta.relerr.sumj <- function(basis, gam1, gam0,
                               true_fun, a0,
                               theta_est, theta_true,
                               Bc = 0,
                               transformer = function(x) {x}, 
                               centers = 0, scales = 1,
                               subtract_mean = FALSE,
                               rm_id = c()) {
  ## calculate the error matrix along a given sequence
  # transformer: the transformer function used before
  # provided to the basis function
  # centers, scales: length=p, scaling the theta_k, k=1,...,p
  p <- dim(gam1)[2]
  nt <- nrow(theta_est)
  if (!is.matrix(theta_est) || ncol(theta_est)!=p) stop("Invalid theta_mat")
  if (!is.matrix(theta_true) || ncol(theta_true)!=p) stop("Invalid theta_mat")
  if (nrow(theta_est)!=nrow(theta_true)) stop("Invalid theta_mat")
  if (length(centers)==1) centers <- rep(centers, p)
  if (length(scales)==1) scales <- rep(scales, p)
  if (is.function(transformer)) {
    transformer <- rep(list(transformer),p)
  } else if (!is.list(transformer) || !all(sapply(transformer,is.function)) || length(transformer)!=p) {
    stop("transformer should be 1 or list of p functions.")
  }
  f_est <- matrix(0,nrow=nt,ncol=p)
  f_true <- matrix(0,nrow=nt,ncol=p)
  f_est <- sweep(f_est,2,gam0,"+") 
  f_true <- sweep(f_true,2,a0,"+")
  for (k in 1:p) {
    thetak_est <- theta_est[,k]
    thetak_true <- theta_true[,k]
    trans_thetak_est <- scale(thetak_est,center=centers[k],scale=scales[k])
    trans_thetak_est <- transformer[[k]](trans_thetak_est)
    nbasis_k <- basis[[k]]$nbasis
    nL <- nbasis_k - length(rm_id)
    if (length(Bc)==1) {Bc_k <- rep(Bc, nL)
    } else { Bc_k <- Bc[((k-1)*nL+1):(k*nL)] }
    tmp <- numeric(nbasis_k)
    if (length(rm_id)!=0) {tmp[-rm_id] <- Bc_k
    } else { tmp <- Bc_k }
    Bc_k <- tmp
    for (j in 1:p) {
      gamjk <- gam1[,k,j]
      if (length(gamjk)==basis[[k]]$nbasis-length(rm_id)) {
        tmp <- numeric(length(rm_id)+length(gamjk))
        if (length(rm_id)!=0) {tmp[-rm_id] <- gamjk
        } else { tmp <- gamjk }
        gamjk <- tmp
      } else { stop("Length error") }
      Bmat_k <- fda::eval.basis(c(trans_thetak_est), basis[[k]])
      Bmat_k <- scale(Bmat_k, center = Bc_k, scale = FALSE)
      fjk_est <- as.vector(Bmat_k %*% gamjk)
      fjk_true <- true_fun[[j]][[k]](thetak_true)
      f_est[,j] <- f_est[,j] + fjk_est
      f_true[,j] <- f_true[,j] + fjk_true
    }
  }
  return(f_est - f_true)
}






get.bcd.ode.fidelity <- function(res, type = c("all", "data", "ode")) {
  repRbind <- function(v, R) {
    return(do.call(rbind,replicate(R,v,simplify=FALSE)))
  }
  keepDim <- function(a, dims) {
    return(apply(a,dims,as.vector))
  }
  neg_log_likelihood <- function(ydata, theta_est, fam) {
    b <- fam$b
    # negative log-likelihood for exp.family
    theta_est <- as.matrix(theta_est); ydata <- as.matrix(ydata)
    nt <- dim(theta_est)[1]; p <- dim(theta_est)[2]
    if (dim(ydata)[1]!=nt || dim(ydata)[2]!=p) stop("preds, y, df not matched")
    btheta_est <- b(theta_est)
    return(-1/(nt*p) * sum(ydata * theta_est - btheta_est))  # or -1/(nt)
  }
  source("R/utils/transformer.R", local = TRUE)
  type <- match.arg(type)
  ydata <- res$inputs$ydata
  Hmat <- res$inputs$Hmat
  H4int <- res$inputs$H4int
  dH4int <- res$inputs$dH4int
  p  <- dim(ydata)[2]; nt <- dim(ydata)[1]
  R <- res$inputs$R
  rm_id <- res$inputs$rm_id
  
  theta_bspl <- Hmat %*% res$Cmat
  theta_bspl_4int <- H4int %*% res$Cmat
  dtheta_bspl_4int <- dH4int %*% res$Cmat
  Bmat <- form_Bmat(theta_bspl_4int, res$inputs$basis_f, normalized = FALSE,
                    center = res$inputs$center_th, scale = res$inputs$scale_th)
  dtheta_ode_4int <- 
    do.call(cbind,
            lapply(1:p, function(j) {
              res$gam0[j] + (Bmat %*% as.vector(res$gam1[,,j]))
            } )
    )
  out <- list()
  if (type %in% c("all", "data")) {
    fidelity.data <-
      neg_log_likelihood(ydata,
                         do.call(rbind, replicate(R, theta_bspl, FALSE)), res$inputs$fam)
    out <- append(out, list(fidelity.data = fidelity.data))
  }
  if (type %in% c("all", "ode")){
    fidelity.ode <- Metrics::mse(dtheta_ode_4int, dtheta_bspl_4int)
    out <- append(out, list(fidelity.ode = fidelity.ode))
  }
  return(out)
}





#### performance recorder ====

perf.calculator <- function(Cmat, gam0, gam1, inputs, truth, subtract_mean=FALSE,
                            transformer = pracma::sigmoid, rm_id=1) {
  source("R/utils/utils.R",local = TRUE)
  theta_mat <- inputs$H4int %*% Cmat
  
  neglik <- get.fidelity(Cmat, gam0, gam1, inputs, type = "data")$fidelity.data
  theta.mse <- mean((fda::eval.basis(truth$td, inputs$basis_th) %*% Cmat - truth$etadense)^2)
  dtheta.mse <- mean((fda::eval.basis(truth$td, inputs$basis_th, 1) %*% Cmat - truth$d.etadense)^2)
  
  graph <- apply(gam1!=0,c(2,3),any)
  
  err_mat <- f.relerr(inputs$basis_f, gam1,
                      true_fun = truth$true_fun,
                      range1 = apply(theta_mat,2,min),
                      range2 = apply(theta_mat,2,max),
                      transformer = transformer,
                      normalize = TRUE,
                      centers = inputs$center_th,
                      scales = inputs$scale_th,
                      subtract_mean = subtract_mean,
                      rm_id = rm_id)
  f.mse = mean(err_mat)
  f.mse.t = mean(err_mat[truth$graph])
  f.mse.f = mean(err_mat[!truth$graph])
  
  err_mat <- ftheta.relerr(inputs$basis_f, gam1,
                           true_fun = truth$true_fun,
                           xmat1 = theta_mat,
                           xmat2 = truth$etadense,
                           transformer = transformer,
                           normalize = TRUE,
                           centers = inputs$center_th,
                           scales = inputs$scale_th,
                           subtract_mean = subtract_mean,
                           rm_id = rm_id)
  ftheta.mse = mean(err_mat)
  ftheta.mse.t = mean(err_mat[truth$graph])
  ftheta.mse.f = mean(err_mat[!truth$graph])
  
  err_mat.sumj <-
    ftheta.relerr.sumj(inputs$basis_f, gam1, gam0,
                       truth$true_fun, truth$a0,
                       theta_est = theta_mat,
                       theta_true = truth$etadense,
                       transformer = transformer,
                       centers = inputs$center_th,
                       scales = inputs$scale_th,
                       rm_id = rm_id)
  ftheta.mse.sumj <- mean(err_mat.sumj^2)
  
  NT <- sum(truth$graph); NF <- sum(!truth$graph)
  tpr <- sum(graph[truth$graph])/NT; fpr <- sum(graph[!truth$graph])/NF
  return(list(
    neglik = neglik, theta.mse = theta.mse, dtheta.mse = dtheta.mse,
    f.mse = f.mse, f.mse.t = f.mse.t, f.mse.f = f.mse.f,
    ftheta.mse = ftheta.mse, ftheta.mse.t = ftheta.mse.t, ftheta.mse.f = ftheta.mse.f,
    ftheta.mse.sumj = ftheta.mse.sumj,
    tpr = tpr, fpr = fpr
  ))
}

record.perf <- function(perf.recorder, Cmat, gam0, gam1, inputs, truth) {
  res <- perf.calculator(Cmat, gam0, gam1, inputs, truth)
  if (!is.list(perf.recorder) || is.null(names(perf.recorder))) stop("perf.recorder: wrong format")
  for (v in names(perf.recorder)) {
    perf.recorder[[v]] <- c(perf.recorder[[v]], res[[v]])
  }
  return(perf.recorder)
}
plot.record <- function(perf.recorder, nmax = NULL) {
  if (is.null(nmax)) {nmax <- min(sapply(perf.recorder, length))}
  for (i in seq_along(perf.recorder)) perf.recorder[[i]] <- perf.recorder[[i]][1:nmax]
  dat <- data.frame(do.call(cbind, perf.recorder))
  names(dat) <- names(perf.recorder)
  dat$ii <- 1:nrow(dat)
  caps <- c("Negative log-likelihood", "Process MSE", "Deriv MSE", "f MSE")
  names(caps) <- c("neglik", "theta.mse", "dtheta.mse", "f.mse")
  dat <- dat %>% pivot_longer(cols = names(perf.recorder),
                              names_to = "Var",
                              values_to = "Val")
  dat$Var <- factor(dat$Var, levels = c("neglik", "theta.mse", "dtheta.mse", "f.mse"))
  fig <- dat %>%
    ggplot(aes(x = ii, y = Val)) +
    geom_line(alpha=0.6, color = "#0BABBD", size=1)+
    facet_wrap(vars(Var), scales = "free", labeller = labeller(Var = caps)) +
    theme_bw() + labs(x="",y="")
  print(fig)
}



perf.calculator.2stage <- function(basis_f, gam1, gam0, graph, theta, dtheta,
                                   theta.center, theta.scale,
                                   fidelity.data,
                                   truth, subtract_mean=FALSE,
                                   Bc = 0,
                                   transformer = function(x) {x}, rm_id=1) {
  neglik <- fidelity.data
  theta.mse <- mean((theta - truth$etadense)^2)
  dtheta.mse <- mean((dtheta - truth$d.etadense)^2)
  err_mat <- f.relerr(basis_f, gam1,
                      true_fun = truth$true_fun,
                      range1 = apply(theta,2,min),
                      range2 = apply(theta,2,max),
                      normalize = TRUE,
                      transformer = transformer,
                      centers = theta.center,
                      scales = theta.scale,
                      subtract_mean = TRUE,
                      rm_id = rm_id)
  f.mse <- mean(err_mat)
  f.mse.t <- mean(err_mat[truth$graph])
  f.mse.f <- mean(err_mat[!truth$graph])
  
  err_mat <- ftheta.relerr(basis_f, gam1,
                       true_fun = truth$true_fun,
                       xmat1 = theta,
                       xmat2 = truth$etadense,
                       normalize = TRUE,
                       transformer = transformer,
                       centers = theta.center,
                       scales = theta.scale,
                       subtract_mean = TRUE,
                       rm_id = rm_id)
  ftheta.mse <- mean(err_mat)
  ftheta.mse.t <- mean(err_mat[truth$graph])
  ftheta.mse.f <- mean(err_mat[!truth$graph])
  
  err_mat.sumj <-
    ftheta.relerr.sumj(basis_f, gam1, gam0,
                       truth$true_fun, truth$a0,
                       theta_est = theta,
                       theta_true = truth$etadense,
                       transformer = transformer,
                       centers = theta.center,
                       scales = theta.scale,
                       Bc = Bc,
                       rm_id = rm_id)
  ftheta.mse.sumj <- mean(err_mat.sumj^2)
  
  tpr = sum(graph[truth$graph])/sum(truth$graph)
  fpr = sum(graph[!truth$graph])/sum(!truth$graph)
  
  return(
    list(
      neglik = neglik, theta.mse = theta.mse, dtheta.mse = dtheta.mse,
      f.mse = f.mse, f.mse.t = f.mse.t, f.mse.f = f.mse.f,
      ftheta.mse = ftheta.mse, ftheta.mse.t = ftheta.mse.t, ftheta.mse.f = ftheta.mse.f,
      ftheta.mse.sumj = ftheta.mse.sumj,
      tpr = tpr, fpr = fpr
    )
  )
}

get.JadeModelScore <- function(jadefit) {
  return(
    with(
      jadefit,
      log(fidelity$fidelity.ode) +
        sum(graph)*log(inputs$p*inputs$nt4int)/(inputs$p*inputs$nt4int)
    )
  )
}
