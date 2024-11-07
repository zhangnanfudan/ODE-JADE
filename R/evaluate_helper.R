get.fjkval.jade <- function(basisk, gamjk,
                            transformer = pracma::sigmoid,
                            range1=0, range2=1, by=0.01,
                            centerk = 0, scalek = 1,
                            subtract_mean = TRUE,
                            rm_id = 1) {
  tmp <- numeric(length(rm_id)+length(gamjk))
  if (length(rm_id)!=0) {tmp[-rm_id] <- gamjk
  } else { tmp <- gamjk }
  gamjk <- tmp
  theta <- seq(range1,range2,by=by)
  sigma_theta <- transformer(scale(theta,centerk,scalek))
  y <- as.vector(fda::eval.fd(sigma_theta,fda::fd(gamjk,basisk)))
  return(list(
    x = theta,
    y = if (subtract_mean) {y - mean(y)} else {y}
  ))
}

get.fjkval.2stage <- function(basisk, gamjk,
                              range1=0, range2=1, by=0.01,
                              centerk = 0, scalek = 1,
                              subtract_mean = TRUE,
                              rm_id = 1) {
  tmp <- numeric(length(rm_id)+length(gamjk))
  if (length(rm_id)!=0) {tmp[-rm_id] <- gamjk
  } else { tmp <- gamjk }
  gamjk <- tmp
  theta <- seq(range1,range2,by=by)
  thetabar <- scale(theta,centerk,scalek)
  y <- as.vector(fda::eval.fd(thetabar,fda::fd(gamjk,basisk)))
  return(list(
    x = theta,
    y = if (subtract_mean) {y - mean(y)} else {y}
  ))
}

get.fjkval.truth <- function(true_fjk,
                             range1=0, range2=1, by=0.01,
                             subtract_mean = TRUE) {
  theta <- seq(range1,range2,by=by)
  y <- true_fjk(theta)
  return(list(
    x = theta,
    y = if (subtract_mean) {y - mean(y)} else {y}
  ))
}


get.fjk_theta.jade <- function(basisk, gamjk, theta,
                               transformer = pracma::sigmoid,
                               centerk = 0, scalek = 1,
                               subtract_mean = TRUE,
                               rm_id = 1) {
  tmp <- numeric(length(rm_id)+length(gamjk))
  if (length(rm_id)!=0) {tmp[-rm_id] <- gamjk
  } else { tmp <- gamjk }
  gamjk <- tmp
  sigma_theta <- transformer(scale(theta,centerk,scalek))
  y <- as.vector(fda::eval.fd(sigma_theta,fda::fd(gamjk,basisk)))
  return(list(
    x = theta,
    y = if (subtract_mean) {y - mean(y)} else {y}
  ))
}

get.fjk_theta.2stage <- function(basisk, gamjk, theta,
                                 centerk = 0, scalek = 1,
                                 subtract_mean = TRUE,
                                 Bc_k = 0,
                                 rm_id = 1) {
  tmp <- numeric(length(rm_id)+length(gamjk))
  if (length(rm_id)!=0) {tmp[-rm_id] <- gamjk
  } else { tmp <- gamjk }
  if (length(Bc_k)==1) Bc_k <- rep(Bc_k, length(gamjk))
  gamjk <- tmp
  thetabar <- scale(theta,centerk,scalek)
  Bk <- fda::eval.basis(c(thetabar), basisk)
  if (length(rm_id)!=0) {
    Bk[,-rm_id] <- scale(Bk[,-rm_id], center = Bc_k, scale = FALSE)
  } else { Bk <- scale(Bk, center= Bc_k, scale = FALSE) }
  y <- as.vector(Bk %*% gamjk)
  return(list(
    x = theta,
    y = if (subtract_mean) {y - mean(y)} else {y}
  ))
}

get.fjk_theta.truth <- function(true_fjk, theta,
                                subtract_mean = TRUE) {
  y <- true_fjk(theta)
  return(list(
    x = theta,
    y = if (subtract_mean) {y - mean(y)} else {y}
  ))
}

get.fjk_theta.mse.2stage <- function(basis, gamj_mat, theta_mat,
                                     true_fj,
                                     center = 0, scale = 1,
                                     subtract_mean = TRUE,
                                     rm_id = 1) {
  # this funtion intends to calculate the fjk.mse
  # for a sequence of gamj
  p <- length(basis)
  nlambda <- dim(gamj_mat)[3]
  if (p != length(true_fj) || p != dim(gamj_mat)[2]) stop("dim error")
  MSE_list <- sapply(
    1:nlambda,
    function(l) {
      MSE <- 0
      for (k in 1:p) {
        residue <- 
          get.fjk_theta.2stage(basis[[k]], gamj_mat[,k,j], theta_mat[,k],
                               centerk = center[k],
                               scalek = scale[k],
                               subtract_mean = subtract_mean,
                               rm_id = rm_id) - 
          get.fjk_theta.truth(true_fj[[k]], theta_mat[,k],
                              subtract_mean = subtract_mean)
        MSE <- MSE + mean(residue^2)
      }
      return(MSE/p)
    }
  )
  return(MSE_list)
}


get.Intfjk_theta.2stage <- 
  function(basisk, gamjk, theta,
           centerk = 0, scalek = 1,
           subtract_mean = TRUE,
           Bc_k = 0,
           rm_id = 1) {
  tmp <- numeric(length(rm_id)+length(gamjk))
  if (length(rm_id)!=0) {tmp[-rm_id] <- gamjk
  } else { tmp <- gamjk }
  if (length(Bc_k)==1) rep(Bc, length(gamjk))
  gamjk <- tmp
  thetabar <- scale(theta,centerk,scalek)
  Bk <- fda::eval.basis(c(thetabar), basisk)
  if (length(rm_id)!=0) {
    Bk[,-rm_id] <- scale(Bk[,-rm_id], center = Bc_k, scale = FALSE)
  } else { Bk <- scale(Bk, center= Bc_k, scale = FALSE) }
  y <- as.vector(Bk %*% gamjk)
  return(list(
    x = theta,
    y = if (subtract_mean) {y - mean(y)} else {y}
  ))
}


which.best_graph <- function(graph.all, graph_true) {
  p <- dim(graph.all)[1]
  ng <- sum(apply(!is.na(graph.all), 3, all))
  out <- list(
    best.lam_id = numeric(p),
    tp = 0, fp = 0
  )
  for (j in 1:p) {
    if (sum(graph_true[, j]) != 0) {
      idx_best_tp <- apply(
        graph.all,
        MARGIN = 3,
        \(graph) all(graph[graph_true[, j], j])
      ) |> which()
      if (length(idx_best_tp) != 0) {
        idx_best_tp <- min(idx_best_tp)
        out$tp <- out$tp + sum(graph_true[, j])
      } else {
        tp <- apply(
          graph.all,
          MARGIN = 3,
          \(graph) sum(graph[graph_true[, j], j])
        )
        idx_best_tp <- min(which(tp == max(tp)))
        out$tp <- out$tp + max(tp)
      }
    } else {
      idx_best_tp <- 1
    }
    fp <- apply(
      graph.all,
      MARGIN = 3,
      \(graph) sum(graph[!graph_true[, j], j])
    )
    idx_best_graph_j <- which.min(fp[idx_best_tp:ng]) + idx_best_tp - 1
    out$fp <- out$fp + fp[idx_best_graph_j]
    out$best.lam_id[j] <- idx_best_graph_j
  }
  return(out)
}

