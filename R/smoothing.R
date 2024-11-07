
smoothing <- function(family = c("gaussian", "poisson", "binomial"),
                      smooth_package = c("stats", "fda", "gss"),
                      inputs, loglam = NULL,
                      by_basis = TRUE, eval_arg = NULL){
  # smoothing B-spline
  # loglam should be a vector containing values of log(lambda)
  # loglam is only used for the case family=="gaussian" & package=="fda"
  # other package with choose the parameter automatically
  # by_basis: 
  #   TRUE: turning the smoothing results to basis expansion
  if (!by_basis && is.null(eval_arg)) stop("wrong setting")
  if (!by_basis && smooth_package=="fda") smooth_package = "stats" 
  if (by_basis) {
    if (family=="gaussian") {
      Cmat <- smooth.gaussian(as.matrix(inputs$ydata), inputs$tvals,
                              basis_obj = inputs$basis_th, eval_arg = NULL,
                              smooth_package, loglam) }
    else {
      Cmat <- smooth.exp(as.matrix(inputs$ydata), inputs$tvals,
                         family, inputs$dist_args,
                         basis_obj = inputs$basis_th, eval_arg = NULL,
                         smooth_package)
    }
    return(Cmat)
  } else {
    if (family=="gaussian") {
      smth <- smooth.gaussian(as.matrix(inputs$ydata), inputs$tvals,
                              basis_obj = NULL, eval_arg = eval_arg,
                              smooth_package, loglam) }
    else {
      smth <- smooth.exp(as.matrix(inputs$ydata), inputs$tvals,
                         family, inputs$dist_args,
                         basis_obj = NULL, eval_arg = eval_arg,
                         smooth_package)
    }
  }
}


smooth.gaussian <- function(Y, x, basis_obj = NULL, eval_arg = NULL,
                            smooth_package = c("stats", "fda"),
                            loglam){
  p <- ncol(Y)
  by_basis = !is.null(basis_obj)
  if (smooth_package=="fda") {
    gcvsave <- matrix(NA, length(loglam), p)
    for (i in 1:length(loglam)){
      lambda_i <- 10^loglam[i]
      fdPar_i <- fda::fdPar(basis_obj, Lfdobj=NULL, lambda_i)
      gcvsave[i,] <- fda::smooth.basis(x, Y, fdPar_i)$gcv
    }
    mingcv.idx = apply(gcvsave, 2, which.min)
    Cmat <- do.call(
      cbind, lapply(1:p,
                    function(j) {
                      fdPar_j <- fda::fdPar(basis_obj, Lfdobj=NULL, 10^loglam[mingcv.idx[j]])
                      return(fda::smooth.basis(x, Y[,j], fdPar_j)$fd$coefs)
                    })
    )
    return(Cmat)
  } else if (smooth_package=="stats") {
    if (by_basis) {
      theta_smooth <- do.call(cbind, lapply(1:p, function(j) {
        smth <- smooth.spline(x, Y[,j], all.knots = T)
        return(predict(smth, x)$y)
      }) )
      smth_obj <- fda::smooth.basis(x, theta_smooth, fdParobj = basis_obj)
      Cmat <- smth_obj$fd$coefs
      return(Cmat)
    } else {
      smth_obj = lapply(1:p, function(j) {
        smth = smooth.spline(x, Y[,j], all.knots = T)
        return(smth)
      })
      theta_smooth = sapply(1:p, function(j) {
        predict(smth_obj[[j]], eval_arg)$y})
      dtheta_smooth = sapply(1:p, function(j) {
        predict(smth_obj[[j]], eval_arg, deriv=1)$y})
      return(list(theta = theta_smooth, dtheta = dtheta_smooth))
    }
    
  }
  
}


smooth.exp <- function(Y, x, basis_obj = NULL, eval_arg = NULL,
                       family = c("poisson", "binomial"), dist_args = NULL,
                       smooth_package = c("gss")){
  # stop("not implemented")
  p <- ncol(Y); by_basis = !is.null(basis_obj)

  get.y4gss <- function(y, family) {
    if (family %in% c("poisson")) {
      gssfamily <- family
      invfun <- function(x) exp(x)
      response <- y
    } else if (family=="binomial") {
      gssfamily <- "binomial"
      invfun <- function(x) 1-1/(1+exp(x))
      response <- y  # n should be defined in ...
    } else if (family=="nbinomial") {
      n <- dist_args$n
      gssfamily <- "nbinomial"
      invfun <- function(x) n*(1-1/(1+exp(x)))
      response <- cbind(y, n)  # n should be defined in ...
    } else if (family=="negbinomial") {
      nu <- dist_args$nu
      gssfamily <- "nbinomial"
      invfun <- function(x) {nu/exp(x)}  # nu should be defined in ...
      response <- cbind(y, nu)
    }
    return(list(
      response = response,
      gssfamily = gssfamily,
      invfun = invfun
    ))
  }
  
  if (by_basis) {
    Cmat <- do.call(
      cbind,
      lapply(1:ncol(Y), function(j) {
        gss.param <- get.y4gss(Y[,j], family)
        response <- gss.param$response; gssfamily <- gss.param$gssfamily
        invfun <- gss.param$invfun
        fit <- gss::gssanova(response ~ x, family = gssfamily)
        est <- predict(fit, data.frame(x=x), se=FALSE)
        cj <- fda::smooth.basis(x, est, basis_obj)$fd$coefs
        return(cj)
      })
    )
    return(Cmat)
  } else {
    smth_obj <- lapply(1:ncol(Y), function(j) {
        gss.param <- get.y4gss(Y[,j], family)
        response <- gss.param$response; gssfamily <- gss.param$gssfamily
        invfun <- gss.param$invfun
        fit <- gss::gssanova(response ~ x, family = gssfamily)
        return(fit)
      })
    fn <- function(x,j) {predict(smth_obj[[j]], data.frame(x=x), se=FALSE)}
    theta_smooth = dtheta_smooth = matrix(NA,length(eval_arg),p)
    for (j in 1:p) {
      theta_smooth[,j] = fn(eval_arg, j)
      dtheta_smooth[,j] = pracma::fderiv(fn, x=eval_arg, j=j)
    }
    return(list(theta = theta_smooth, dtheta = dtheta_smooth))
  }
  
}

