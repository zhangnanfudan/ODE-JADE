get.family <- function(family = c("gaussian", "poisson",
                                  "binomial", "nbinomial",
                                  "gamma", "negbinomial"),
                       ...) {
  family <- match.arg(family)
  if (family=="gaussian") {
    b <- function(x) return(0.5*x^2)
    db <- function(x) return(x)
    d2b <- function(x) return(structure(rep(1, length(x)), dim = dim(x)))
    dbinv <- function(theta) return(theta)
    linkfun <- function(mu) return(mu)
    dlinkfun <- function(mu) return(structure(rep(1, length(mu)), dim = dim(mu)))
    linkinv <- function(x) return(x)
    thetafun <- function(mu) return(mu)
    dthetafun <- function(mu) return(structure(rep(1, length(mu)), dim = dim(mu)))
  } else if (family=="poisson") {
    b <- function(x) return(exp(x))
    db <- function(x) return(exp(x))
    d2b <- function(x) return(exp(x))
    dbinv <- function(theta) return(log(theta))
    linkfun <- function(mu) return(log(mu))
    dlinkfun <- function(mu) return(1/mu)
    linkinv <- function(x) return(exp(x))
    thetafun <- function(mu) return(log(mu))
    dthetafun <- function(mu) return(1/mu)
  } else if (family=="binomial") {
    b <- function(x) return(log(1+exp(x)))
    db <- function(x) return(exp(x)/(1+exp(x)))
    d2b <- function(x) {
      tmp <- exp(x)/(1+exp(x))
      return(tmp * (1-tmp))
    }
    dbinv <- function(theta) return(log(theta) - log(1-theta))
  } else if (family=="nbinomial") {
    stop('NOT FINISHED')
    # n should be in global environment
    b <- function(x) return(n*log(1+exp(x)))
    db <- function(x) return(n*exp(x)/(1+exp(x)))
    d2b <- function(x) {
      tmp <- exp(x)/(1+exp(x))
      return(n * tmp * (1-tmp))
    }
    dbinv <- function(theta) return(1/n*(log(theta)-log(1-theta)))
    linkfun <- function(mu) return(mu/(n-mu))
    dlinkfun <- function(mu) return(n/(n-mu)^2)
    linkinv <- function(x) return(n*x/(x+1))
    thetafun <- function(mu) return(log(linkfun(mu)))
    dthetafun <- function(mu) return(dlinkfun(mu)/linkfun(mu))
  } else if (family=="negbinomial") {
    stop('NOT FINISHED')
    # nu should be in global environment
    b <- function(x) return( -nu*log(1-exp(x)) )
    db <- function(x) return( nu/(exp(-x)-1) )
    d2b <- function(x) return( nu/(exp(-x)+exp(x)-2) )
    linkinv <- function(eta) return(nu/exp(eta))
    # stop("unimplemented")
  } else if (family=="gamma") {
    stop('NOT FINISHED')
    b <- function(x) return(-log(-x))
    db <- function(x) return(-1/x)
    d2b <- function(x) return(1/x^2)
    stop("unimplemented")
  } 
  
  return(list(
    familyname=family,
    b = b,
    db = db,
    d2b = d2b,
    dbinv = dbinv
  ))
}

