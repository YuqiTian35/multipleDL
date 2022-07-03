#' Calculate the covariance matrix
#'
#' This functions calculates the covariance matrix based on the point estimates
#'
#' @param coef coefficients (alpha, beta)
#' @param n number of subjects
#' @param x original covariate matrix
#' @param y ranks of code values
#' @param delta censoring indicators
#' @param k the number of unique code values
#' @param p the number of covariates
#' @param fam a list of functions subject to the link function
#' @return A covariance matrix of coefficients
func_V <- function(coef, n, x, y, delta, k, p, fam){
  kint <- k - 1L # k-1
  allp <- as.integer(kint + p)

  f   <- fam$cumprob
  fp  <- fam$deriv
  fpp <- fam$deriv2

  xb <- c(x %*% coef[-(1L : kint)])
  ints <- c(-1e100, coef[1:kint], 1e100) # a0, a1,.., an
  xby <- ints[y] - xb #a_{i-1}+beta*x_i (starting from i=1)
  xby1 <- ints[y + 1L] - xb #a_{i}-beta*x_i
  fa <- f(xby)
  fb <- f(xby1)
  fpa  <- fp(xby, fa)
  fpb  <- fp(xby1, fb)
  fppa <- fpp(xby, fa, fpa)
  fppb <- fpp(xby1, fb, fpb)
  P  <- (fb - fa)^(delta %in% c(1, 12, 13)) * fb^(delta == 2) * (1 - fa)^(delta == 3)

  l <- if(kint == 1){
    allp ^ 2
  }else{
    as.integer(p * p + 2L * kint * p+ 3L * kint - 2L)
  }
  lia <- as.integer(allp + 1L)

  w <- hessian(n, allp, kint, p, x, y, delta, P, fa, fb, fpa, fpb, fppa, fppb, l, lia)
  V <- if(kint == 1L){
    matrix(w$v, nrow=allp, ncol=allp)
  }else new('matrix.csr', ra=w$v, ja=as.integer(w$ja), ia=as.integer(w$ia), dimension=c(allp, allp))
  V <- (V + SparseM::t(V))/2   # force symmetry; chol() complains if 1e-15 off
  inv_V <- SparseM::as.matrix(SparseM::solve(V, tol=1e-7))

  return(inv_V)
}
