#' QR Decomposition Preserving Selected Columns
#'
##' Runs a matrix through the QR decomposition and returns the transformed matrix and the forward and inverse transforming matrices
##' `R, Rinv`.  If columns of the input matrix `X` are centered the QR transformed matrix will be orthogonal.
##'  This is helpful in understanding the transformation and in scaling prior distributions on the transformed scale.
##' `not` can be specified to keep selected columns as-is.
##' `cornerQr` leaves the last column of `X` alone (possibly after centering).
##'  When `not` is specified, the square transforming matrices have appropriate identity submatrices inserted
##' so that recreation of original `X` is automatic.
##' @param X a numeric matrix
##' @param not an integer vector specifying which columns of `X` are to be kept with their original values
##' @param corner set to `FALSE` to not treat the last column specially.  You may not specify both `not` and `corner`.
##' @param center set to `FALSE` to not center columns of `X` first
##' @return list with elements `X, R, Rinv, xbar` where `xbar` is the vector of means (vector of zeros if `center=FALSE`)
#'  @export
selectedQr <- function(X, not=NULL, corner=FALSE, center=TRUE) {
  if(center) {
    X    <- scale(X, center=TRUE, scale=FALSE)
    xbar <- as.vector(attr(X, 'scaled:center'))
  } else xbar <- rep(0., ncol(X))

  if(length(not)) {
    if(corner) stop('may not specify both not and corner=TRUE')
    p <- ncol(X)
    # Handle the case where no variables are orthogonalized
    if(length(not) == p)
      return(list(X=X, R=diag(p), Rinv=diag(p), xbar=xbar))
    Xo <- X
    X  <- Xo[, -not, drop=FALSE]
  }
  p     <- ncol(X)
  QR    <- qr(X)
  Q     <- qr.Q(QR)
  R     <- qr.R(QR)
  sgns  <- sign(diag(R))
  X     <- sweep(Q, MARGIN = 2, STATS = sgns, FUN = `*`)
  R_ast <- sweep(R, MARGIN = 1, STATS = sgns, FUN = `*`)
  cornr <- if(corner) R_ast[p, p] else 1.
  R_ast_inverse <- backsolve(R_ast, diag(p))
  X <- X * cornr
  R_ast <- R_ast / cornr
  R_ast_inverse <- R_ast_inverse * cornr

  if(length(not)) {
    Xo[, -not]       <- X
    R <- Rinv        <- diag(p + length(not))
    R[-not,    -not] <- R_ast
    Rinv[-not, -not] <- R_ast_inverse
    return(list(X=Xo, R=R, Rinv=Rinv, xbar=xbar))
  }
  list(X = X, R = R_ast, Rinv = R_ast_inverse, xbar=xbar)
}
