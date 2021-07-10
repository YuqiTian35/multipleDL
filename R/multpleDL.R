#' CPMs for multiple DLs
#'
#' This function build the CPM for multiple DLs
#'
#' @param formula a R formula object
#' @param data a data frame
#' @param delta_lower (optional) indicators of lower DLs censoring (1: observed; 0:censored). If not specified, treat as observed.
#' @param delta_upper (optional) indicators of upper DLs censoring(1: observed; 0:censored). If not specified, treat as observed.
#' @param link the link function
#' @return A list of results: coefficients, covariance matrix, unique response values, number of alphas, number of betas, link functions
#' @export
multipleDL <- function(formula, data, delta_lower = NULL, delta_upper = NULL, link){

  mf <- model.frame(formula=formula, data=data)
  terms <- attr(mf, "terms")
  x <- as.matrix(model.matrix(terms, data=mf)[,-1]) # exclude intercept
  z <- model.response(mf)

  ###########
  # DL indicators
  if(is.null(delta_lower)){
    # if no delta_lower => all observed
    delta_lower <- rep(1, length(z))
  }else{
    if(length(delta_lower) != nrow(data)){
      stop("delta_lower: length of delta_lower should equal the number of rows of data")
    }
  }

  if(is.null(delta_upper)){
    # if no delta_upper => all observed
    delta_upper <- rep(1, length(z))
  }else{
    if(length(delta_upper) != nrow(data)){
      stop("delta_upper: length of delta_upper should equal the number of rows of data")
    }
  }

  # link functions
  links <- c('logit', 'probit', 'loglog', 'cloglog')
  ilinks <- as.integer(match(link, links, -1))
  if (ilinks < 1){
    stop("link function: set to be logit, probit, loglog or cloglog")
  }

  # prepare data
  h <- 0.1 # to indicate the smallest/largest level
  N <- nrow(x) # number of subjects
  p <- ncol(x) # number of parameters
  inf <- 1e6 # represent infinity

  # QR decomposition
  qrX <- selectedQr(x)
  X <- qrX$X
  xbar <- qrX$xbar
  Rinv <- qrX$Rinv

  # lower DLs (dl^l_1 < dl^l_2 < ...)
  dl_l <- sort(unique(z[delta_lower == 0]))
  # uppder DLs (dl^u_1 > dl^l_2 > ...)
  dl_u <- sort(unique(z[delta_upper == 0]), decreasing = TRUE)
  # set S: all observed values
  S <- sort(unique(z[delta_lower & delta_upper]))

  # check if there is no observations between two DLs
  if(length(dl_l) > 1){
    ind_l <- which(table(cut(S, c(-inf, dl_l, inf)))[2:length(dl_l)] == 0)
    sprintf("No observed values between lower detection limits %f and %f. Merging them into one detection limit %f",
          dl_l[ind_l], dl_l[ind_l+1L], dl_l[ind_l+1L])
    # keep the larger lower DL
    if(length(ind_l) > 0) dl_l <- dl_l[-ind_l]

  }

  # check if there is no observations between two DLs
  if(length(dl_u) > 1){
    ind_u <- which(table(cut(S, c(inf, dl_u, -inf)))[2:length(dl_u)] == 0)
    sprintf("No observed values between upper detection limits %f and %f. Merging them into one detection limit %f",
            dl_u[ind_u-1L], dl_u[ind_u], dl_u[ind_u])
    # keep the smaller upper DL
    if(length(ind_u) > 0) dl_u <- dl_u[-(ind_u-1)]
  }


  # check if z_(0) is needed
  if(length(dl_l) > 0 & dl_l[1] <= S[1]){
    S <- c(S[1] - h, S)
    z_0 <- TRUE
  }else{
    z_0 <- FALSE
  }

  # check if z_(J+1) is needed
  if(length(dl_u) > 0 & dl_u[1] >= S[length(S)]){
    S <- c(S, S[length(S)] + h)
    z_J1 <- TRUE
  }else{
    z_J1 <- FALSE
  }

  # number of unique values
  J <- length(S)


  ## indicator
  delta <- ifelse(delta_lower & delta_upper, 1,
                  ifelse(!delta_lower & z_0 & z == dl_l[1], 12,
                         ifelse(!delta_upper & z_J1 & z == dl_u[1], 13,
                                ifelse(!delta_lower, 2, 3))))

  ## code value
  code_value <- sapply(1:N, function(i){
    ## [at least one of delta_l and delta_u == 1]
    if(delta_lower[i] & delta_upper[i]){
      return(z[i])
    }else if(!delta_lower[i] & z_0 & z[i] == dl_l[1]){
      return(S[1])
    }else if(!delta_upper[i] & z_J1 & z[i] == dl_u[1]){
      return(S[length(S)])
    }else if(!delta_lower[i]){
      return(S[max(which(S < z[i]))])
    }else if(!delta_upper[i]){
      S[min(which(S > z[i]))]
    }else{
      return(NA)
    }
  })



  # rank the code_value
  j <- match(code_value, sort(unique(code_value)))

  # data for stan code
  link_num <- func_link_num(link)
  sds <- 1e2
  data_stan <- list(N = N, p = p, X = X, J = J, j = j, delta = delta,
                    link_num = link_num, sds = sds)

  # stan code
  mod <- stanmodels[['multipe_dls_cpm']]
  stancode <- rstan::get_stancode(mod)
  res.stan <- rstan::optimizing(mod,
                                init = '0',
                                iter = 1e5,
                                tol_obj = 1e-5,
                                data = data_stan)

  # coefficiencts
  beta <- c(matrix(res.stan$par[grep("beta", names(res.stan$par))], nrow=1) %*% t(Rinv))
  names(beta) <- attr(terms, 'term.labels')
  alpha <- res.stan$par[grep("alpha", names(res.stan$par))] + sum(beta * xbar)
  # alpha <- res.stan$par[grep("alpha", names(res.stan$par))] - sum(beta * xbar) # orm version
  coef <- c(alpha, beta)

  # covariance
  fam <- func_link(link)
  var <- func_V(coef = coef, n = N, x = x, y = j, delta = delta,
                k = J, p = p, fam = fam)
  rownames(var) <- colnames(var) <-  names(coef)

  return(list(coef = coef, var = var,
               yunique = sort(unique(code_value)),
               kint = J-1, p = p, fam = fam))

}

