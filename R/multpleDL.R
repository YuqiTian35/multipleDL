#' CPMs for multiple detection limits
#'
#' This function build the CPM for multiple detection limits (DLs).
#'
#' @param formula an R formula object
#' @param data a data frame including response data and covariates
#' @param delta_lower (optional) indicators of lower DLs censoring (1: observed; 0:censored). If not specified, treat as observed.
#' @param delta_upper (optional) indicators of upper DLs censoring(1: observed; 0:censored). If not specified, treat as observed.
#' @param link the link function (probit, logit, loglog, cloglog)
#' @return A list containing the following components:
#' @return \item{coef}{a numeric vector of estimated coeffiencts}
#' @return \item{var}{covariance matrix of estimated coeffiencts}
#' @return \item{yunique}{a numeric vector of unique response values}
#' @return \item{kint}{number of alphas (intercept terms)}
#' @return \item{p}{number of betas (regression coeffiencts)}
#' @return \item{fam}{a list of functions associated with the specified link function}
#' @export
#'
#' @details When there are multiple DLs, we appropriately modify the CPM likelihood.
#' If a value is below a lower DL, set the censored value as the lower DL and set the
#' lower DL indicator `delta_lower` to be 0. Similarly, if a value is above an upper DL,
#' set the censored value as the upper DL and set the upper DL indicator `delta_upper` to be 0.
#' This function also works when there is only a single lower and/or upper DL.
#'
#' Conditional quantiles and CDFs and corresponding 95% confidence intervals can be calculated
#' from the model fit.
#'
#' @seealso \code{\link{cdf_dl}, \link{quantile_dl}}
#'
#' @references
#' Tian et al. "Addressing detection limits by semiparametric cumulative probability models." (2022) (to be submitted)
#' @references
#' Stan Development Team (2020). RSroxygen2::roxygenize()tan: the R interface to Stan. R package version 2.19.3. https://mc-stan.org
#' @references
#' Harrell, F. (2020). rms: Regression modeling strategies. R package version 6.1.0. https://CRAN.R-project.org/package=rms
#'
#'
#' @examples
#' ## Multiple DLs
#' ## generate a small example data: 3 sites with different lower and upper DLs
#' ## lower DLs: site 1: - 0.2; site 2: 0.3; site 3: no lower DL
#' ## upper DLs: site 1: no upper DL; site 2: 4; site 3: 3.5
#' ## each site includes 100 subjects
#' n <- 100
#' x <- rnorm(n * 3)
#' e <- rnorm(n * 3)
#' y <- exp(x + e)
#' no_dl <- 1e6
#' data <- data.frame(y = y, x = x, subset = rep(c(1, 2, 3), each=n))
#' data$dl_l <- ifelse(data$subset == 1, 0.2, ifelse(data$subset == 2, 0.3, -no_dl))
#' data$dl_u <- ifelse(data$subset == 1, no_dl, ifelse(data$subset == 2, 4, 3.5))
#' data$delta_l <- ifelse(data$y >= data$dl_l, 1, 0)
#' data$delta_u <- ifelse(data$y <= data$dl_u, 1, 0)
#' data$z <- ifelse(data$delta_l == 0, data$dl_l, ifelse(data$delta_u == 0, data$dl_u, data$y))
#' # model
#' mod <- multipleDL(formula = z ~ x, data = data,
#'                  delta_lower = data$delta_l, delta_upper = data$delta_u, link='probit')
#' # new data
#' new.data <- data.frame(x = c(0, 1))
#' conditional_median <- quantile_dl(mod, new.data, probs = 0.5)
#' conditional_cdf <- cdf_dl(mod, new.data, at.y = 1.5) # P(y <= 1.5 | new.data)
#'
#'
#' ## Single DL: lower DL at 0.5
#' n <- 100
#' x <- rnorm(n)
#' e <- rnorm(n)
#' y <- exp(x + e)
#' lower_dl <- 0.5
#' data <- data.frame(y = y, x = x)
#' data$delta_lower <- ifelse(data$y >= lower_dl, 1, 0)
#' data$z <- ifelse(data$delta_lower == 0, lower_dl, data$y)
#' mod <- multipleDL(formula = z ~ x, data = data,
#'                   delta_lower = data$delta_l, link='probit')
multipleDL <- function(formula, data, delta_lower = NULL, delta_upper = NULL, link){

  mf <- model.frame(formula=formula, data=data)
  terms <- attr(mf, "terms")
  x <- as.matrix(model.matrix(terms, data=mf)[,-1]) # exclude intercept
  coef_name <- colnames(model.matrix(terms, data=mf))[-1]
  z <- model.response(mf)
  yunique <- sort(unique(z))

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
    if(length(ind_l) > 0){
      dl_l <- dl_l[-ind_l]
      z[z == dl_l[ind_l]] <- dl_l[ind_l+1L]
    }
  }

  # check if there is no observations between two DLs
  if(length(dl_u) > 1){
    ind_u <- which(table(cut(S, c(inf, dl_u, -inf)))[2:length(dl_u)] == 0)
    sprintf("No observed values between upper detection limits %f and %f. Merging them into one detection limit %f",
            dl_u[ind_u-1L], dl_u[ind_u], dl_u[ind_u])
    # keep the smaller upper DL
    if(length(ind_u) > 0){
      dl_u <- dl_u[-(ind_u-1)]
      z[z == dl_u[ind_u-1L]] <- dl_u[ind_u]
    }
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
      if(length(which(S < z[i])) == 0){
        return(S[1])
      }else{
        return(S[max(which(S < z[i]))])
      }
    }else if(!delta_upper[i]){
      if(length(which(S > z[i])) == 0){
        return(S[length(S)])
      }else{
        S[min(which(S > z[i]))]
      }
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

  # coefficients
  beta <- c(matrix(res.stan$par[grep("beta", names(res.stan$par))], nrow=1) %*% t(Rinv))
  names(beta) <- coef_name # attr(terms, 'term.labels')
  alpha <- res.stan$par[grep("alpha", names(res.stan$par))] + sum(beta * xbar)
  # alpha <- res.stan$par[grep("alpha", names(res.stan$par))] - sum(beta * xbar) # orm version
  coef <- c(alpha, beta)

  # covariance
  fam <- func_link(link)
  var <- func_V(coef = coef, n = N, x = x, y = j, delta = delta,
                k = J, p = p, fam = fam)
  rownames(var) <- colnames(var) <-  names(coef)

  # log-likelihood
  px <- fam$cumprob(outer(alpha, as.vector(x %*% beta), "-"))
  ll <- sum(sapply(1:N, function(i){
    if(delta[i] == 1){
      return(log(c(px[,i], 1)[j[i]] - c(0, px[,i])[j[i]]))
    }else if(delta[i] == 12){
      return(0)
    }else if(delta[i] == 13){
      return(0)
    }else if(delta[i] == 2){
      return(log(c(px[,i], 1)[j[i]]))
    }else if(delta[i] == 3){
      return(log(1 - c(0, px[,i])[j[i]]))
    }else{
      return(NA)
    }
  }))

  return(list(coef = coef, var = var,
               yunique = yunique,
               kint = J-1, p = p, fam = fam,
               x = x, log_likelihood = ll))

}

