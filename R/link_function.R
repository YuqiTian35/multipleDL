#' Link functions (number)
#'
#' This function faciliates the stan code (used as an internal function)
#'
#' @param link the link function
#' @return An integer representing corresponding link function
func_link_num <- function(link){
  if(link == 'logit'){
    return(1)
  }else if(link == 'probit'){
    return(2)
  }else if(link == 'loglog'){
    return(3)
  }else if(link == 'cloglog'){
    return(4)
  }
}

#' Link functions
#'
#' This function includes necessary functions related to each link function
#'
#' @param link the link function
#' @return A list of functions subject to a link function
#' @export
func_link <- function(link){
  families <-
    list(logit =
           list(cumprob=function(x)    1 / (1 + exp(-x)),
                inverse=function(x)    log(x / (1 - x)),
                deriv  =function(x, f) f * (1 - f),
                deriv2 =function(x, f, deriv) f * (1 - 3*f + 2*f*f),
                name = 'logit'),

         probit =
           list(cumprob=pnorm,
                inverse=qnorm,
                deriv  =function(x, ...)      dnorm(x),
                deriv2 =function(x, f, deriv) - deriv * x,
                name = 'probit'),
         loglog =
           list(cumprob=function(x)      exp(-exp(-x)),
                inverse=function(x)      -log(-log(x    )),
                deriv  =function(x, ...) exp(-x - exp(-x)),
                deriv2 =function(x, ...) ifelse(abs(x) > 200, 0,
                                                exp(-x - exp(-x)) * (-1 + exp(-x))),
                name = 'loglog'),
         cloglog =
           list(cumprob=function(x)      1 - exp(-exp(x)),
                inverse=function(x)      log(-log(1 - x)),
                deriv  =function(x, ...) exp( x - exp( x)),
                deriv2 =function(x, f, deriv) ifelse(abs(x) > 200, 0,
                                                     deriv * ( 1 - exp( x))),
                name = 'cloglog'))

  return(families[[link]])
}
