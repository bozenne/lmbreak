### mixedbreak.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jul 14 2026 (10:40) 
## Version: 
## Last-Updated: jul 15 2026 (10:24) 
##           By: Brice Ozenne
##     Update #: 117
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:


## * logLik_mixed10 (documentation)
##' @title Mixed Model Likelihood for 10 pattern
##' @description Evaluate the log-likelihood of a breakpoint model with 1 random slope, 1 random breakpoint, and one plateau (10 pattern) using numerical integration.
##'
##' @param Y [numeric vector] outcome values
##' @param t [numeric vector] timepoints at which the outcome value were measured.
##' @param beta [numeric] value of the population slope parameter.
##' @param psi [numeric] value of the population breakpoint.
##' @param tau_uu [numeric,>0] variance of the subject specific deviations to the population slope.
##' @param tau_vv [numeric] variance of the subject specific deviations to the population breakpoint.
##' @param tau_uv [numeric,>0] covariance between the subject specific deviations slope and breapoint.
##' @param sigma2 [numeric,>0] residual variance.
##' @param lower [numeric] lower bound for the numerical integration: below this value the density of the random effect distribution is ignored.
##' @param upper [numeric] upper bound for the numerical integration: above this value the density of the random effect distribution is ignored.
##' @param indiv [logical] should the log-likelihood be output for each observation or only its sum across all observations
##' @param integration [character] method use to compute the likelihood: \code{"bivariate"} or \code{"univariate"}.
##' The latter is much faster.
##'
##' @return A numeric value for the value of the total log-likelihood (\code{indiv=FALSE}) or a \code{data.frame} with as many rows as observations (\code{indiv=TRUE}).

## * logLik_mixed10 (examples)
##' @examples
##' #### EXAMPLE 1 ####
##' logLik_mixed10(Y = 0, t = 0, beta = 0, psi = 0, tau_uu = 1, tau_uv = 0, tau_vv = 1, sigma2 = 1, lower = -5, upper = 5, integration = "univariate")
##' ##   integral1      error1 integral2       error2  integral error    logLik
##' ## 1  0.199471 1.72908e-10 0.1575104 1.015364e-06 0.3569814    NA -1.030072
##'
##' ## first integral
##' pnorm(5)*(1-pnorm(0)) / sqrt(2*pi)
##' ## second integral
##' pnorm(0) * integrate(f = function(u){return(exp(-u^2/2) / sqrt(1+u^2))}, lower = -5, upper = 5)$value / (2*pi)
##'
##' logLik_mixed10(Y = 0, t = 0, beta = 0, psi = 0, tau_uu = 1, tau_uv = 0, tau_vv = 1, sigma2 = 1, lower = -5, upper = 5, integration = "bivariate")
##' ## integral1 error1 integral2 error2  integral        error    logLik
##' ##        NA     NA        NA     NA 0.3569812 3.522105e-06 -1.030072
##' 
##' #### EXAMPLE 2 ####
##' logLik_mixed10(Y = 0, t = 0, beta = 1, psi = 1, tau_uu = 1, tau_uv = 0, tau_vv = 1, sigma2 = 1, lower = -5, upper = 5, integration = "univariate")
##' ##   integral1       error1  integral2       error2  integral error     logLik
##' ## 1 0.3356478 2.909504e-10 0.05056137 4.592356e-06 0.3862092    NA -0.9513762
##' 
##' logLik_mixed10(Y = 0, t = 0, beta = 1, psi = 1, tau_uu = 1, tau_uv = 0, tau_vv = 1, sigma2 = 1, lower = -5, upper = 5, integration = "bivariate")
##' ##  integral1 error1 integral2 error2  integral        error     logLik
##' ##         NA     NA        NA     NA 0.3862079 3.845569e-06 -0.9513795
##' 
##' #### EXAMPLE 3 ####
##' logLik_mixed10(Y = 0.5, t = 2.2, beta = 1, psi = 3, tau_uu = 2, tau_uv = 0.3, tau_vv = 0.25, sigma2 = 1, lower = -5, upper = 5, integration = "univariate")
##' ##    integral1       error1   integral2       error2  integral error    logLik
##' ## 1 0.09976617 1.909237e-06 0.007462216 6.232982e-09 0.1072284    NA -2.232794
##' 
##' logLik_mixed10(Y = 0.5, t = 2.2, beta = 1, psi = 3, tau_uu = 2, tau_uv = 0.3, tau_vv = 0.25, sigma2 = 1, lower = -5, upper = 5, integration = "bivariate")
##' ##   integral1 error1 integral2 error2  integral        error    logLik
##' ## 1        NA     NA        NA     NA 0.1072284 1.061953e-06 -2.232794


## * logLik_mixed10 (code)
##' @export
logLik_mixed10 <- function(Y, t,
                           beta, psi, tau_uu, tau_uv, tau_vv, sigma2,
                           lower = -5, upper = 5, indiv = TRUE, integration){

    ## ** normalize user input

    ## *** integration
    integration <- match.arg(integration, c("univariate","bivariate"))
    if(is.matrix(tau_uu)){
        Mtau <- tau_uu
        tau_vv <- Mtau[2,2]
        tau_uv <- Mtau[1,2]
        tau_uu <- Mtau[1,1]
    }else{
        Mtau <- matrix(c(tau_uu,tau_uv,tau_uv,tau_vv), nrow = 2, ncol = 2)
    }
    if((tau_uu*tau_vv - tau_uv^2)<=0){
        if(indiv==FALSE){
            warning("Variance-covariance matrix of the random effects does not seem positive definite. \n")
            return(Inf)
        }else{
            stop("Variance-covariance matrix of the random effects does not seem positive definite. \n")
        }
    }
    w_uu <- tau_vv/(tau_uu*tau_vv - tau_uv^2)
    w_uv <- -tau_uv/(tau_uu*tau_vv - tau_uv^2)
    w_vv <- tau_uu/(tau_uu*tau_vv - tau_uv^2)
    
    ## ** prepare
    n.obs <- length(Y)
    out <- data.frame(matrix(NA, nrow = n.obs, ncol = 7, dimnames = list(NULL,c("integral1","error1","integral2","error2","integral","error","logLik"))))  

    if(integration == "bivariate"){
        if(length(lower)==1){
            lower <- rep(lower, 2)
        }
        if(length(upper)==1){
            upper <- rep(upper, 2)
        }
    }else if(integration == "univariate"){
        normalizing.log <- -log(2*pi)-0.5*log(sigma2*(tau_uu*tau_vv-tau_uv^2))
        normalizing <- 1/(2*pi*sqrt(sigma2*(tau_uu*tau_vv-tau_uv^2)))
    }

    ## ** loop
    for(iId in 1:n.obs){ ## iId <- 1
        iY <- Y[iId]
        iT <- t[iId]

        if(integration == "bivariate"){
            iInt <- cubature::cubintegrate(f = function(x){ ## x <- c(1,1)
                term1 <- stats::dnorm(iY,
                                      mean = (beta+x[1])*iT*(iT<=(psi+x[2])) + (beta+x[1])*(psi+x[2])*(iT>(psi+x[2])),
                                      sd = sqrt(sigma2))
                term2 <- mvtnorm::dmvnorm(x, mean = c(0,0), sigma = Mtau)
                return(term1*term2)
            }, lower = lower, upper = upper)

            out[iId,c("integral","error","logLik")] <- c(iInt$integral,iInt$error,log(iInt$integral))
            
        }else if(integration == "univariate"){

            ## *** first term
            iInt1 <- stats::integrate(f = function(u){ ## u <- 1
                term1.log <- (iY - (beta+u)*iT)^2/(2*sigma2) + u^2*(w_uu - w_uv^2/w_vv)/2
                term2 <- (w_vv*(iT-psi)+u*w_uv)/sqrt(w_vv)
                return( exp(-term1.log) * (1 - stats::pnorm(term2)) / sqrt(w_vv) )
            }, lower = lower, upper = upper)

            ## when y=beta=t=w_uv=0 and tau_vv=tau_uu=1            
            ## sqrt(2*pi)*(1 - pnorm(-psi))*normalizing
         
            ## *** second term
            iInt2 <- stats::integrate(f = function(u){ ## u <- 1
                A <- w_vv + (beta+u)^2/sigma2
                B <- u*w_uv - (beta+u)*(iY-psi*(beta+u))/sigma2
                term1.log <- (iY - psi*(beta+u))^2/(2*sigma2) + (u^2*w_uu - B^2/A)/2
                term2 <- (A*(iT-psi)+B)/sqrt(A)
                return(exp(-term1.log) * stats::pnorm(term2)/sqrt(A) )
            }, lower = lower, upper = upper)

            ## when y=beta=t=w_uv=0 and tau_vv=tau_uu=1
            ## stats::integrate(f = function(u){ ## u <- 1
            ##     return(exp(-u^2*(1 + psi^2 + u^2)/(2*(1+u^2))) / sqrt(1 + u^2))
            ## }, lower = lower, upper = upper)$value * stats::pnorm(-psi)

            out[iId,c("integral1","error1","integral2","error2","integral","error","logLik")] <- c(iInt1$value*normalizing,
                                                                                                   iInt1$abs.error*normalizing,
                                                                                                   iInt2$value*normalizing,
                                                                                                   iInt2$abs.error*normalizing,
                                                                                                   (iInt1$value+iInt2$value)*normalizing,
                                                                                                   NA,
                                                                                                   log(iInt1$value+iInt2$value)+normalizing.log)
        }
    }

    ## ** export
    if(indiv==TRUE){
        return(out)
    }else{
        return(sum(out$logLik))
    }  
}




##----------------------------------------------------------------------
### mixedbreak.R ends here
