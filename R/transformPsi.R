### transformPsi.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Apr 14 2024 (11:28) 
## Version: 
## Last-Updated: jul 18 2024 (14:16) 
##           By: Brice Ozenne
##     Update #: 21
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * transformPsi & backtransformPsi (documentation)
##' @description Transform breakpoints to an unconstrained scale.
##' 
##' @param psi [numeric vector] brekapoints.
##' @param min [numeric] small value in the dataset.
##' @param max [numeric] largest value in the dataset.
##' @param mindiff [numeric] small difference between two observations in the dataset.
##' @param jacobian [logical] should the jacobian of the transformation be output?n
##' 
##' @noRd
##' @examples
##'
##' #### single breakpoint ####
##' res1 <- transformPsi(50, min = 0.5, max = 99.2, mindiff = 2, jacobian = TRUE)
##' ires1 <- backtransformPsi(res1, min = 0.5, max = 99.2, mindiff = 2, jacobian = TRUE)
##' attr(res1,"jacobian") * attr(ires1,"jacobian") 
##' 
##' #### two breakpoints ####
##' res2 <- transformPsi(c(50,75), min = 0.5, max = 99.2, mindiff = 2, jacobian = TRUE)
##' attr(res2,"jacobian") - numDeriv::jacobian(c(50,75), func = transformPsi, min = 0.5, max = 99.2, mindiff = 2)
##' ires2 <- backtransformPsi(res2, min = 0.5, max = 99.2, mindiff = 2, jacobian = TRUE)
##' ires2
##' attr(ires2, "jacobian") - numDeriv::jacobian(res2, func = backtransformPsi, min = 0.5, max = 99.2, mindiff = 2)
##' 
##' #### three breakpoints ####
##' res3 <- transformPsi(c(20,50,75), min = 0.5, max = 99.2, mindiff = 2, jacobian = TRUE)
##' attr(res3,"jacobian") - numDeriv::jacobian(c(20,50,75), func = transformPsi, min = 0.5, max = 99.2, mindiff = 2)
##' 
##' ires3 <- backtransformPsi(res3, min = 0.5, max = 99.2, mindiff = 2, jacobian = TRUE)
##' ires3
##' attr(ires3, "jacobian") - numDeriv::jacobian(res3, func = backtransformPsi, min = 0.5, max = 99.2, mindiff = 2)

## * transformPsi (code)
transformPsi <- function(psi, min, max, mindiff, jacobian = FALSE){

    ## ** normalize user input
    min2 <- min + mindiff/2
    max2 <- max - mindiff/2
    n.breakpoint <- length(psi)

    ## ** transformation
    out <- rep(NA, n.breakpoint)
    if(jacobian){
        attr(out, "jacobian") <- matrix(0, n.breakpoint, n.breakpoint)
    }

    out1 <- transformFree(psi[1], min = min2, max = max2, jacobian = jacobian)
    out[1] <- out1
    if(jacobian){
        attr(out, "jacobian")[1,1] <- attr(out1, "jacobian")
    }

    if(n.breakpoint>1){
        out[-1] <- log(diff(psi))
        if(jacobian){
            for(iBreakpoint in 2:n.breakpoint) {
                attr(out, "jacobian")[iBreakpoint,c(iBreakpoint-1,iBreakpoint)] <- c(-1,1)*1/diff(psi)
            }
        }
    }

    ## ** export
    return(out)
}

## * backtransformPsi (code)
backtransformPsi <- function(psi, min, max, mindiff, jacobian = FALSE){

    ## ** normalize user input
    min2 <- min + mindiff/2
    max2 <- max - mindiff/2
    n.breakpoint <- length(psi)

    ## ** transformation
    out <- rep(NA, n.breakpoint)
    if(jacobian){
        attr(out, "jacobian") <- matrix(0, n.breakpoint, n.breakpoint)
    }

    out1 <- backtransformFree(psi[1], min = min2, max = max2, jacobian = jacobian)
    out[1] <- out1
    if(jacobian){
        attr(out, "jacobian")[1,1] <- attr(out1, "jacobian")
    }

    if(n.breakpoint>1){
        out[-1] <- out[1] + cumsum(exp(psi[-1]))
        if(jacobian){
            for(iBreakpoint in 2:n.breakpoint) {
                attr(out, "jacobian")[iBreakpoint,1:iBreakpoint] <- c(attr(out, "jacobian")[1,1],exp(psi[2:iBreakpoint]))
            }
        }
    }

    ## ** export
    return(out)
}

##----------------------------------------------------------------------
### transformPsi.R ends here
