### utils.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr 12 2024 (10:19) 
## Version: 
## Last-Updated: Apr 20 2024 (18:39) 
##           By: Brice Ozenne
##     Update #: 63
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * combn2
##' @title Generate All Combination of m Elements with n at a Time
##' @description same as utils::combn but handle certain special case instead of returning an error
##' called by \code{lmbreak} when initializing with gam or quantiles
##' @param x [vector] vector source for combinations, or integer.
##' @param m [integer,>0] number of elements to choose.
##' @param name [character vector] name for the columns.
##' Should either has length 1 (same name for all columns) or length the number of possibilities.
##' @noRd
##' @examples
##' combn2(NULL,2)
##' combn2(1:5,2)
##' combn2(1:5,2, name = "auto")

combn2 <- function(x, m, name = NULL){
    if(is.null(x)||length(x)<m){
        out <- NULL
    }else{
        out <- utils::combn(x=x, m=m)
        if(!is.null(name)){
            if(length(name)==1){
                name <- rep(name,NCOL(out))
            }else if(length(name)!=NCOL(out)){
                stop("Incorrect length for argument \'name\': should have length ",NCOL(out),"\n")
            }
            colnames(out) <- name
        }
    }

    return(out)
}

## * transformFree
##' @title Transformation to an unconstrained scale
##' @description Transform a value defined on an interval to a value defined on the real Line.
##' 
##' @param value [numeric vector] values to be transformed.
##' @param min [numeric] minimum value. Can be \code{-Inf}.
##' @param max [numeric] maximum value. Can be \code{Inf}.
##' @param jacobian [logical] should the jacobian of the transformation be output.
##' 
##' @noRd
##' @examples
##' res <- transformFree(1:99, min = 0, max = 100, jacobian = TRUE)
##' attr(res, "jacobian") - numDeriv::grad(1:99, func = transI2R, min = 0, max = 100)
##' 
##' resInf <- transformFree(6:10, min = 5, max = Inf, jacobian = TRUE)
##' attr(resInf, "jacobian") - numDeriv::grad(6:10, func = transI2R, min = 5, max = Inf)
##' 
##' resMinf <- transformFree(-1:5, min = -Inf, max = 10, jacobian = TRUE)
##' attr(resMinf, "jacobian") - numDeriv::grad(-1:5, func = transI2R, min = -Inf, max = 10)
transformFree <- function(value, min, max, jacobian = FALSE){

    ## ** normalize user input
    if(any(is.na(value))){
        stop("Argument \'value\' cannot be NA. \n")
    }else if(any(is.infinite(value))){
        stop("Argument \'value\' cannot be infinite. \n")
    }

    ## ** transformation
    if(!is.infinite(min) && !is.infinite(max)){
        out <- atanh(2*(value - min)/(max - min) - 1)
        if(jacobian){
            attr(out,"jacobian") <- 2/((1 - (2*(value - min)/(max - min) - 1)^2)*(max - min))
        }
    }else if(!is.infinite(min)){
        out <- log(value - min)
        if(jacobian){
            attr(out,"jacobian") <- 1/(value - min)
        }
    }else if(!is.infinite(max)){
        out <- log(-value + max)
        if(jacobian){
            attr(out,"jacobian") <- -1/(-value + max)
        }
    }else{
        out <- value
        if(jacobian){
            attr(out,"jacobian") <- rep(1, length(out))
        }
    }
    return(out)

}

## * backtransformFree
##' @title Backtransformation from an unconstrained scale
##' @description Backtransform a value defined on an interval form a value defined on the real Line.
##' @noRd
##' @examples
##' res <- transformFree(1:99, min = 0, max = 100)
##' ires <- backtransformFree(res, min = 0, max = 100, jacobian = TRUE)
##' range(ires - 1:99)
##' range(attr(ires, "jacobian") - numDeriv::grad(res, func = transR2I, min = 0, max = 100))
##' 
##' resInf <- transformFree(6:10, min = 5, max = Inf)
##' iresInf <- backtransformFree(resInf, min = 5, max = Inf, jacobian = TRUE)
##' range(iresInf - 6:10)
##' range(attr(iresInf, "jacobian") - numDeriv::grad(resInf, func = transR2I, min = 5, max = Inf))
##' 
##' resMinf <- transformFree(-1:5, min = -Inf, max = 10)
##' iresMinf <- backtransformFree(resMinf, min = -Inf, max = 10, jacobian = TRUE)
##' range(iresMinf - (-1:5))
##' range(attr(iresMinf, "jacobian") - numDeriv::grad(resMinf, func = transR2I, min = -Inf, max = 10))
backtransformFree <- function(value, min, max, jacobian = FALSE){

    ## ** normalize user input
    if(any(is.na(value))){
        stop("Argument \'value\' cannot be NA. \n")
    }

    ## ** backtransformation
    if(!is.infinite(min) && !is.infinite(max)){
        out <- (tanh(value) + 1)*(max - min)/2
        if(jacobian){
            attr(out,"jacobian") <- (1-tanh(value)^2)*(max - min)/2
        }
    }else if(!is.infinite(min)){
        out <- exp(value) + min
        if(jacobian){
            attr(out,"jacobian") <- exp(value)
        }
    }else if(!is.infinite(max)){
        out <- -(exp(value) - max)
        if(jacobian){
            attr(out,"jacobian") <- -exp(value)
        }
    }else{
        out <- value
        if(jacobian){
            attr(out,"jacobian") <- rep(1, length(out))
        }
    }
    return(out)

}

##----------------------------------------------------------------------
### utils.R ends here
