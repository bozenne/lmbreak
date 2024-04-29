### RRS.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Apr 20 2024 (09:27) 
## Version: 
## Last-Updated: Apr 20 2024 (18:54) 
##           By: Brice Ozenne
##     Update #: 25
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * anova.lmbreak (documentation)
##' @title Residual Sum of Squares
##' @description Evaluate the residual sum of squares (RSS) at given breakpoints.
##'
##' @param object [lmbreak] breakpoint model
##' @param psi [numeric vector] breakpoint values.
##' @param indiv [logical] should the contribution of each observation to the RSS be output.
##' Otherwise outputs the sum of the contributions.
##' @param ... Not used. For compatibility with the generic method.
##'
##' @return A numeric value (\code{indiv=FALSE}) or a numeric vector (\code{indiv=TRUE}) where each element corresponds to an individual.
##' 
##' @keywords methods
##' 
##' @examples
##' set.seed(10)
##' df <- simBreak(c(1,100), breakpoint = c(0,50,100), slope = c(1,-1))
##' anova(lm(Y ~ X, data = df))
##' 
##' anova(lmbreak(Y ~ bp(X, "11"), data = df))


## * anova.lmbreak (code)
##' @export
anova.lmbreak <- function(object, indiv = FALSE, psi = NULL, ...){

    ## ** normalize user input
    if(is.null(psi)){
        psi <- object$breakpoint$value
    }else{ 
        if(!is.numeric(psi)){
            stop("Argument \'psi\' must be a numeric vector. \n")
        }
        if(length(psi) != NROW(object$breakpoint)){
            stop("Argument \'psi\' must be a vector of length ",NROW(object$breakpoint),". \n")
        }
    }

    ## data
    if(length(object$index.NA)==0){
        data <- object$data
    }else{
        data <- object$data[-object$index.NA,,drop=FALSE]
    }


    ## ** evalate residual sum of squares
    out <- .RSS.lmbreak(psi, indiv = indiv, transform = FALSE,
                        formula = attr(object$formula.pattern,"noVs"), data = data, var.bp = object$args$breakpoint.var, var.response = object$args$response.var)

    ## ** export
    return(out)

}


## * .RSS.lmbreak
.RSS.lmbreak <- function(psi, indiv, transform,
                         formula, data, var.bp, var.response){


    n.breakpoint <- length(psi)

    ## ** transform
    if(transform){
        psi <- backtransformPsi(psi, min = attr(transform,"range")[1], max = attr(transform,"range")[2], mindiff = attr(transform,"mindiff")[1], jacobian = FALSE)
    }

    ## ** update dataset & design matrix
    if("Us0" %in% names(data) == FALSE){
        data$Us0 <- data[[var.bp]]
    }
    for(iPoint in 1:n.breakpoint){ ## iPoint <- 1
        data[[paste0("Us",iPoint)]] <- (data[[var.bp]] - psi[iPoint])*(data[[var.bp]] > psi[iPoint])
    }
    iX <- stats::model.matrix(formula, data = data)

    ## ** evaluate residual sum of squares (RSS)
    if(det(crossprod(iX))>0){
        iOut <- stats::lm.fit(y = data[[var.response]], x = iX)$residuals^2
    }else{
        iOut <- data[[var.response]]^2
    }

    ## ** export
    if(indiv){
        return(iOut)
    }else{
        return(sum(iOut))
    }
}

##----------------------------------------------------------------------
### RRS.R ends here
