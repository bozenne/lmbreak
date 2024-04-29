### score.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Apr 14 2024 (10:41) 
## Version: 
## Last-Updated: Apr 20 2024 (18:56) 
##           By: Brice Ozenne
##     Update #: 105
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * score.lmbreak (documentation)
##' @title Breakpoint Score Equation
##' @description First derivative of the residual sum of squares with respect to the breakpoints.
##'
##' @param x [lmbreak] breakpoint model
##' @param psi [numeric vector] breakpoint values.
##' @param indiv [logical] should the contribution of each observation to the score be output.
##' Otherwise outputs the sum of the contributions.
##' @param ... Not used. For compatibility with the generic method.
##'
##' @return A numeric vector (\code{indiv=FALSE}) where each element corresponds to a breakpoint.
##' or a matrix (\code{indiv=TRUE}) where each row corresponds to an individual and each column to a breakpoint.
##' 
##' @keywords methods
##' 
##' @examples
##' if(require(lava)){
##' set.seed(10)
##' df <- simBreak(c(1,100), breakpoint = c(0,50,100), slope = c(1,-1))
##' score(lmbreak(Y ~ bp(X, "11"), data = df))
##' }

## * score.lmbreak (code)
##' @export
score.lmbreak <- function(x, indiv = FALSE, psi = NULL, ...){

    ## ** normalize user input
    ## psi
    if(is.null(psi)){
        psi <- x$breakpoint$value
    }else{ 
        if(!is.numeric(psi)){
            stop("Argument \'psi\' must be a numeric vector. \n")
        }
        if(length(psi) != NROW(x$breakpoint)){
            stop("Argument \'psi\' must be a vector of length ",NROW(x$breakpoint),". \n")
        }
    }

    ## data
    if(length(x$index.NA)==0){
        data <- x$data
    }else{
        data <- x$data[-x$index.NA,,drop=FALSE]
    }

    ## ** evalate the score
    out <- .score.lmbreak(psi, indiv = indiv, transform = FALSE,
                          formula = attr(x$formula.pattern,"noVs"), data = data, var.bp = x$args$breakpoint.var, var.response = x$args$response.var,
                          dX.skeleton = NULL, tol = x$opt$tol[[1]][1])

    ## ** export
    return(out)

}

## * .score.lmbreak
.score.lmbreak <- function(psi, indiv, transform,
                           formula, data, var.bp, var.response, dX.skeleton,
                           tol){ ## psi <- c(40,70)


    n.breakpoint <- length(psi)
    n.data <- NROW(data)
    
    ## ** transform
    if(transform){
        psi <- backtransformPsi(psi, min = attr(transform,"range")[1], max = attr(transform,"range")[2], mindiff = attr(transform,"mindiff")[1], jacobian = TRUE)
    }

    ## ** update dataset & design matrix
    if(is.null(dX.skeleton)){
        dX.skeleton <-lapply(1:n.breakpoint, function(iPoint){ ## iPoint <- 2
            ddata <- data
            ddata[,paste0("Us",0:n.breakpoint)] <- matrix(as.numeric(0:n.breakpoint == iPoint), byrow = TRUE, nrow = NROW(data), ncol = n.breakpoint+1)
            iX <- stats::model.matrix(formula, ddata)
            if("(Intercept)" %in% colnames(iX)){
                iX[,"(Intercept)"] <- 0
            }
            return(iX)
        })
    }
    if("Us0" %in% names(data) == FALSE){
        data$Us0 <- data[[var.bp]]
    }
    M.dpsi <- matrix(NA, nrow = n.data, ncol = n.breakpoint)
    for(iPoint in 1:n.breakpoint){ ## iPoint <- 1
        M.dpsi[,iPoint] <- data[[var.bp]] - psi[iPoint]
        data[[paste0("Us",iPoint)]] <- M.dpsi[,iPoint]*(M.dpsi[,iPoint]>0)
    }

    ## ** model fit
    iLM <- stats::lm(formula, data = data)
    iXX.M1 <- stats::vcov(iLM) / stats::sigma(iLM)^2
    iY <- data[[var.response]]
    iX <- stats::model.matrix(iLM)
    if(indiv){
        iXXX.M1 <- iX %*% iXX.M1
        iXXX.M1_XmI <- iXXX.M1 %*% t(iX) - diag(1,n.data,n.data)
        itXXX.M1 <- t(iXXX.M1)
    }else{
        iYXXX.M1 <- (iY %*% iX) %*% iXX.M1
        iXXXY.M1_XmI <- iYXXX.M1 %*% t(iX) - iY
        itYXXX.M1 <- t(iYXXX.M1)
    }
    
    ## ** derivative
    ## RSS = \sum_i (Yi-Xi\beta)^2 = (Y - X(XX)^-1XY)^2 = Y(I-H)(I-H)Y = Y(I-H)Y
    ## d{RSS}/d\psi = - Y d{H}/d\psi Y
    ##              = - Y d{X(psi)(X(psi)X(psi))^{-1}X(psi)}/d\psi Y
    ## 
    ## dX = d{X(\psi)}/d{psi} = d{(X-psi)\Ind[X>psi]}/d{psi} = - \Ind[X>psi] or an element between [0,1] if X==psi
    ## (With plateau the design matrix would be Z = CX with deterministic C, e.g. [(X-psi[1]) - (X-psi[3]), (X-psi[2]) - (X-psi[3])].
    ##  Then the derivative is dZ = CdX )
    ## 
    ## term1 = - dX (t(X)X)^{-1}X
    ## term2 = X (t(X)X)^{-1} (dX X + X dX) (t(X)X)^{-1}X
    ## term3 = - X (t(X)X)^{-1}dX = t(term1)
    ls.iOut <- lapply(1:n.breakpoint, function(iPoint){ ## iPoint <- 1
        iDD <- - dX.skeleton[[iPoint]] * (M.dpsi[,iPoint]>0) * stats::runif(n.data, min = 1-exp(-abs(M.dpsi[,iPoint])/tol), max = 1) ## if very close to breakpoint use subgradient
        ## ## version 1
        ## term1 <- - iDD %*% itXXX.M1
        ## term2 <- iXXX.M1 %*% (t(iDD) %*% iX + t(iX) %*% iDD) %*% itXXX.M1
        ## iScore <- (iY %*% (term1 + term2 + t(term1))) * iY        
        if(indiv){
            ## ## version 2 (faster)
            term <- iXXX.M1_XmI %*% iDD %*% itXXX.M1
            iScore <- iY * (term + t(term)) %*% iY
        }else{
            ## ## version 3 (even faster)
            iScore <- 2 * iXXXY.M1_XmI %*% iDD %*% itYXXX.M1
        }
        return(iScore)
    })
    if(indiv){
        iOut <- do.call(cbind,ls.iOut)
    }else{
        iOut <- unlist(ls.iOut)
    }

    if(transform){
        return(as.double(iOut %*% attr(psi,"jacobian")))
    }else{        
        return(iOut)
    }
}
##----------------------------------------------------------------------
### score.R ends here
