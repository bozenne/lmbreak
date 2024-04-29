### hessian.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Apr 14 2024 (20:18) 
## Version: 
## Last-Updated: Apr 20 2024 (09:12) 
##           By: Brice Ozenne
##     Update #: 55
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * .hessian.lmbreak
.hessian.lmbreak <- function(psi, transform,
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

    ## model fit
    iLM <- stats::lm(formula, data = data)
    iXX.M1 <- stats::vcov(iLM) / stats::sigma(iLM)^2
    iX <- stats::model.matrix(iLM)
    iXXX.M1 <- iX %*% iXX.M1
    iBeta <- stats::coef(iLM)
    iRes <- stats::residuals(iLM)
    iBetaRes <- tcrossprod(iRes, iBeta)

    ## diagonal terms
    iOut <- matrix(NA, nrow = n.breakpoint, ncol = n.breakpoint)
    ls.term3.score <- vector(mode = "list", length = n.breakpoint)
    ls.term123.score <- vector(mode = "list", length = n.breakpoint)
    ls.dterm2.score <- vector(mode = "list", length = n.breakpoint)
    ls.dterm3.score <- vector(mode = "list", length = n.breakpoint)
    ls.dbetares <- vector(mode = "list", length = n.breakpoint)

    for(iPoint in 1:n.breakpoint){ ## iPoint <- 2
        ## score sum((iDD + iXXX.M1 %*% (t(iDD) %*% (iX-1) + t(iX) %*% iDD)) * iBetaRes)
        ## term1: iDD
        ## term2: iXXX.M1
        ## term3: t(iDD) %*% (iX-1) + t(iX) %*% iDD)

        ## derivative of the contrast matrix
        iDD <- dX.skeleton[[iPoint]] * (M.dpsi[,iPoint]>0) * stats::runif(n.data, min = 1-exp(-abs(M.dpsi[,iPoint])/tol), max = 1)
        iDD2 <- t(iDD) %*% iX + t(iX) %*% iDD
        ## score terms (not yet defined, i.e. involving derivatives)
        iTerm3.score <- t(iDD) %*% (iX-1) + t(iX) %*% iDD
        iTerm123.score <- iDD + iXXX.M1 %*% iTerm3.score
        ## derivative of the score terms (except iBetaRes)
        iDD.term2.score <- (- iDD + iXXX.M1 %*% iTerm3.score) %*% iXX.M1
        iDD.term3.score <- (t(iDD) %*% (iDD-1) + t(iDD) %*% iDD)
        ## derivative of residuals times regression coefficients
        iDD.XXX.M1 <- iXX.M1 %*% (iDD2 %*% iXX.M1 %*% t(iX) - t(iDD)) 
        iDD.BetaRes <- tcrossprod(iDD %*% iBeta - iX %*% iDD.XXX.M1 %*% data[[var.response]], iBeta) + tcrossprod(iRes, iDD.XXX.M1 %*% data[[var.response]])
        ## assemble
        iOut[iPoint,iPoint] <- 2 * sum( (iDD.term2.score %*% iTerm3.score + iXXX.M1 %*% iDD.term3.score) * iBetaRes + iTerm123.score * iDD.BetaRes)

        ls.term3.score[[iPoint]] <- iTerm3.score
        ls.term123.score[[iPoint]] <- iTerm123.score
        ls.dterm2.score[[iPoint]] <- iDD.term2.score
        ls.dterm3.score[[iPoint]] <- iDD.term3.score
        ls.dbetares[[iPoint]] <- iDD.BetaRes
        ## term1 <- - 2 * sum((iDD %*% iXX.M1 %*% (t(iDD) %*% (iX-1) + t(iX) %*% iDD)) * iBetaRes)
        ## term2 <- 2 * sum((iX %*% iXX.M1 %*% (t(iDD) %*% (iX-1) + t(iX) %*% iDD) %*% iXX.M1 %*% (t(iDD) %*% (iX-1) + t(iX) %*% iDD)) * iBetaRes)
        ## term3 <- - 2 * sum((iX %*% iXX.M1 %*% (t(iDD) %*% (iDD-1) + t(iDD) %*% iDD)) * iBetaRes)
        ## term4 <- 2 * sum((iDD + iX %*% iXX.M1 %*% (t(iDD) %*% (iX-1) + t(iX) %*% iDD)) * iDD.BetaRes)
        ## term1 + term2 + term3 + term4
    }

    ## off-diagonal terms
    if(n.breakpoint>1){
        for(iPoint1 in 2:n.breakpoint){ ## iPoint1 <- 2
            for(iPoint2 in 1:(iPoint1-1)){ ## iPoint2 <- 1
                iOut[iPoint1,iPoint2] <- 2 * sum( (ls.dterm2.score[[iPoint2]] %*% ls.term3.score[[iPoint1]] + iXXX.M1 %*% ls.dterm3.score[[iPoint2]]) * iBetaRes + ls.term123.score[[iPoint1]] * ls.dbetares[[iPoint2]])
                iOut[iPoint2,iPoint1] <- iOut[iPoint1,iPoint2]
                ## 2 * sum( (ls.dterm2.score[[iPoint1]] %*% ls.term3.score[[iPoint2]] + iXXX.M1 %*% ls.dterm3.score[[iPoint1]]) * iBetaRes + ls.term123.score[[iPoint2]] * ls.dbetares[[iPoint1]])            
            }
        }
    }

    if(transform){
        return(as.double(iOut %*% attr(psi,"jacobian")))
    }else{
        return(iOut)
    }
}


##----------------------------------------------------------------------
### hessian.R ends here
