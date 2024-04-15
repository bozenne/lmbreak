### score.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Apr 14 2024 (10:41) 
## Version: 
## Last-Updated: Apr 15 2024 (22:52) 
##           By: Brice Ozenne
##     Update #: 42
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * .score.lmbreak
.score.lmbreak <- function(psi, transform,
                           formula, data, var.bp, var.response, dX.skeleton,
                           tol){ ## psi <- c(40,70)


    n.breakpoint <- length(psi)
    n.data <- NROW(data)
    
    ## ** transform
    if(transform){
        psi <- backtransformPsi(psi, min = attr(transform,"range")[1], max = attr(transform,"range")[2], mindiff = attr(transform,"mindiff")[1], jacobian = TRUE)
    }

    ## ** update dataset & design matrix
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
    iBetaRes <- tcrossprod(stats::residuals(iLM), stats::coef(iLM))

    ## derivative of the model matrix
    ## minus times minus -> +
    iOut <- 2 * sapply(1:n.breakpoint, function(iPoint){ ## iPoint <- 1
        iDD <- dX.skeleton[[iPoint]] * (M.dpsi[,iPoint]>0) * stats::runif(n.data, min = 1-exp(-abs(M.dpsi[,iPoint])/tol), max = 1) ## if very close to breakpoint use subgradient
        idRSS <- sum((iDD + iXXX.M1 %*% (t(iDD) %*% (iX-1) + t(iX) %*% iDD)) * iBetaRes)
        return(idRSS)
    })

    if(transform){
        return(as.double(iOut %*% attr(psi,"jacobian")))
    }else{
        return(iOut)
    }
}
##----------------------------------------------------------------------
### score.R ends here
