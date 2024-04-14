### logLik.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr 11 2024 (14:16) 
## Version: 
## Last-Updated: Apr 14 2024 (16:22) 
##           By: Brice Ozenne
##     Update #: 34
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * logLik.lmbreak
##' @export
logLik.lmbreak <- function(object,...){

    out <- as.double(stats::logLik(object$model))
    return(out)

}

## * logLik.mlmbreak
##' @export
logLik.mlmbreak <- function(object, cluster = FALSE,...){

    out <- unlist(lapply(object$model,stats::logLik))
    if(identical(cluster,FALSE)){
        out <- sum(out)
    }else if(all(cluster %in% object$args$U.cluster)){
        out <- out[match(cluster,object$args$U.cluster)]
    }
    return(out)

}


## * .RRS.lmbreak
.RRS.lmbreak <- function(psi, transform,
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
        iOut <- sum(stats::lm.fit(y = data[[var.response]], x = iX)$residuals^2)
    }else{
        iOut <- sum(data[[var.response]]^2)
    }

    ## ** export
    return(iOut)
}

##----------------------------------------------------------------------
### logLik.R ends here
