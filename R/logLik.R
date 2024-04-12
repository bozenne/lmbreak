### logLik.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr 11 2024 (14:16) 
## Version: 
## Last-Updated: apr 11 2024 (14:29) 
##           By: Brice Ozenne
##     Update #: 5
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
logLik.mlmbreak <- function(object,...){

    ls.logLik <- lapply(object$model,stats::logLik)
    out <- sum(unlist(ls.logLik))
    return(out)

}
##----------------------------------------------------------------------
### logLik.R ends here
