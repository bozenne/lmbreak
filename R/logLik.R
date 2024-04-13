### logLik.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr 11 2024 (14:16) 
## Version: 
## Last-Updated: apr 13 2024 (20:20) 
##           By: Brice Ozenne
##     Update #: 11
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
##----------------------------------------------------------------------
### logLik.R ends here
