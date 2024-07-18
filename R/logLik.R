### logLik.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr 11 2024 (14:16) 
## Version: 
## Last-Updated: jul 18 2024 (10:30) 
##           By: Brice Ozenne
##     Update #: 40
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
logLik.lmbreak <- function(object, continuity = NULL, ...){

    ## ** normalize user input
    ## *** dots
    dots <- list(...)
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }

    ## *** continuity
    if(is.null(continuity)){
        continuity <- (object$opt$continuity==FALSE)
    }

    ## ** extract from object
    if(continuity && !is.null(attr(object$model,"continuity"))){
        model <- attr(object$model,"continuity")
    }else{
        model <- object$model
    }

    out <- as.double(stats::logLik(model))

    ## ** export
    return(out)

}

## * logLik.mlmbreak
##' @export
logLik.mlmbreak <- function(object, cluster = FALSE, ...){

    out <- unlist(lapply(object$model,stats::logLik, ...))
    if(identical(cluster,FALSE)){
        out <- sum(out)
    }else if(all(cluster %in% object$args$U.cluster)){
        out <- out[match(cluster,object$args$U.cluster)]
    }
    return(out)

}



##----------------------------------------------------------------------
### logLik.R ends here
