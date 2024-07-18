### logLik.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr 11 2024 (14:16) 
## Version: 
## Last-Updated: jul 18 2024 (16:45) 
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

        ## ** extract from object
    var.cluster <- object$args$cluster
    U.cluster <- object$args$U.cluster

    ## ** normalize user input
    ## *** cluster
    if(is.null(cluster) || identical(cluster,FALSE)){
        cluster.level <- U.cluster
    }else if(is.numeric(cluster)){
        if(cluster %in% 1:length(U.cluster) == FALSE){
            stop("When a numeric value \'cluster\' should be an integer between 1 and the number of cluster (here ",length(U.cluster),"). \n")
        }else{
            cluster.level <- U.cluster[cluster]
        }
    }else if(is.character(cluster) || is.factor(cluster)){
        cluster <- as.character(cluster) 
        if(cluster %in% U.cluster == FALSE){
            stop("When a character or factor value, \'cluster\' should be one of the strings representing the clusters (here \"",utils::head(U.cluster,1),"\", ... \"",utils::tail(U.cluster,1),"\"). \n")
        }else{
            cluster.level <- cluster
        }
    }else{
        stop("\'cluster\' should either be indexing the cluster (1 or 2 or ..., i.e. numeric) \n",
             "or the character string identifying the cluster (character or factor).")
    }

    ## ** extract logLik
    vec.logLik <- unlist(lapply(object$model,stats::logLik, ...))
    if(is.null(cluster) || identical(cluster,FALSE)){
        ## total log-likelihood
        out <- sum(vec.logLik)
    }else{
        ## individual log-likelihood
        out <- vec.logLik[match(cluster.level,U.cluster)]
    }

    ## ** export
    return(out)

}



##----------------------------------------------------------------------
### logLik.R ends here
