### as.lmbreak.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr 10 2024 (14:05) 
## Version: 
## Last-Updated: jul 18 2024 (16:37) 
##           By: Brice Ozenne
##     Update #: 21
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * as.lmbreak (documentation)
##' @title From Multiple Breakpoint models to a Single Breakpoint Model
##' @description Retrive the output of \code{\link{lmbreak}} from \code{\link{mlmbreak}} for a specific cluster.
##'
##' @param object output of \code{\link{mlmbreak}}.
##' @param cluster [numeric or character] cluster relative to which the breakpoint model should be converted.
##'
##' @keywords utilities

## * as.lmbreak (code)
##' @export
as.lmbreak <- function(object, cluster){

    ## ** normalize user input
    if(!inherits(object,"mlmbreak")){
        stop("Argument \'object\' must inherits from mlmbreak. \n")
    }

    ## *** cluster
    if(length(cluster)!=1){
        stop("Argument \'cluster\' should have length 1.")
    }
    U.cluster <- object$args$U.cluster
    if(is.numeric(cluster)){
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

    ## ** re-generate object
    out <- list(model = object$model[[cluster]],
                breakpoint = object$breakpoint[[cluster]],
                phase = object$phase[[cluster]],
                opt = object$opt[[cluster]],
                call = object$call,
                args = object$args,
                data = object$data[as.character(object$data[[object$args$cluster]])==cluster.level,,drop=FALSE])
    out$call$cluster <- paste0(out$call$cluster,"==",cluster)
    class(out) <- "lmbreak"

    ## ** export
    return(out)
}


##----------------------------------------------------------------------
### as.lmbreak.R ends here
