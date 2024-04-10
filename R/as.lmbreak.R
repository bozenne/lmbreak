### as.lmbreak.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr 10 2024 (14:05) 
## Version: 
## Last-Updated: apr 10 2024 (16:42) 
##           By: Brice Ozenne
##     Update #: 14
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
##' @param cluster [numeric or character] cluster for which the breakpoint model should be converted.
##'
##' @keywords utilities

## * as.lmbreak (code)
##' @export
as.lmbreak <- function(object, cluster){

    ## ** normalize user input
    if(!inherits(object,"mlmbreak")){
        stop("Argument \'object\' must inherits from mlmbreak. \n")
    }
    if(length(cluster)!=1){
        stop("Argument \'cluster\' should have length 1.")
    }

    index.cluster <- which(cluster == object$args$U.cluster)
    if(length(index.cluster)!=1){
        stop("Unknown value for argument \'cluster\'.")
    }

    ## ** re-generate object
    out <- list(model = object$model[[index.cluster]],
                breakpoint = object$breakpoint[object$breakpoint[[object$args$cluster]]==cluster,,drop=FALSE],
                opt = object$opt[index.cluster,,drop=FALSE],
                call = object$call,
                args = object$args,
                data = object$data[object$data[[object$args$cluster]]==cluster,,drop=FALSE])
    attr(out$breakpoint,"all") <- attr(object$breakpoint,"all")[[cluster]]
    attr(out$opt,"all") <- attr(object$opt,"all")[[cluster]]
    out$call$cluster <- paste0(out$call$cluster,"==",cluster)
    class(out) <- "lmbreak"

    ## ** export
    return(out)
}


##----------------------------------------------------------------------
### as.lmbreak.R ends here
