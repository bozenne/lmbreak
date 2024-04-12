### model.tables.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr  9 2024 (11:38) 
## Version: 
## Last-Updated: apr 11 2024 (20:00) 
##           By: Brice Ozenne
##     Update #: 38
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * model.tables.lmbreak
##' @title Extract Key Summary Statistics From Breakpoint Model
##' @description Extract key summary statistics from the breakpoint model
##'
##' @param x output of \code{\link{lmbreak}}
##' @param ... Not used. For compatibility with the generic method.
##' 
##' @keywords methods
##' @export
model.tables.lmbreak <- function(x, ...){

    ## ** normalize user input
    dots <- list(...)
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }

    ## ** extract
    allcoef <- coef(x, type = c("breakpoint.range","duration","intercept","slope"))
    allcoef$duration <- c(allcoef$duration,NA)
    allcoef$slope <- c(allcoef$slope,NA)
    out <- as.data.frame(allcoef)
    names(out)[1] <- x$args$breakpoint.var

    ## ** export
    return(out)

}

## * model.tables.mlmbreak
##' @title Extract Key Summary Statistics From Multiple Breakpoint Model
##' @description Extract key summary statistics from each breakpoint model
##'
##' @param x output of \code{\link{mlmbreak}}.
##' @param cluster [vector] cluster relative to which the summary statistics should be extracted..
##' @param format [character] should the output be a data.frame (with cluster as a column) or an array (with cluster as the third dimension)?
##' @param ... Not used. For compatibility with the generic method.
##' 
##' @keywords methods
##' @export
model.tables.mlmbreak <- function(x, cluster = NULL, format = "data.frame", ...){

    ## ** extract from object
    var.cluster <- x$args$cluster
    U.cluster <- x$args$U.cluster

    ## ** normalize user input
    if(is.null(cluster)){
        cluster <- U.cluster
    }else if(any(cluster %in% U.cluster == FALSE)){
        stop("Unknown value for argument \'cluster\'.")
    }
    format <- match.arg(format, c("data.frame","array","list"))

    ## ** extract
    ls.table <- lapply(cluster, function(iC){ ## iC <- cluster[6]
        iOut <- model.tables(as.lmbreak(x, cluster = iC))
        if(format == "data.frame"){
            return(cbind(iC,iOut))
        }else{
            return(iOut)
        }
    })    

    ## ** export
    if(format == "data.frame"){
        out <- do.call(rbind,ls.table)
        names(out)[1] <- var.cluster
    }else if(format == "array"){
        if(length(unique(sapply(ls.table,NROW)))>1){
            stop("Cannot convert the output in an array.\n",
                 "Not all clusters have the same number of breakpoints. \n")
        }
        out <- array(unlist(ls.table), dim = c(dim(ls.table[[1]]), length(cluster)),
                     dimnames = c(dimnames(ls.table[[1]]), list(cluster)))
    }else if(format == "list"){
        out <- ls.table
        names(out) <- cluster
    }
    return(out)

}

##----------------------------------------------------------------------
### model.tables.R ends here
