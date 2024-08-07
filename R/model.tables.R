### model.tables.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr  9 2024 (11:38) 
## Version: 
## Last-Updated: jul 18 2024 (16:53) 
##           By: Brice Ozenne
##     Update #: 46
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
##' @param continuity [logical] should coefficients be extracted from a breakpoint model ensuring continuity (i.e. no Vs terms)?
##' Often not relevant as Vs term should be 0 when proper convergence has been reached.
##' @param ... Not used. For compatibility with the generic method.
##' 
##' @keywords methods
##' @export
model.tables.lmbreak <- function(x, continuity = NULL, ...){

    ## ** normalize user input
    dots <- list(...)
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }

    ## ** extract
    allcoef <- coef(x, type = c("breakpoint.range","duration","intercept","slope"), continuity = continuity)
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
##' @param cluster [numeric or character] cluster relative to which the summary statistics should be extracted.
##' @param format [character] should the output be a data.frame (with cluster as a column) or an array (with cluster as the third dimension)?
##' @param ... Passed to \code{model.tables.lmbreak}.
##' 
##' @keywords methods
##' @export
model.tables.mlmbreak <- function(x, cluster = NULL, format = "data.frame", ...){

    ## ** extract from object
    var.cluster <- x$args$cluster
    U.cluster <- x$args$U.cluster

    ## ** normalize user input
    ## *** cluster
    if(is.null(cluster)){
        cluster.level <- U.cluster
    }else if(is.numeric(cluster)){
        if(any(cluster %in% 1:length(U.cluster) == FALSE)){
            stop("When a numeric value \'cluster\' should be an integer between 1 and the number of cluster (here ",length(U.cluster),"). \n")
        }else{
            cluster.level <- U.cluster[cluster]
        }
    }else if(is.character(cluster) || is.factor(cluster)){
        cluster <- as.character(cluster) 
        if(any(cluster %in% U.cluster == FALSE)){
            stop("When a character or factor value, \'cluster\' should be one of the strings representing the clusters (here \"",utils::head(U.cluster,1),"\", ... \"",utils::tail(U.cluster,1),"\"). \n")
        }else{
            cluster.level <- cluster
        }
    }else{
        stop("\'cluster\' should either be indexing the cluster (1 or 2 or ..., i.e. numeric) \n",
             "or the character string identifying the cluster (character or factor).")
    }

    ## *** format
    format <- match.arg(format, c("data.frame","array","list"))

    ## ** extract
    ls.table <- lapply(cluster.level, function(iC){ ## iC <- cluster[6]
        iOut <- model.tables(as.lmbreak(x, cluster = iC), ...)
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
        out <- array(unlist(ls.table), dim = c(dim(ls.table[[1]]), length(cluster.level)),
                     dimnames = c(dimnames(ls.table[[1]]), list(cluster.level)))
    }else if(format == "list"){
        out <- ls.table
        names(out) <- cluster.level
    }
    return(out)

}

##----------------------------------------------------------------------
### model.tables.R ends here
