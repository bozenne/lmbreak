### mlmbreak.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr  9 2024 (12:22) 
## Version: 
## Last-Updated: apr 13 2024 (19:41) 
##           By: Brice Ozenne
##     Update #: 95
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * mlmbreak (documentation)
##' @title Fit Multiple Breakpoint models
##' @description Fit a linear regression with breakpoints separately on each cluster of data.
##' Observations are independent between cluster but can be correlated within cluster.
##' @name mlmbreak
##'
##' @param formula a formula where the breakpoint variable appears on the right hand side
##' with an argument pattern specifying the number of breakpoints and possible constrains. See details section of \code{\link{lmbreak}}.
##' @param data [data.frame] dataset
##' @param cluster [character] variable in the dataset indicating the cluster level.
##' @param trace [0,1,2,3] trace the execution of the function.
##' @param cpus [integer, >0] the number of CPU to use. When greater than 1, calculations are done parallel for different clusters. 
##' @param ... additional arguments passed to \code{\link{lmbreak}}.
##' 
##' @seealso
##' \code{\link{as.lmbreak}} to extract the results specific to a single cluster of data. \cr
##' \code{\link{lmbreak}} for fitting the breakpoint model on a single cluster of data. \cr
##' \code{\link{model.tables.mlmbreak}} for extracting the breakpoints, slopes, intercept, and duration. \cr
##' \code{\link{plot.mlmbreak}} for a graphical display of the fitted breakpoint model. \cr
##'
##' @keywords models
##' 
##' @examples
##' 
##' ####  simulate data ####
##' set.seed(10)
##' df1 <- simBreak(c(10,25), breakpoint = c(0,1.5,3,4,5.5), slope = c(1,0,-1.2,0), sigma = 0.05)
##'
##' #### fit breakpoint regression ####
##' e.mlmbreak1010 <- mlmbreak(Y ~ 0 + bp(X, "1010"), cluster = "id", data = df1)
##' ## expect 2 warnings
##' plot(e.mlmbreak1010, scales = "free")
##' plot(e.mlmbreak1010, scales = "free", subtitle = 0)
##'
##' ## inspect results for a specific breakpoint
##' e4.lmbreak1010 <- as.lmbreak(e.mlmbreak1010, cluster = 4)
##' summary(e4.lmbreak1010)
##'
##' ## re-run with adpative step 
##' e4.lmbreak1010 <- lmbreak(Y ~ 0 + bp(X, "1010"), data = df1[df1$id==4,],
##'                           control = list(optimize.step = TRUE))
##' plot(e4.lmbreak1010)

## * mlmbreak (code)
##' @export
mlmbreak <- function(formula, data, cluster, trace = 1, cpus = 1, ...){

    ## ** normalize user input
    data <- as.data.frame(data)
    if(!is.character(cluster) || length(cluster)!=1){
        stop("Argument \'cluster\' should be character string of length 1. \n")
    }
    if(cluster %in% names(data) == FALSE){
        stop("Variable defined by argument \'cluster\' not found in argument \'data\'. \n")
    }
    if(cluster %in% all.vars(formula)){
        stop("Argument \'cluster\' cannot take value \"",cluster,"\" as it is a variable used in the mean model. \n")
    }
    if(cluster %in% c("breakpoint","estimate","se","df","lower","upper")){
        stop("Argument \'cluster\' cannot take value \"",cluster,"\" as it is a reserved name. \n")
    }

    ## ** run lmbreak on each cluster
    if(is.factor(data[[cluster]])){
        U.cluster <- levels(droplevels(data[[cluster]]))
    }else{
        U.cluster <- levels(as.factor(data[[cluster]]))
    }
    n.cluster <- length(U.cluster)
    trace2 <- ifelse(trace >= 2, trace-1, FALSE)
    
    if(cpus==1){

        if(trace>0 & trace <=1){
            ls.lmbreak <- pbapply::pblapply(U.cluster, function(iC){ ## iC <- U.cluster[1]
                lmbreak(formula, data = data[data[[cluster]]==iC,,drop=FALSE], trace = trace2, ...)
            })
        }else{
            ls.lmbreak <- lapply(U.cluster, function(iC){
                if(trace>1){
                    cat(cluster,": ",iC,"\n")
                }
                lmbreak(formula, data = data[data[[cluster]]==iC,,drop=FALSE], trace = trace2, ...)
            })
        }

    }else{

        cl <- parallel::makeCluster(cpus)
        if(trace>0){
            pb <- utils::txtProgressBar(max = n.cluster, style = 3)          
            progress <- function(n){utils::setTxtProgressBar(pb, n)}
            opts <- list(progress = progress)
        }else{
            opts <- list()
        }
        ## link to foreach
        doSNOW::registerDoSNOW(cl)

        ## export package
        parallel::clusterCall(cl, fun = function(x){
            suppressPackageStartupMessages(library(lmbreak, quietly = TRUE, warn.conflicts = FALSE, verbose = FALSE))
        })
        ## export functions
        toExport <- c("lmbreak.fit")
        iC <- NULL ## [:forCRANcheck:] foreach        
        ls.lmbreak <- foreach::`%dopar%`(
                                   foreach::foreach(iC=1:n.cluster,
                                                    .export = toExport,
                                                    .options.snow = opts),                                            
                                   {
                                       lmbreak(formula, data = data[data[[cluster]]==iC,,drop=FALSE], trace = FALSE, ...)
                                   })

        parallel::stopCluster(cl)
        if(trace>0){close(pb)}
    }
    
    ## ** export
    col.cluster <- as.data.frame(stats::setNames(list(U.cluster), cluster))
    out <- list(model = lapply(ls.lmbreak,"[[","model"),
                breakpoint = do.call(rbind,lapply(1:n.cluster, function(iC){cbind(U.cluster[iC],ls.lmbreak[[iC]]$breakpoint)})),
                opt = cbind(col.cluster,do.call(rbind,lapply(ls.lmbreak,"[[","opt"))),
                call = match.call(),
                args = c(ls.lmbreak[[1]]$args,list(cluster.var = cluster, U.cluster = U.cluster)),
                data = data)
    attr(out$breakpoint,"all") <- stats::setNames(lapply(1:n.cluster, function(iC){attr(ls.lmbreak[[iC]]$breakpoint,"all")}),U.cluster)
    attr(out$opt,"all") <- stats::setNames(lapply(1:n.cluster, function(iC){attr(ls.lmbreak[[iC]]$opt,"all")}),U.cluster)
    names(out$breakpoint)[1] <- cluster
    class(out) <- "mlmbreak"
    return(out)
}


##----------------------------------------------------------------------
### mlmbreak.R ends here
