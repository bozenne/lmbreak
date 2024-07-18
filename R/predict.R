### predict.lmbreak.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Apr  8 2024 (08:49) 
## Version: 
## Last-Updated: jul 18 2024 (16:54) 
##           By: Brice Ozenne
##     Update #: 73
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * predict.lmbreak (documentation)
#' @title Fitted Values of the Breakpoint Model
#' @description Extract fitted values from the breakpoint model.
#' 
#' @param object [lmbreak] output of \code{\link{lmbreak}}.
#' @param newdata [data.frame] dataset containing the covariate to condition on when evaluated the expected outcome.
#' @param extrapolate [logical] should fitted values before the first or beyond the last observation be set to NA?
#' @param continuity [logical] should predictions be extracted from a breakpoint model ensuring continuity (i.e. no Vs terms)?
#' Often not relevant as Vs term should be 0 when proper convergence has been reached.
#' @param keep.newdata [logical] should the dataset be added as additional columns in the output?
#' @param ... Not used. For compatibility with the generic method.

## * predict.lmbreak
##' @export
predict.lmbreak <- function(object, newdata = NULL, extrapolate = FALSE, continuity = NULL, keep.newdata = TRUE, ...){
    
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
    breakpoint.var <- object$args$breakpoint.var

    ## ** predictions
    newdata.frame <- model.frame(object, newdata = newdata)
    out <- stats::predict(model, newdata = newdata.frame)
    
    if(extrapolate==FALSE){
        response.var <- object$args$response.var
        data <- object$data
        xlim.NNA <- range(data[!is.na(data[[response.var]]),breakpoint.var],na.rm=TRUE)
        out[newdata.frame[[breakpoint.var]] < xlim.NNA[1] | newdata.frame[[breakpoint.var]] > xlim.NNA[2]] <- NA
    }
    
    ## ** reformat
    if(keep.newdata){
        if(is.matrix(out)){
            out <- apply(out, MARGIN = 2, FUN = identity, simplify = FALSE)
        }
        if(is.numeric(out)){
            out <- cbind(newdata, estimate = out)
        }else if(is.list(out)){
            old2new <- c("fit" = "estimate","se.fit" = "se", "df" = "df", "lwr" = "lower", "upr" = "upper")
            out2 <- do.call(cbind, out)
            out2.red <- out2[,colnames(out2) %in% names(old2new),drop=FALSE] ## only keep specific columns
            colnames(out2.red) <- old2new[match(colnames(out2.red),names(old2new))] ## rename columns
            out <- cbind(newdata, out2.red[,intersect(old2new,colnames(out2.red)),drop=FALSE]) ## reorder columns and combine
        }else{
            stop("Unknown output format for prediction.lm(). \n")
        }
    }

    ## ** export
    return(out)
}

## * predict.mlmbreak (documentation)
#' @title Fitted Values of the Multiple Breakpoint Model
#' @description Extract fitted values from each breakpoint model.
#' 
#' @param object [mlmbreak] output of \code{\link{mlmbreak}}.
#' @param cluster [numeric or character] cluster relative to which the fit of the breakpoint model should be output.
#' @param ... additional arguments (newdata, extrapolate, continuity, keep.newdata) passed to \code{predict.lmbreak}.

## * predict.mlmbreak
##' @export
predict.mlmbreak <- function(object, cluster = NULL, ...){

    ## ** extract from object
    var.cluster <- object$args$cluster
    U.cluster <- object$args$U.cluster

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
    
    ## ** extract prediction
    ls.pred <- lapply(cluster.level, function(iC){
        iOut <- predict(as.lmbreak(object, cluster = iC), ...)
        iOutA <- cbind(stats::setNames(list(rep(iC,NROW(iOut))),var.cluster),iOut)
        return(iOutA)
    })
    out <- do.call(rbind,ls.pred)

    ## ** export
    rownames(out) <- NULL
    return(out)
}
##----------------------------------------------------------------------
### predict.lmbreak.R ends here
