### coef.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Apr  6 2024 (12:23) 
## Version: 
## Last-Updated: apr 10 2024 (16:28) 
##           By: Brice Ozenne
##     Update #: 89
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * coef.lmbreak (documentation)
##' @title Extract Summary Statistics From Breakpoint Model
##' @description Extract summary statistics from the breakpoint model
##'
##' @param object output of \code{\link{lmbreak}}.
##' @param type [vector of character] summary statistic to be output: \itemize{
##' \item breakpoint: position of the breakpoint
##' \item breakpoint.range: position of the breakpoint, minimum and maximum X values.
##' \item duration: duration (in term of X values) of each broken line.
##' \item Us: slope terms in the linear regression model used to estimates the parameters.
##' \item slope: slope associated to each breakpoint.
##' \item Vs: discontinuity terms in the linear regression model used to estimates the parameters.
##' \item intercept: expected outcome value at each breakpoint value.
##' \item lm: regression coefficients of the linear regression model used to estimates the parameters.
##' }
##' @param simplify [logical] simplify the data format from a list to a vector or a data.frame.
##' @param ... Not used. For compatibility with the generic method.
##'
##' @return A vector when a single type is requested. A list or a data.frame otherwise.
##'  
##' @keywords methods
##' @export
coef.lmbreak <- function(object, type = "breakpoint", simplify = TRUE, ...){

    ## ** normalize user input
    dots <- list(...)
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }

    valid.types <- c("breakpoint","breakpoint.range","duration",
                     "Us","slope",
                     "Vs","intercept",
                     "lm",
                     "pattern","cv","continuity")
    if(any(type %in% valid.types == FALSE)){
        stop("Incorrect value for argument \'type\'. \n",
             "Invalid values: \"",paste(setdiff(type,valid.types), collapse = "\", \""),"\". \n",
             "Valid values: \"",paste(setdiff(valid.types,type), collapse = "\", \""),"\". \n")
    }

    ## ** extract coefficient
    data <- object$data
    breakpoint.var <- object$args$breakpoint.var
    response.var <- object$args$response.var
    table.breakpoint <- object$breakpoint
    continuity <- object$opt$continuity
    model <- object$model
    object.opt <- object$opt
    vec.pattern <- object$args$pattern[[object.opt$pattern]]
    out <- list()

    if(any(c("pattern") %in% type)){
        out$pattern <- object.opt$pattern
    }
    if(any(c("cv") %in% type)){
        out$cv <- object.opt$cv
    }
    if(any(c("continuity") %in% type)){
        out$continuity <- object.opt$continuity
    }
    if(any(c("bp","breakpoint","breakpoint-range","duration","intercept") %in% type)){
        out$breakpoint <- table.breakpoint$value
    }
    if(any(c("breakpoint.range","duration","intercept") %in% type)){
        out$breakpoint.range <- c(min(data[[breakpoint.var]][!is.na(data[[response.var]])], na.rm = TRUE), out$breakpoint, max(data[[breakpoint.var]][!is.na(data[[response.var]])], na.rm = TRUE))
    }
    if(any(c("duration","intercept") %in% type)){
        out$duration <- diff(out$breakpoint.range)
    }
    if(any(c("Us","slope","intercept") %in% type)){
        name.Us <- stats::na.omit(c(ifelse("Us0" %in% names(coef(model)),"Us0",NA),table.breakpoint$Us))
        if(continuity == FALSE && !is.null(attr(model,"continuity"))){
            out$Us <- coef(attr(model,"continuity"))[name.Us]
        }else{
            out$Us <- coef(model)[name.Us]
        }
    }
    if(any(c("slope","intercept") %in% type)){
        Us2slope <- stats::setNames(table.breakpoint$sign,table.breakpoint$Us)
        if(vec.pattern[1]=="1"){
            Us2slope <- c(stats::setNames(1,names(out$Us)[1]),Us2slope)
        }
        out$slope <- unname(cumsum(Us2slope * out$Us[names(Us2slope)]))
    }
    if(any(c("Vs") %in% type)){
        out$Vs <- coef(model)[table.breakpoint$Vs]
    }
    if(any(c("intercept") %in% type)){
        out$intercept <- cumsum(c(0,out$slope*out$duration))
        if(attr(model$terms,"intercept")){
            out$intercept <-  out$intercept + coef(model)["(Intercept)"]
        }
    }
    if(any(c("lm") %in% type)){
        out$beta <- coef(model)
    }

    ## ** output
    if(simplify){
        if(length(type) == 1){
            return(out[[type]])
        }else if(length(unique(lengths(out[type])))==1){
            return(as.data.frame(out))
        }else{
            return(out[type])
        }
    }else{
        return(out[type])
    }
}

## * coef.lmbreak (documentation)
##' @title Extract Summary Statistics From Multiple Breakpoint Model
##' @description Extract summary statistics from each breakpoint model
##'
##' @param object output of \code{\link{mlmbreak}}.
##' @param cluster [vector] cluster relative to which the summary statistics should be extracted..
##' @param format [character] should the output be a data.frame (with cluster as a column) or a list (with one element for each cluster)?
##' @param ... additional arguments passed to \code{\link{coef.lmbreak}}.
##'
##' @return A data.frame or a list.
##'  
##' @keywords methods
##' @export
coef.mlmbreak <- function(object, cluster = NULL, format = "data.frame", ...){

    ## ** extract from object
    var.cluster <- object$args$cluster
    U.cluster <- object$args$U.cluster

    ## ** normalize user input
    if(is.null(cluster)){
        cluster <- U.cluster
    }else if(any(cluster %in% U.cluster == FALSE)){
        stop("Unknown value for argument \'cluster\'.")
    }
    format <- match.arg(format, c("data.frame","list"))


    ## ** extract
    ls.table <- lapply(cluster, function(iC){ ## iC <- cluster[1]
        iOut <- coef(as.lmbreak(object, cluster = iC), simplify = FALSE,...)
        if(length(unique(lengths(iOut)))!=1){
            iMax.size <- max(lengths(iOut))
            iOut <- lapply(iOut, function(ii){c(ii,rep(NA, iMax.size-length(ii)))})
        }
        iOut <- as.data.frame(c(iC,iOut))
        names(iOut)[1] <- var.cluster
        return(iOut)
    })    
    ## ** export
    if(format=="list"){
        return(ls.table)
    }else if(format=="data.frame"){
        return(do.call(rbind,ls.table))
    }
}
##----------------------------------------------------------------------
### coef.R ends here
