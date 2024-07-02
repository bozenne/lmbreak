### coef.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Apr  6 2024 (12:23) 
## Version: 
## Last-Updated: jul  2 2024 (15:33) 
##           By: Brice Ozenne
##     Update #: 125
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
##' \item auc: area under the curve
##' \item breakpoint: position of the breakpoint
##' \item breakpoint.range: position of the breakpoint, minimum and maximum X values.
##' \item duration: duration (in term of X values) of each broken line.
##' \item Us: slope terms in the linear regression model used to estimates the parameters.
##' \item slope: slope associated to each breakpoint.
##' \item Vs: discontinuity terms in the linear regression model used to estimates the parameters.
##' \item intercept: expected outcome value at each breakpoint value.
##' \item lm: regression coefficients of the linear regression model used to estimates the parameters.
##' \item R2: coefficient of determination of the model fit.
##' \item pattern: pattern used.
##' \item cv: convergence of the optimization algorithm.
##' \item continuity: continuity at the breakpoints of the model fit.
##' }
##' @param interval [numeric vector of length 2] values from where to start and to stop calculating the area under the curve.
##' Only relevant when \code{type="auc"}.
##' @param simplify [logical] simplify the data format from a list to a vector or a data.frame.
##' @param ... Not used. For compatibility with the generic method.
##'
##' @return A vector when a single type is requested. A list or a data.frame otherwise.
##'  
##' @keywords methods
##' @export
coef.lmbreak <- function(object, type = "breakpoint", interval, simplify = TRUE, ...){


    ## ** extract from object
    data <- object$data
    breakpoint.var <- object$args$breakpoint.var
    response.var <- object$args$response.var
    table.breakpoint <- object$breakpoint
    continuity <- object$opt$continuity
    model <- object$model
    object.opt <- object$opt
    vec.pattern <- object$args$pattern[[object.opt$pattern]]
    
    ## ** normalize user input
    dots <- list(...)
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }

    valid.types <- c("auc",
                     "breakpoint","breakpoint.range","duration",
                     "Us","slope",
                     "Vs","intercept",
                     "lm", "R2",
                     "pattern","cv","continuity")
    names(valid.types) <- tolower(valid.types)
    type <- tolower(type)
    if(any(type %in% names(valid.types) == FALSE)){
        stop("Incorrect value for argument \'type\'. \n",
             "Invalid values: \"",paste(setdiff(type,valid.types), collapse = "\", \""),"\". \n",
             "Valid values: \"",paste(setdiff(valid.types,type), collapse = "\", \""),"\". \n")
    }

    if("auc" %in% type){
        if(missing(interval)){
            stop("Argument \'interval\' should be specifying when argument \'type\' contains \"auc\". \n")
        }
        if(length(interval)!=2 & !is.numeric(interval)){
            stop("Argument \'interval\' should be a numeric vector of length 2. \n")
        }

        range.obs <- range(data[[breakpoint.var]], na.rm=TRUE)
        if(min(interval)<range.obs[1] || max(interval)>range.obs[2]){
            stop("Argument \'interval\' should take value between ",range.obs[1]," and ",range.obs[2],". \n")
        }
    }

    ## ** extract coefficient
    out <- list()

    if(any(c("auc") %in% type)){
        ggdata <- autoplot(object)$data
        out$auc <- AUC(x = ggdata[[breakpoint.var]], y = ggdata$estimate, from = interval[1], to = interval[2],
                       method = "trapezoid", na.rm = TRUE)
    }
    if(any(c("pattern") %in% type)){
        out$pattern <- object.opt$pattern
    }
    if(any(c("r2") %in% type)){
        out$R2 <- object.opt$r2
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
    if(any(c("us","slope","intercept") %in% type)){
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
            out$slope <- unname(cumsum(Us2slope * out$Us[names(Us2slope)]))
        }else{
            out$slope <- c(0,cumsum(Us2slope * out$Us[names(Us2slope)]))
        }
    }
    if(any(c("vs") %in% type)){
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
    ## use valid.types[type] to handle upper/lower case
    if(simplify){
        if(length(type) == 1){
            return(out[[valid.types[type]]])
        }else if(length(unique(lengths(out[valid.types[type]])))==1){
            return(as.data.frame(out))
        }else{
            return(out[valid.types[type]])
        }
    }else{
        return(out[valid.types[type]])
    }
}

## * coef.mlmbreak (documentation)
##' @title Extract Summary Statistics From Multiple Breakpoint Model
##' @description Extract summary statistics from each breakpoint model
##'
##' @param object output of \code{\link{mlmbreak}}.
##' @param type [character vector] summary statistic to be output. See argument \code{type} in \code{\link{coef.lmbreak}}.
##' @param cluster [vector] cluster relative to which the summary statistics should be extracted..
##' @param format [character] should the output be a data.frame (with cluster as a column) or a list (with one element for each cluster)?
##' @param ... additional arguments passed to \code{\link{coef.lmbreak}}.
##'
##' @return A data.frame or a list.
##'  
##' @keywords methods
##' @export
coef.mlmbreak <- function(object, type = "breakpoint", cluster = NULL, format = "data.frame", ...){

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
        iOut <- coef(as.lmbreak(object, cluster = iC), type = type, simplify = FALSE,...)
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
