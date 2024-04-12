### autoplot.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Apr  5 2024 (15:33) 
## Version: 
## Last-Updated: apr 12 2024 (13:07) 
##           By: Brice Ozenne
##     Update #: 145
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * autoplot.lmbreak (documentation)
#' @title Graphical Display of the Breakpoint Model
#' @description Display the fitted response value as a function of the breakpoint values
#' @name plot.lmbreak
#'
#' @param object,x output of \code{\link{lmbreak}}.
#' @param xlim [numeric vector of length 2] range of displayed values for the x-axis (breakpoint variable).
#' @param ylim [numeric vector of length 2] range of displayed values for the y-axis (response variable).
#' @param color [character vector of length 2] colors used to display the observations and the model fit.
#' @param size [numeric vector of length 2] size of the points representing the observations and size of the points representing the breakpoints.
#' @param linewidth [numeric] width of the line representing the model fit
#' @param alpha [numeric] transparency of the points representing the observations. Can be set to 0 to hide the observations and only display the fitted broken lines.
#' @param scales [character] argument passed to \code{ggplot2::facet_wrap} when displaying covariate specific model fit.
#' @param title [character] character string to be displayed at the top of the graphical window.
#' @param ... not used. For compatibility with the generic function.
#' 
#' @return A list containing a ggplot2 object (element \code{plot}) and the dataset with the fitted values (element \code{data})
#' 
#' @keywords hplot
#' @export 
autoplot.lmbreak <- function(object, xlim = NULL, ylim = NULL,
                             color = grDevices::palette.colors(2), size = c(1.5,2.5), linewidth = 1, alpha = 1, scales = "fixed", title = NULL, ...){

    ## ** extract from model
    breakpoint <- object$breakpoint$value
    breakpoint.var <- object$args$breakpoint.var
    data <- object$data
    response.var <- object$args$response.var
    Z.var <- object$args$covariate
    cv <- object$opt$cv
    continuity <- object$opt$continuity
    pattern <- object$opt$pattern

    ## ** check user input
    dots <- list(...)
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }

    ## xlim
    if(is.null(xlim)){
        xlim <- range(data[[breakpoint.var]], na.rm = TRUE)        
    }else{
        if(!is.numeric(xlim) || length(xlim)!=2){
            stop("Argument \'xlim\' should be a numeric vector of length 2. \n")
        }
        if(xlim[1]>=xlim[2]){
            stop("The first element of argument \'xlim\' should be strictly smaller than the second. \n")
        }
    }    

    ## ylim
    if(is.null(ylim) && ("ylim" %in% names(match.call()) == FALSE)){
        ylim <- range(data[[response.var]], na.rm = TRUE)
    }
    
    ## ** design matrix
    if(all(!is.na(breakpoint))){
        ls.X <- mapply(xxx = c(xlim[1], breakpoint),
                       yyy = c(breakpoint-1e-10, xlim[2]),
                       FUN = function(xxx, yyy){seq(xxx,yyy, length.out = 3)},
                       SIMPLIFY = FALSE)

        if(length(Z.var)>0){
            newdata <- expand.grid(c(list(unlist(ls.X)),unique(data[Z.var])))
            names(newdata)[1] <- breakpoint.var
        }else{
            newdata <- data.frame(unlist(ls.X))
            names(newdata) <- breakpoint.var
        }
    }

    ## ** fitted values
    if(all(!is.na(breakpoint))){
        newdataA <- stats::predict(object, newdata = newdata, keep.newdata = TRUE)
        newdataA$breakpoint <- newdataA[[breakpoint.var]] %in% breakpoint
        newdataB <- newdataA[newdataA$breakpoint,,drop=FALSE]    
    }else{
        newdataA <- NULL
    }

    ## ** graphical display
    out <- ggplot2::ggplot()
    out <- out + ggplot2::geom_point(data = data, ggplot2::aes(x = .data[[breakpoint.var]], y = .data[[response.var]], color = "observation"), alpha = alpha, size = size[1])
    if(!is.na(size[2]) && all(!is.na(breakpoint))){
        out <- out + ggplot2::geom_point(data = newdataB, ggplot2::aes(x = .data[[breakpoint.var]], y = .data$estimate, color = "fit"), size = size[2])
    }
    if(all(!is.na(breakpoint))){
        out <- out + ggplot2::geom_line(data = newdataA, ggplot2::aes(x = .data[[breakpoint.var]], y = .data$estimate, color = "fit"), linewidth = linewidth)
    }
    out <- out + ggplot2::scale_colour_manual(name = c("observation","fit"), values = stats::setNames(color,c("observation","fit")))
    if(length(Z.var)>0){
        out <- out + ggplot2::facet_wrap(stats::as.formula(paste0("~",paste0(Z.var,collapse="+"))), scales = scales)
    }
    if(!is.null(ylim)){
        out <- out + ggplot2::coord_cartesian(ylim = ylim)
    }else if(alpha<1e-10 && "ylim" %in% names(match.call()) == FALSE){
        out <- out + ggplot2::coord_cartesian(ylim = range(newdataA$estimate))
    }
    if(is.null(title)){
        out <- out + ggplot2::ggtitle(label = paste0("Pattern ",paste(pattern, collapse="")," (convergence: ",cv,", continuity: ",continuity,")"))
    }else if(any(!is.na(title))){
        out <- out + ggplot2::ggtitle(label = title)
    }

    ## ** export
    return(list(plot = out,
                data = newdataA))
}


## * autoplot.mlmbreak (documentation)
#' @title Graphical Display for Multiple Breakpoint Model
#' @description Display the fitted response value for each breapoint model as a function of the breakpoint values
#' @name plot.mlmbreak
#'
#' @param object,x output of \code{\link{mlmbreak}}.
#' @param cluster [vector] cluster relative to which the breakpoint model should be displayed.
#' @param xlim [numeric vector of length 2] range of displayed values for the x-axis (breakpoint variable).
#' @param ylim [numeric vector of length 2] range of displayed values for the y-axis (response variable).
#' @param color [character vector of length 2] colors used to display the observations and the model fit.
#' @param size [numeric vector of length 2] size of the points representing the observations and size of the points representing the breakpoints.
#' @param linewidth [numeric] width of the line representing the model fit
#' @param alpha [numeric] transparency of the points representing the observations. Can be set to 0 to hide the observations and only display the fitted broken lines.
#' @param scales [character] argument passed to \code{ggplot2::facet_wrap} when displaying covariate specific model fit.
#' @param labeller [character] function used to label each facet. Passed to \code{ggplot2::facet_wrap}.
#' @param subtitle [logical] should the subtitle include information about convergence of the model.
#' @param ... not used. For compatibility with the generic function.
#' 
#' @return A list containing a ggplot2 object (element \code{plot}) and the dataset with the fitted values (element \code{data})
#' 
#' @keywords hplot
#' @export 
autoplot.mlmbreak <- function(object, cluster = NULL, xlim = NULL, ylim = NULL,
                              color = grDevices::palette.colors(2), size = c(1.5,2.5), linewidth = 1, alpha = 1, scales = "fixed", labeller = "label_both",
                              subtitle = 2, ...){
    
    ## ** extract from object
    var.cluster <- object$args$cluster
    U.cluster <- object$args$U.cluster
    data <- object$data
    breakpoint.var <- object$args$breakpoint.var
    response.var <- object$args$response.var
    opt <- object$opt

    ## ** normalize user input
    dots <- list(...)
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }
    if(is.null(cluster)){
        cluster <- U.cluster
    }else if(any(cluster %in% U.cluster == FALSE)){
        stop("Unknown value for argument \'cluster\'.")
    }

    if(length(object$args$covariate)>0){
        stop("Cannot provide a graphical disply for mlmbreak object with covariate(s). \n")
    }

    ## ylim
    if(is.null(ylim) && ("ylim" %in% names(match.call()) == FALSE)){
        ylim <- range(data[[response.var]], na.rm = TRUE)
    }
    if(is.null(xlim) && ("xlim" %in% names(match.call()) == FALSE)){
        xlim <- range(data[[breakpoint.var]], na.rm = TRUE)
    }
    

    ## ** fitted values
    ls.ggdata <- lapply(cluster, function(iC){
        iOut <- autoplot(as.lmbreak(object, cluster = iC), xlim = xlim)$data
        if(!is.null(iOut)){
            return(cbind(iC,iOut))
        }else{
            return(NULL)
        }        
    })    
    iOut <- autoplot(as.lmbreak(object, cluster = cluster[1]), xlim = xlim)$data

    newdataA <- do.call(rbind,ls.ggdata)
    names(newdataA)[1] <- var.cluster

    if(subtitle>1){
        old2new <- paste0(U.cluster," (cv = ",opt$cv,",",opt$continuity,")")
        newdataA[[var.cluster]] <- factor(newdataA[[var.cluster]], levels = U.cluster, labels = old2new)
        data[[var.cluster]] <- factor(data[[var.cluster]], levels = U.cluster, labels = old2new)
    }else if(subtitle>0){
        old2new <- paste0(U.cluster," (cv=",as.numeric(opt$cv),")")
        newdataA[[var.cluster]] <- factor(newdataA[[var.cluster]], levels = U.cluster, labels = old2new)
        data[[var.cluster]] <- factor(data[[var.cluster]], levels = U.cluster, labels = old2new)
    }else{
        newdataA[[var.cluster]] <- factor(newdataA[[var.cluster]], levels = U.cluster)
    }
    newdataB <- newdataA[newdataA$breakpoint,,drop=FALSE]

    ## ** graphical display
    out <- ggplot2::ggplot()
    out <- out + ggplot2::geom_point(data = data, ggplot2::aes(x = .data[[breakpoint.var]], y = .data[[response.var]], color = "observation"), alpha = alpha, size = size[1])
    if(!is.na(size[2])){
        out <- out + ggplot2::geom_point(data = newdataB, ggplot2::aes(x = .data[[breakpoint.var]], y = .data$estimate, color = "fit"), size = size[2])
    }
    out <- out + ggplot2::geom_line(data = newdataA, ggplot2::aes(x = .data[[breakpoint.var]], y = .data$estimate, color = "fit"), linewidth = linewidth)
    out <- out + ggplot2::scale_colour_manual(name = c("observation","fit"), values = stats::setNames(color,c("observation","fit")))
    out <- out + ggplot2::facet_wrap(stats::as.formula(paste0("~",var.cluster)), scales = scales, labeller = labeller)

    if(!is.null(ylim)){
        out <- out + ggplot2::coord_cartesian(xlim = xlim, ylim = ylim)
    }else if(alpha<1e-10 && "ylim" %in% names(match.call()) == FALSE){
        out <- out + ggplot2::coord_cartesian(xlim = xlim, ylim = range(newdataA$estimate))
    }
    
    ## ** export
    return(list(plot = out,
                data = newdataA))


}

##----------------------------------------------------------------------
### autoplot.R ends here
