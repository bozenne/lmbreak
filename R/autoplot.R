### autoplot.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Apr  5 2024 (15:33) 
## Version: 
## Last-Updated: jun 10 2024 (16:12) 
##           By: Brice Ozenne
##     Update #: 209
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
#' @param object,x [lmbreak] output of \code{\link{lmbreak}}.
#' @param y [lmbreak] another model whose fit is to be compared with the first one.
#' @param xlim [numeric vector of length 2] range of displayed values for the x-axis (breakpoint variable).
#' @param ylim [numeric vector of length 2] range of displayed values for the y-axis (response variable).
#' @param color [character vector of length 2 or 3] colors used to display the observations and the model fit.
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
autoplot.lmbreak <- function(object, y = NULL, xlim = NULL, ylim = NULL,
                             color = grDevices::palette.colors(3), size = c(1.5,2.5), linewidth = 1, alpha = 1, scales = "fixed", title = NULL, ...){

    ## ** extract from model
    breakpoint <- stats::coef(object, type = "breakpoint")
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

    ## names
    if(any(c(object$args$response.var,object$args$breakpoint.var) %in% "model")){
        stop("The response and breakpoint variable should not be called \"model\". \n",
             "This name is used internally for buidling the graphical display. \n")
    }
    if(any(names(data) %in% "model")){
        stop("The column of argument \'data\' (used when calling lmbreak or mlmbreak) should not be called \"model\". \n",
             "This name is used internally for buidling the graphical display. \n")
    }

    ## y
    if(!is.null(y)){
        if(!inherits(y,"lmbreak")){
            stop("Incorrect argument \'y\': should inherit of lmbreak. \n")
        }
        attr(breakpoint,"y") <- stats::coef(y,"breakpoint")
        allBreakpoint <- unique(sort(c(breakpoint,attr(breakpoint,"y"))))
    }else{
        allBreakpoint <- breakpoint
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
        ls.X <- mapply(xxx = c(xlim[1], allBreakpoint),
                       yyy = c(allBreakpoint-1e-10, xlim[2]),
                       FUN = function(xxx, yyy){seq(xxx,yyy, length.out = 25)},
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
        newdataA$model <- "fit"
        if(!is.null(y)){
            if(y$opt$optimizer[1]!=object$opt$optimizer[1]){
                newdataA$model <- paste0(newdataA$model," (",object$opt$optimizer[1],")")
            }else{
                newdataA$model <- paste0(newdataA$model," (1)")
            }
        }
        newdataB <- newdataA[newdataA$breakpoint,,drop=FALSE]    
    }else{
        newdataA <- NULL
    }

    if(!is.null(y) && all(!is.na(attr(breakpoint,"y")))){
        newdataA.bis <- stats::predict(y, newdata = newdata, keep.newdata = TRUE)
        newdataA.bis$breakpoint <- newdataA.bis[[breakpoint.var]] %in% attr(breakpoint,"y")
        if(y$opt$optimizer[1]!=object$opt$optimizer[1]){
            newdataA.bis$model <- paste0("fit (",y$opt$optimizer[1],")")
        }else{
            newdataA.bis$model <- paste0("fit (2)")
        }

        newdataB.bis <- newdataA.bis[newdataA.bis$breakpoint,,drop=FALSE]

        newdataA <- rbind(newdataA, newdataA.bis)
        newdataB <- rbind(newdataB, newdataB.bis)
    }
    
    
    ## ** graphical display
    out <- ggplot2::ggplot()
    out <- out + ggplot2::geom_point(data = cbind(data, model = "observation"), ggplot2::aes(x = .data[[breakpoint.var]], y = .data[[response.var]], color = .data$model), alpha = alpha, size = size[1])
    if(!is.na(size[2]) && all(!is.na(breakpoint))){
        out <- out + ggplot2::geom_point(data = newdataB, ggplot2::aes(x = .data[[breakpoint.var]], y = .data$estimate, color = .data$model), size = size[2])
    }
    if(all(!is.na(breakpoint))){
        out <- out + ggplot2::geom_line(data = newdataA, ggplot2::aes(x = .data[[breakpoint.var]], y = .data$estimate, group = .data$model, color = .data$model), linewidth = linewidth)
    }
    out <- out + ggplot2::scale_colour_manual(values = stats::setNames(color[1:(2+!is.null(y))], c("observation",unique(newdataA$model))))
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
#' @param object,x [mlmbreak] output of \code{\link{mlmbreak}}.
#' @param y [mlmbreak] another model whose fit is to be compared with the first one.
#' @param cluster [vector] cluster relative to which the breakpoint model should be displayed.
#' @param xlim [numeric vector of length 2] range of displayed values for the x-axis (breakpoint variable).
#' @param ylim [numeric vector of length 2] range of displayed values for the y-axis (response variable).
#' @param color [character vector of length 2 or 3] colors used to display the observations and the model fit.
#' @param size [numeric vector of length 2] size of the points representing the observations and size of the points representing the breakpoints.
#' @param linewidth [numeric] width of the line representing the model fit
#' @param alpha [numeric] transparency of the points representing the observations. Can be set to 0 to hide the observations and only display the fitted broken lines.
#' @param scales [character] argument passed to \code{ggplot2::facet_wrap} when displaying covariate specific model fit.
#' Can be \code{"none"} to display all trajectories on the same plot.
#' @param labeller [character] function used to label each facet. Passed to \code{ggplot2::facet_wrap}.
#' @param subtitle [logical] should the subtitle include information about convergence of the model.
#' @param ... not used. For compatibility with the generic function.
#' 
#' @return A list containing a ggplot2 object (element \code{plot}) and the dataset with the fitted values (element \code{data})
#' 
#' @keywords hplot
#' @export 
autoplot.mlmbreak <- function(object, cluster = NULL, y = NULL, xlim = NULL, ylim = NULL,
                              color = grDevices::palette.colors(3), size = c(1.5,2.5), linewidth = 1, alpha = 1, scales = "fixed", labeller = "label_both",
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

    ## y
    if(!is.null(y) && !inherits(y,"mlmbreak")){
        stop("Incorrect argument \'y\': should inherit of mlmbreak. \n")
    }

    ## scales
    if(is.null(scales) || any(is.na(scales))){
        scales <- "none"
    }
    scales <- match.arg(scales, c("none","fixed","free","free_x","free_y"))
    
    ## ylim
    if(is.null(ylim) && ("ylim" %in% names(match.call()) == FALSE)){
        if(scales == "fixed"){
            ylim <- range(data[[response.var]], na.rm = TRUE)
        }else{
            ylim <- NULL
        }
    }
    if(is.null(xlim) && ("xlim" %in% names(match.call()) == FALSE)){
        if(scales == "fixed"){
            xlim <- range(data[[breakpoint.var]], na.rm = TRUE)
        }else{
            xlim <- NULL
        }
    }
    

    ## ** fitted values
    ls.ggdata <- lapply(cluster, function(iC){
        if(!is.null(y)){
            iOut <- autoplot(as.lmbreak(object, cluster = iC), y = as.lmbreak(y, cluster = iC), xlim = xlim)$data
        }else{
            iOut <- autoplot(as.lmbreak(object, cluster = iC), xlim = xlim)$data
        }
        if(!is.null(iOut)){
            return(cbind(iC,iOut))
        }else{
            return(NULL)
        }        
    })    
    newdataA <- do.call(rbind,ls.ggdata)
    names(newdataA)[1] <- var.cluster

    if(is.null(y) && subtitle>1){
        old2new <- paste0(U.cluster," (cv = ",opt$cv,",",opt$continuity,")")
        newdataA[[var.cluster]] <- factor(newdataA[[var.cluster]], levels = U.cluster, labels = old2new)
        data[[var.cluster]] <- factor(data[[var.cluster]], levels = U.cluster, labels = old2new)
    }else if(is.null(y) && subtitle>0){
        old2new <- paste0(U.cluster," (cv=",as.numeric(opt$cv),")")
        newdataA[[var.cluster]] <- factor(newdataA[[var.cluster]], levels = U.cluster, labels = old2new)
        data[[var.cluster]] <- factor(data[[var.cluster]], levels = U.cluster, labels = old2new)
    }else{
        newdataA[[var.cluster]] <- factor(newdataA[[var.cluster]], levels = U.cluster)
    }
    newdataB <- newdataA[newdataA$breakpoint,,drop=FALSE]

    ## ** graphical display
    out <- ggplot2::ggplot()
    out <- out + ggplot2::geom_point(data = cbind(data, model = "observation"), ggplot2::aes(x = .data[[breakpoint.var]], y = .data[[response.var]], color = .data$model),
                                     alpha = alpha, size = size[1])
    if(!is.na(size[2])){
        out <- out + ggplot2::geom_point(data = newdataB, ggplot2::aes(x = .data[[breakpoint.var]], y = .data$estimate, color = .data$model), size = size[2])
    }
    if(scales == "none"){
        out <- out + ggplot2::geom_line(data = newdataA, ggplot2::aes(x = .data[[breakpoint.var]], y = .data$estimate, color = .data$model, group = .data[[var.cluster]]), linewidth = linewidth)
        out <- out + ggplot2::scale_colour_manual(values = stats::setNames(color[1:(2+!is.null(y))], c("observation",unique(newdataB$model))))
    }else{
        out <- out + ggplot2::geom_line(data = newdataA, ggplot2::aes(x = .data[[breakpoint.var]], y = .data$estimate, color = .data$model), linewidth = linewidth)
        out <- out + ggplot2::scale_colour_manual(values = stats::setNames(color[1:(2+!is.null(y))], c("observation",unique(newdataB$model))))
        out <- out + ggplot2::facet_wrap(stats::as.formula(paste0("~",var.cluster)), scales = scales, labeller = labeller)
    }

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
