### autoplot.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Apr  5 2024 (15:33) 
## Version: 
## Last-Updated: Apr  8 2024 (11:47) 
##           By: Brice Ozenne
##     Update #: 54
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
#'
#' @param object output of \code{lmbreak}
#' @param xlim [numeric vector of length 2] range of displayed values for the x-axis (breakpoint variable).
#' @param ylim [numeric vector of length 2] range of displayed values for the y-axis (response variable).
#' @param color [character vector of length 2] colors used to display the observations and the model fit.
#' @param size [numeric vector of length 2] size of the points representing the observations and size of the points representing the breakpoints.
#' @param linewidth [numeric] width of the line representing the model fit
#' @param alpha [numeric] transparency of the points representing the observations. Can be set to 0 to hide the observations and only display the fitted broken lines.
#' @param title [character] character string to be displayed at the top of the graphical window.
#' @param ... not used. For compatibility with the generic function.
#' 
#' @return A list containing a ggplot2 object (element \code{plot}) and the dataset with the fitted values (element \code{data})
#' 

## * autoplot.lmbreak (code)
#' @method autoplot lmBreak
#' @export 
autoplot.lmbreak <- function(object, xlim = NULL, ylim = NULL,
                             color = grDevices::palette.colors(2), size = c(1.5,2.5), linewidth = 1, alpha = 1, title = NULL, ...){

    ## ** extract from model
    breakpoint <- object$breakpoint$value
    breakpoint.var <- object$breakpoint.var
    data <- object$data
    all.var <- all.vars(object$model$terms)
    response.var <- setdiff(all.var, all.vars(delete.response(object$model$terms)))
    
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

    ## ** design matrix
    ls.X <- mapply(xxx = c(xlim[1], breakpoint),
                   yyy = c(object$breakpoint$value-1e-10, xlim[2]),
                   FUN = function(xxx, yyy){seq(xxx,yyy, length.out = 3)},
                   SIMPLIFY = FALSE)

    Z.var <- setdiff(all.var, c(response.var,"Us0",paste0("Us",1:length(breakpoint)),paste0("Vs",1:length(breakpoint))))
    if(length(Z.var)>0){
        browser()
    }else{
        newdata <- data.frame(unlist(ls.X))
        names(newdata) <- breakpoint.var
    }
    
    ## ** fitted values
    newdataA <- predict(object, newdata = newdata, keep.newdata = TRUE)
    newdataB <- newdataA[newdataA[[breakpoint.var]] %in% breakpoint,,drop=FALSE]    

    ## ** graphical display
    out <- ggplot2::ggplot()
    out <- out + ggplot2::geom_point(data = data, aes(x = .data[[breakpoint.var]], y = .data[[response.var]], color = "observation"), alpha = alpha, size = size[1])
    if(!is.na(size[2])){
        out <- out + ggplot2::geom_point(data = newdataB, aes(x = .data[[breakpoint.var]], y = .data$estimate, color = "fit"), size = size[2])
    }
    out <- out + ggplot2::geom_line(data = newdataA, aes(x = .data[[breakpoint.var]], y = .data$estimate, color = "fit"), linewidth = linewidth)
    out <- out + ggplot2::scale_colour_manual(name = c("observation","fit"), values = stats::setNames(color,c("observation","fit")))

    if(!is.null(ylim)){
        out <- out + ggplot2::coord_cartesian(ylim = ylim)
    }else if(alpha<1e-10 && "ylim" %in% names(match.call()) == FALSE){
        out <- out + ggplot2::coord_cartesian(ylim = range(newdataA$estimate))
    }
    if(is.null(title)){
        out <- out + ggplot2::ggtitle(label = paste0("Pattern ",paste(object$pattern, collapse="")," (convergence: ",object$cv,", continuity: ",object$continuity,")"))
    }else if(any(!is.na(title))){
        out <- out + ggplot2::ggtitle(label = title)
    }

    ## ** export
    return(list(plot = out,
                data = newdataB))
}



##----------------------------------------------------------------------
### autoplot.R ends here
