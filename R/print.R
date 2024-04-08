### print.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Apr  8 2024 (10:00) 
## Version: 
## Last-Updated: Apr  8 2024 (10:07) 
##           By: Brice Ozenne
##     Update #: 13
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * print.lmbreak (code)
#' @method print lmBreak
#' @export 
print.lmbreak <- function(x, digits = options()$digits, ...){

    object <- x$model
    object$call <- x$call
    print(object)
    cat("Breakpoints: \n", paste0(round(x$breakpoint$value, digits = digits), collapse = ", "))
    cat("\n\nConvergence: ", x$cv, ", continuity: ",x$continuity,"\n",sep="")

    return(invisible(TRUE))
}

##----------------------------------------------------------------------
### print.R ends here
