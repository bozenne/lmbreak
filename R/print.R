### print.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Apr  8 2024 (10:00) 
## Version: 
## Last-Updated: apr  8 2024 (16:14) 
##           By: Brice Ozenne
##     Update #: 14
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * print.lmbreak (code)
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
