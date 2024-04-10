### print.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Apr  8 2024 (10:00) 
## Version: 
## Last-Updated: apr 10 2024 (15:35) 
##           By: Brice Ozenne
##     Update #: 36
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
    cat("\n\nConvergence: ", x$opt$cv, ", continuity: ",x$opt$continuity,"\n",sep="")

    return(invisible(TRUE))
}

## * print.mlmbreak (code)
#' @export 
print.mlmbreak <- function(x, digits = options()$digits-2, ...){
    
    df.res <- coef(x, type = c("pattern","cv","continuity"))
    df.res$breakpoint <- lapply(coef(x, type = c("breakpoint"), format = "list"),"[[","breakpoint")
    df.res$maxVs <- lapply(coef(x, type = c("Vs"), format = "list"),"[[","Vs")

    ## ** display
    cat("\nCall:\n")
    print(x$call)
    cat("\n")
    cat("Breakpoints:\n")
    df.res2print <- df.res
    df.res2print$maxVs <- lapply(df.res$maxVs,function(iDf){format.pval(max(abs(iDf)), eps = 10^-digits)})
    print(df.res2print, digits = digits)

    ## ** export
    return(invisible(df.res))
}
##----------------------------------------------------------------------
### print.R ends here
