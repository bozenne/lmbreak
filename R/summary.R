### summary.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr  8 2024 (17:58) 
## Version: 
## Last-Updated: apr  8 2024 (18:37) 
##           By: Brice Ozenne
##     Update #: 28
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * summary.lmbreak (code)
#' @export 
summary.lmbreak <- function(x, digits = c(options()$digits,1), ...){

    cat("\t\tLinear regression with estimated breakpoints \n")
    
    object <- x$model
    object$call <- x$call
    print(summary(object))
    
    cat("Optmisation: convergence ", x$cv, ", continuity: ",x$continuity,"\n",sep="")
    if(!is.null(attr(x$breakpoint.init,"all"))){
        n.init <- NCOL(attr(x$breakpoint,"all"))
        M.diff <- attr(x$breakpoint,"all")- matrix(x$breakpoint$value, nrow = NROW(x$breakpoint), ncol = n.init, byrow = FALSE)
        cat("             ",sum(attr(x$cv,"all"))," (",round(100*sum(attr(x$cv,"all"))/n.init,digits[2]),"%) sucessfull convergence \n", sep="")
        if(x$cv){
            cat("             with difference in breakpoints between ",round(min(M.diff[,attr(x$cv,"all"),drop=FALSE]),digits[1]),
                " and ",round(max(M.diff[,attr(x$cv,"all"),drop=FALSE]),digits[1]),"\n",sep="")
        }
    }

    cat("\n")
    cat("Breakpoints: \n")
    df.breakpoint <- data.frame(estimate = x$breakpoint$value,
                                se = x$breakpoint$se,
                                lower = x$breakpoint$value + qnorm(0.025)*x$breakpoint$se,
                                upper = x$breakpoint$value + qnorm(0.975)*x$breakpoint$se)
    print(df.breakpoint, digits = digits[1])

    return(invisible(TRUE))
}


##----------------------------------------------------------------------
### summary.R ends here
