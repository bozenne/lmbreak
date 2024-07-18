### print.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Apr  8 2024 (10:00) 
## Version: 
## Last-Updated: jul 18 2024 (11:32) 
##           By: Brice Ozenne
##     Update #: 55
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
print.lmbreak <- function(x, digits = options()$digits, continuity = NULL, ...){

    ## ** normalize user input
    ## *** dots
    dots <- list(...)
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }

    ## *** continuity
    if(is.null(continuity)){
        continuity <- (x$opt$continuity==FALSE)
    }

    ## ** extract from object
    if(continuity && !is.null(attr(x$model,"continuity"))){
        model <- attr(x$model,"continuity")
    }else{
        model <- x$model
    }

    ## ** print
    print(model)
    cat("Breakpoints: \n", paste0(round(x$breakpoint$value, digits = digits), collapse = ", "))
    cat("\n\nConvergence: ", x$opt$cv, ", continuity: ",x$opt$continuity,", regularity: ",x$opt$regularity,"\n",sep="")

    ## ** export
    return(invisible(TRUE))
}

## * print.mlmbreak (code)
#' @export 
print.mlmbreak <- function(x, digits = options()$digits-2, ...){

    df.res <- coef(x, type = c("pattern","cv","continuity"))
    ls.pattern <- tapply(df.res[[x$args$cluster]],df.res$pattern,function(iID){ ## iID <- df.res$PatientID[2:3]
        iChar <- paste0(ifelse(df.res$cv[match(iID,df.res[[x$args$cluster]])],"","*"),ifelse(df.res$continuity[match(iID,df.res[[x$args$cluster]])],"","|"))
        paste0(iID,ifelse(nchar(iChar)==0,"",paste0("(",iChar,")")))
    })
    
    ## ** display
    cat("\nCall:\n")
    print(x$call)
    cat("\n")
    if(all(df.res$cv==TRUE) & all(df.res$continuity==TRUE)){
        cat("Pattern:\n")
    }else if(all(df.res$cv==TRUE)){
        cat("Pattern (| no continuity):\n")
    }else if(all(df.res$continuity==TRUE)){
        cat("Pattern (* no cv):\n")
    }else{
        cat("Pattern (* no cv, | no continuity):\n")
    }
    print(ls.pattern, quote = FALSE, row.names = FALSE)
    ## df.res2print <- df.res
    ## df.res2print$maxVs <- lapply(df.res$maxVs,function(iDf){format.pval(max(abs(iDf)), eps = 10^-digits)})
    ## print(df.res2print, digits = digits)

    ## ** export
    return(invisible(TRUE))
}
##----------------------------------------------------------------------
### print.R ends here
