### summary.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr  8 2024 (17:58) 
## Version: 
## Last-Updated: jul 18 2024 (13:51) 
##           By: Brice Ozenne
##     Update #: 132
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
summary.lmbreak <- function(object, digits = c(options()$digits,1), continuity = NULL, ...){

    ## ** normalize user input
    if(is.null(continuity)){
        continuity <- (object$opt$continuity==FALSE)
    }
    
    ## ** extract from object
    object.pattern <- object$opt$pattern
    object.vecpattern <- object$arg$pattern[[object.pattern]]
    object.breakpoint <- object$breakpoint
    n.breakpoint <- length(object.breakpoint$value)
    objectAll.breakpoint <- attr(object$breakpoint,"all")
    nLS.init <- sapply(objectAll.breakpoint,NROW)
    n.init <- sum(nLS.init)
        
    object.opt <- object$opt 
    objectAll.opt <- attr(object.opt,"all")

    if(continuity && !is.null(attr(object$model,"continuity"))){
        model <- attr(object$model,"continuity")
    }else{
        model <- object$model
    }
    
    ## ** print
    cat("\t\tLinear regression with estimated breakpoints \n")

    model$call <- object$call
    print(summary(model))
    
    cat("Optmisation: convergence ", object.opt$cv, ", continuity: ",object.opt$continuity, ", regularity: ",object.opt$regularity,"\n",sep="")
    if(!is.null(objectAll.opt)){        

        index.pattern <- which(sapply(objectAll.breakpoint,NCOL)==n.breakpoint)
        index.pattern2 <- (c(0,cumsum(nLS.init))[index.pattern]+1):cumsum(nLS.init)[index.pattern]

        M.diff <- objectAll.breakpoint[[index.pattern]]- matrix(object.breakpoint$value, nrow = NROW(objectAll.breakpoint[[index.pattern]]), ncol = n.breakpoint, byrow = TRUE)
        cat("             ",sum(objectAll.opt$cv)," (",round(100*sum(objectAll.opt$cv)/n.init,digits[2]),"%) sucessfull convergence", sep="")

        if(object.opt$cv & NROW(objectAll.opt)>1){
            if(all(abs(M.diff[objectAll.opt$cv[index.pattern2],,drop=FALSE])<object.opt$tol[[1]][1])){
                cat(", all to the same breakpoint",ifelse(n.breakpoint>1,"s",""),".\n", sep = "")
            }else{
                cat("\n",
                    "             with difference in breakpoints between ",round(min(M.diff[objectAll.opt$cv[index.pattern2],,drop=FALSE]),digits[1]),
                    " and ",round(max(M.diff[objectAll.opt$cv[index.pattern2],,drop=FALSE]),digits[1]),"\n",sep="")
            }
        }else{
            cat("\n")
        }

    }

    cat("\n")
    txt.breakpoint <- paste0(n.breakpoint," breakpoint",ifelse(n.breakpoint>1,"s",""))
    txt.plateau <- ifelse(all(object.vecpattern==1),"",paste0(", ",sum(object.vecpattern==0)," plateau",ifelse(sum(object.vecpattern==0)>1,"x","")))
    cat("Pattern: ",paste(object.pattern,collapse="")," (",txt.breakpoint,txt.plateau,")\n",sep="")

    object.table <- model.tables(object, continuity = continuity)
    table2print <- format.data.frame(object.table, digits = digits[1], na.encode = FALSE)
    table2print[is.na(object.table)] <- ""
    print(table2print, row.names = FALSE)

    return(invisible(TRUE))
}

## * summary.mlmbreak (code)
#' @export 
summary.mlmbreak <- function(object, digits = c(options()$digits,1), ...){

    df.res <- coef(object, type = c("pattern","cv","continuity","r2"))
    df.res$breakpoint <- lapply(coef(object, type = c("breakpoint"), format = "list"),"[[","breakpoint")
    df.res$maxVs <- lapply(coef(object, type = c("Vs"), format = "list"),"[[","Vs")

    ## ** display
    cat("\nCall:\n")
    print(object$call)
    cat("\n")
    cat("Breakpoints:\n")
    df.res2print <- df.res
    df.res2print$maxVs <- lapply(df.res$maxVs,function(iDf){format.pval(max(abs(iDf)), eps = 10^-digits[1])})
    print(df.res2print, digits = digits, row.names = FALSE)

    ## ** export
    return(invisible(df.res))
}
##----------------------------------------------------------------------
### summary.R ends here
