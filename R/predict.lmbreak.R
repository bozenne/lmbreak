### predict.lmbreak.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Apr  8 2024 (08:49) 
## Version: 
## Last-Updated: Apr  8 2024 (11:53) 
##           By: Brice Ozenne
##     Update #: 32
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * predict.lmbreak
##' @export
predict.lmbreak <- function(object, newdata = NULL, keep.newdata = TRUE, ...){
    
    ## ** predictions
    if(object$continuity == FALSE && !is.null(attr(object$model,"continuity"))){
        out <- predict(attr(object$model,"continuity"), model.frame(object, newdata = newdata), ...)
    }else{
        out <- predict(object$model, model.frame(object, newdata = newdata), ...)
    }
        
    ## ** reformat
    if(keep.newdata){
        if(is.matrix(out)){
            out <- apply(out, MARGIN = 2, FUN = identity, simplify = FALSE)
        }
        if(is.numeric(out)){
            out <- cbind(newdata, estimate = out)
        }else if(is.list(out)){
            old2new <- c("fit" = "estimate","se.fit" = "se", "df" = "df", "lwr" = "lower", "upr" = "upper")
            out2 <- do.call(cbind, out)
            out2.red <- out2[,colnames(out2) %in% names(old2new),drop=FALSE] ## only keep specific columns
            colnames(out2.red) <- old2new[match(colnames(out2.red),names(old2new))] ## rename columns
            out <- cbind(newdata, out2.red[,intersect(old2new,colnames(out2.red)),drop=FALSE]) ## reorder columns and combine
        }else{
            stop("Unknown output format for prediction.lm(). \n")
        }
    }

    ## ** export
    return(out)
}
##----------------------------------------------------------------------
### predict.lmbreak.R ends here
