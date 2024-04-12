### utils.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr 12 2024 (10:19) 
## Version: 
## Last-Updated: apr 12 2024 (10:34) 
##           By: Brice Ozenne
##     Update #: 7
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * combn2
## same as utils::combn but handle certain special case instead of returning an error
## called by lmbreak when initializing with gam or quantiles
combn2 <- function(x, m, name = NULL){
    if(is.null(x)||length(x)<m){
        out <- NULL
    }else{
        out <- utils::combn(x=x, m=m)
        if(!is.null(name)){
            if(length(name)==1){name <- rep(name,NCOL(out))}
            colnames(out) <- name
        }
    }

    return(out)
}

##----------------------------------------------------------------------
### utils.R ends here
