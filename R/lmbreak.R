### lmbreak.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Apr  5 2024 (15:33) 
## Version: 
## Last-Updated: apr 10 2024 (16:36) 
##           By: Brice Ozenne
##     Update #: 546
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * lmbreak (documentation)
#' @title Fit Breakpoint model
#' @description Fit a linear regression with breakpoints.
#' @name lmBreak
#'
#' @param formula a formula where the breakpoint variable appears on the right hand side
#' with an argument pattern specifying the number of breakpoints and possible constrains. See details section.
#' @param data [data.frame] dataset
#' @param n.iter [integer, >0] the maximum number of iterations used to estimates the breakpoints.
#' @param tol [numeric, >0] the maximum accpetable difference between two consecutive estimates of the breakpoints.
#' When reached, the estimation algorithm stops.
#' @param enforce.continuity [logical] in the case where no continuous solution could be found,
#' should a non-continuous breakpoint model be kept (\code{FALSE}) or a continuous breakpoint model be enforced (\code{TRUE}, refiting without the Vs terms).
#' @param init.gam [logical] should a spline model be used to define initialization points (where the first derivative of the spline changes sign).
#' @param init.quantile [logical or integer] number of quantiles used to try different initializations.
#' @param trace [0,1,2] trace the execution of the function.
#' @param digits [integer] how to round values that are displayed in the terminal.
#'
#'
#' @details
#' \strong{formula}: \code{Y~bp(X, pattern = "111")} indicates a two breakpoints model without constrains (one intercept and three slopes).
#' \code{Y~bp(X, pattern = "101")} indicates a two breakpoints model with a plateau (one intercept and two slopes).
#' It can also contain an additional argument to specify where to initalize the breakpoint via a numeric vector: \code{Y~bp(X, pattern = "101", c(1,2))}.
#' 
#' By default initial values for the breakpoints are obtained fitting a spline model and using points where the first derivative of the spline changes sign.
#' If not enought points are found, quantiles of the response value (Y) are used instead.
#' 
#' @seealso
#' \code{\link{coef.lmbreak}} for extracting the estimated model parameters. \cr
#' \code{\link{mlmbreak}} for fitting the breakpoint model over multiple clusters of data. \cr
#' \code{\link{model.tables.lmbreak}} for extracting the breakpoints, slopes, intercept, and duration. \cr
#' \code{\link{plot.lmbreak}} for a graphical display of the fitted breakpoint model. \cr
#' 
#' @references Muggeo, V. M. R. Estimating regression models with unknown break-points.
#' Statistics in medicine 2003; 22:3055-3071.
#' 
#' @keywords models

## * lmbreak (example)
#' @examples
##'
##' ####  simulate data ####
##' set.seed(10)
##' df1 <- simBreak(c(1, 100), breakpoint = c(0,1,3,4), slope = c(1,0,-1), sigma = 0.05)
##'
##' #### fit breakpoint regression ####
##' ## broken line
##' e.lmbreak111 <- lmbreak(Y ~ bp(X, "111"), data = df1)
##' plot(e.lmbreak111)
##' summary(e.lmbreak111)
##' coef(e.lmbreak111, type = "breakpoint")
##' coef(e.lmbreak111, type = "slope")
##' coef(e.lmbreak111, type = "intercept")
##' coef(e.lmbreak111, type = "duration")
##' model.tables(e.lmbreak111)
##'
##' ## broken line with plateau
##' e.lmbreak101 <- lmbreak(Y ~ bp(X, "101"), data = df1)
##' plot(e.lmbreak101)
##' summary(e.lmbreak101)
##' model.tables(e.lmbreak101)
##' 
##' ## broken line with plateau and no intercept
##' e0.lmbreak101 <- lmbreak(Y ~ 0 + bp(X, "101"), data = df1)
##' plot(e0.lmbreak101, xlim = c(0,4))
##' summary(e0.lmbreak101)
##' model.tables(e0.lmbreak101)
##'
##' #### handle covariates ####
##' set.seed(10)
##' df2 <- df1
##' df2$tester <- cut(df1$X,c(0,0.5,2,4))
##' df2$Y <- df1$Y + as.numeric(df2$tester) + rnorm(NROW(df1), sd = 0.05)
##' 
##' eCov.lmbreak101 <- lmbreak(Y ~ tester + bp(X, "101"), data = df2)
##' summary(eCov.lmbreak101)
##' plot(eCov.lmbreak101)

## * lmbreak (code)
#' @rdname breakpoint
#' @export
lmbreak <- function(formula, data,
                    n.iter = 50, tol = 1e-3, enforce.continuity = TRUE, init.gam = TRUE, init.quantile = NULL,
                    trace = FALSE, digits = -log10(tol)){

    ## ** normalize user input
    ## *** formula
    if(!inherits(formula,"formula")){
        stop("Argument \'formula\' should be or inherit from formula. \n")
    }
    response.var <- all.vars(stats::update(formula,".~0"))
    if(length(response.var)!=1){
        stop("The argument \'formula\' should have exactly a single variable on the left hand side. \n")
    }
    terms.formula <- stats::terms(formula, specials = "bp")

    ## *** find breakpoint variable
    formula.bp <- attr(terms.formula,"specials")$bp
    if(is.null(formula.bp)){
        stop("The argument \'formula\' should contain a single variable on the right hand side indicating the breakpoints. \n",
             "Something like: Y ~ bp(X, pattern = \"111\") for three breakpoints.")
    }else if(length(formula.bp)!=1){
        stop("The argument \'formula\' currently only supports a single breakpoint variable. \n")
    }
    
    var.bp <- all.vars(attr(terms.formula,"variables")[[formula.bp+1]])
    if(length(attr(terms.formula,"term.labels"))==1){
        if(attr(terms.formula,"intercept")==1){
            terms.nobp <- stats::terms(stats::update(terms.formula, .~1))
        }else{
            terms.nobp <- stats::terms(stats::update(terms.formula, .~0))
        }
    }else{
        terms.nobp <- stats::drop.terms(terms.formula, dropx = formula.bp-1, keep.response = TRUE)
    }

    ## *** find pattern
    term.bp <- attr(terms.formula,"variables")[[formula.bp+1]]
    if(length(term.bp)<3){
        stop("The number breakpoints and possible constrained should be specified in argument \'formula\'. \n",
             "Something like: Y ~ bp(X, pattern = \"111\") for 2 breakpoints \n",
             "                Y ~ bp(X, pattern = \"101\") for 2 breakpoints with a plateau. \n")
    }
    pattern <- eval(term.bp[[3]])
    if(!is.character(pattern) || any(sapply(pattern,nchar)<2)){
        stop("The number breakpoints and possible constrained should be specified in argument \'formula\' using a character string of length at least 2. \n",
             "Something like: Y ~ bp(X, pattern = \"111\") for 2 breakpoints \n",
             "                Y ~ bp(X, pattern = \"101\") for 2 breakpoints with a plateau. \n")
    }
    n.pattern <- length(pattern)
    n.slope <- sapply(pattern, nchar)
    n.breakpoint <- n.slope-1
    nMax.breakpoint <- max(n.breakpoint)
    vec.pattern <- stats::setNames(strsplit(pattern,split="", fixed=TRUE), pattern)
    if(any("1" %in% unlist(vec.pattern) == FALSE & "0" %in% unlist(vec.pattern) == FALSE)){
        stop("When specifying the breakpoint in argument \'formula\', the pattern should only be specified with 0 or 1. \n",
             "Something like: Y ~ bp(X, pattern = \"111\") for 2 breakpoints \n",
             "                Y ~ bp(X, pattern = \"101\") for 2 breakpoints with a plateau. \n")
    }
    if(any(grepl("00",pattern,fixed=TRUE))){
        stop("When specifying the breakpoint in argument \'formula\', the pattern cannot contain two consecutive 0.")
    }

    ## *** data
    data <- as.data.frame(data)
    reserved.names <- c("breakpoint","estimate","se","df","lower","upper", paste0("Us",1:n.pattern),paste0("Vs",1:n.pattern),paste0("beta",1:n.pattern),paste0("gamma",1:n.pattern))
    if(any(names(data) %in% reserved.names)){
        txt <- names(data[names(data) %in% reserved.names])
        stop("Argument \'data\' contains reserved names: \"",paste0(txt, collapse = "\" \""),"\"\n")
    }
    if(any(all.vars(formula) %in% names(data) == FALSE)){
         stop("Argument \'data\' should contain all variables mentionned in argument \'formula\'. \n",
              "Missing variables: \"",paste0(all.vars(formula)[all.vars(formula) %in% names(data) == FALSE], collapse = "\" \""),"\".\n")
    }
    index.NA <- which(is.na(data[[response.var]]) | is.na(data[[var.bp]]))
    if(length(index.NA)>0){
        data.fit <- data[-index.NA,,drop=FALSE]
    }else{
        data.fit <- data
    }
    if(NROW(data)<3*nMax.breakpoint){
        stop("Insufficient amount of data compared to the number of breakpoints. \n",
             "Should be at least 3 datapoints per breakpoint. \n")
    }
    ## reorder
    data.fit <- data.fit[order(data.fit[[var.bp]]),]
    attr(data.fit, "min.var.bp") <- data.fit[[var.bp]][1]
    attr(data.fit, "max.var.bp") <- data.fit[[var.bp]][NROW(data.fit)]

    ## *** find initialization
    if(length(term.bp)==4){

        term.init <- cbind(as.double(eval(term.bp[[4]])))
        if(length(term.init) < nMax.breakpoint){
            stop("Insufficient number of initialization values: should be at least ",nMax.breakpoint,"\n")
        }
        breakpoint.init <- stats::setNames(lapply(n.breakpoint, function(iN){utils::combn(term.init, m = iN)}), pattern)
        breakpoint.init.rescue <- stats::setNames(lapply(n.breakpoint, function(iN){NULL}), pattern)

    }else{
        
        ## via gam
        if(init.gam){
            formulaS <- stats::as.formula(paste(response.var,"~s(",var.bp,")"))
            e.gam <- mgcv::gam(formulaS, data = data.fit)
            ## plot(e.gam)
            df.spline <- data.frame(seq(attr(data.fit, "min.var.bp"), attr(data.fit, "max.var.bp"), length.out = 1e4))
            names(df.spline) <- var.bp
            df.spline$fit <- stats::predict(e.gam, newdata = df.spline)
            breakpointGam.init <- df.spline[[1]][which(diff(sign(diff(df.spline$fit)))!=0)]
            
            if(length(breakpointGam.init)<nMax.breakpoint && identical(init.quantile,FALSE)){
                stop("Cannot initialize breakpoint position using gam. \n",
                     "Consider setting the argument \'init.quantile\' to TRUE or providing starting values, e.g. ~pb(X, \"111\", c(1,2)).")
            }
        }else{
            breakpointGam.init <- NULL
        }
        
        ## via quantiles
        if(!identical(init.quantile,FALSE)){

            if(is.null(init.quantile) || identical(init.quantile,TRUE)){
                ## +2 because 2 quantiles are removed (0,1)
                if(nMax.breakpoint==1){
                    n.init.quantile <- 2 + nMax.breakpoint + 4
                }else if(nMax.breakpoint==2){
                    n.init.quantile <- 2 + nMax.breakpoint + 3
                }else if(nMax.breakpoint==3){
                    n.init.quantile <- 2 + nMax.breakpoint + 2
                }else{
                    n.init.quantile <- 2 + nMax.breakpoint + 1
                }
            }else if(length(init.quantile)!=1 || init.quantile<(nMax.breakpoint+2) || (init.quantile %% 1 != 0)){
                stop("Argument \'nQuantile.init\' should be TRUE or an integer greater or equal to ",nMax.breakpoint+2,". \n")
            }else{
                n.init.quantile <- init.quantile
            }
            
            probs.breakpoint <- seq(0,1, length.out = n.init.quantile)[2:(n.init.quantile-1)]
            breakpointQuantile.init <- stats::quantile(data.fit[[var.bp]], probs = probs.breakpoint)
        }else{
            breakpointQuantile.init <- NULL
        }

        if(is.null(init.quantile)){
            breakpoint.init <- stats::setNames(lapply(n.breakpoint, function(iN){
                if(length(breakpointGam.init)>=iN){
                    return(utils::combn(c(breakpointGam.init), m = iN))
                }else{
                    return(utils::combn(breakpointQuantile.init, m = iN))
                }                
            }), pattern)
            breakpoint.init.rescue <- stats::setNames(lapply(n.breakpoint, function(iN){
                if(length(breakpointGam.init)>=iN){
                    return(utils::combn(breakpointQuantile.init, m = iN))
                }else{
                    return(NULL)
                }                
            }), pattern)
        }else{
            breakpoint.init <- stats::setNames(lapply(n.breakpoint, function(iN){utils::combn(c(breakpointGam.init,breakpointQuantile.init), m = iN)}), pattern)
            breakpoint.init.rescue <- stats::setNames(lapply(n.breakpoint, function(iN){NULL}), pattern)
        }
    }

    ## ** fit pattern(s)
    ls.fit <- vector(mode = "list", length = n.pattern)
    for(iPattern in 1:n.pattern){ ## iPattern <- 3

        iPrepUsVs <- .pattern2UsVs(pattern = vec.pattern[[iPattern]], terms = terms.nobp)
        ls.fit[[iPattern]] <- apply(breakpoint.init[[iPattern]], MARGIN = 2, FUN = function(iInit){
            lmbreak.fit(formula = iPrepUsVs$formula, pattern = pattern[[iPattern]], termUs.labels = iPrepUsVs$termUs.labels, signUs = iPrepUsVs$signUs, termVs.labels = iPrepUsVs$termVs.labels,
                        variable = var.bp, data = data.fit,
                        n.iter = n.iter, tol = tol, initialization = iInit, enforce.continuity = enforce.continuity,
                        trace = trace, digits = digits)
        }, simplify = FALSE)
        iFit <- do.call(rbind,lapply(ls.fit[[iPattern]], "[[","opt"))

        ## try more initialization points if no satisfactory solution
        if(length(breakpoint.init.rescue[[iPattern]])>0 && all(iFit$cv + iFit$continuity + iFit$continuity2<3)){
            ls.fit[[iPattern]] <- c(ls.fit[[iPattern]], apply(breakpoint.init.rescue[[iPattern]], MARGIN = 2, FUN = function(iInit){
                lmbreak.fit(formula = iPrepUsVs$formula, pattern = pattern[[iPattern]], termUs.labels = iPrepUsVs$termUs.labels, signUs = iPrepUsVs$signUs, termVs.labels = iPrepUsVs$termVs.labels,
                            variable = var.bp, data = data.fit,
                            n.iter = n.iter, tol = tol, initialization = iInit, enforce.continuity = enforce.continuity,
                            trace = trace, digits = digits)
            }, simplify = FALSE))
            iFit <- do.call(rbind,lapply(ls.fit[[iPattern]], "[[","opt"))
        }

        ## do not try other patterns if a satisfactory solution is found
        if(any(iFit$cv + iFit$continuity + iFit$continuity2==3)){
            break
        }
    }

    ## ** select best fitting model
    if(length(ls.fit)==1 && length(ls.fit[[1]])==1){
        out <- ls.fit[[1]][[1]]
    }else{
        lsC.fit <- unlist(ls.fit, recursive = FALSE)
        dfAll.fit <- do.call(rbind,lapply(lsC.fit, "[[","opt"))
        out <- lsC.fit[[which.max(2*dfAll.fit$cv+dfAll.fit$continuity2+dfAll.fit$R2)]]
        attr(out$opt,"all") <- dfAll.fit

        ls.allbreakpoint <- lapply(lsC.fit, function(iM){iM$breakpoint$value})
        attr(out$breakpoint,"all") <- tapply(ls.allbreakpoint,factor(lengths(ls.allbreakpoint),unique(lengths(ls.allbreakpoint))),function(iLS){do.call(rbind, iLS)})
    }

    ## ** standard error
    if(out$opt$cv){
        vcov.bp <- stats::vcov(out$model)
        beta.Us <- coef(out$model)[out$breakpoint$Us]
        beta.Vs <- coef(out$model)[out$breakpoint$Vs]
    
        term1 <- diag(vcov.bp)[out$breakpoint$Us] / beta.Us^2
        term2 <- diag(vcov.bp)[out$breakpoint$Vs] * (beta.Vs/beta.Us^2)^2
        term3 <- -2 * out$breakpoint$sign * diag(vcov.bp[out$breakpoint$Us,out$breakpoint$Vs,drop=FALSE]) * beta.Vs / beta.Us^3
        out$breakpoint$se <- sqrt(term1+term2+term3)
    }else{
        out$breakpoint$se <- NA
    }
    
    ## ** export
    if(out$opt$cv == FALSE){
        warning("The optimizer did not converge to a stable solution. \n")
    }else if(out$opt$continuity == FALSE){
        warning("The optimizer did not converge to a continuous solution (non-0 Vs terms). \n")
    }
    out$call <- match.call()
    out$args <- list(breakpoint.var = var.bp,
                     covariate = setdiff(all.vars(formula),c(var.bp,response.var)),
                     response.var = response.var,
                     pattern = vec.pattern)
    out$data <- data
    out$index.NA <- index.NA
    class(out) <- "lmbreak"
    return(out)
}



## * .lmbreak.fit
lmbreak.fit <- function(formula, pattern, termUs.labels, signUs, termVs.labels,
                        variable, data,
                        n.iter, tol, initialization, enforce.continuity,
                        trace, digits){

    ## ** initialize data
    n.breakpoint <- length(initialization)
    bp.range <- c(attr(data, "min.var.bp"),attr(data, "max.var.bp"))
    ## if(transform){
    ##     beta <- -2/diff(bp.range)
    ##     alpha <- sum(bp.range)/diff(bp.range)
    ##     ## alpha + beta*bp.range
    ## }

    ## ** loop
    Hist.breakpoint <- matrix(NA, nrow = n.iter, ncol = n.breakpoint)
    iBreakpoint <- initialization
    cv <- FALSE
    escape <- FALSE
    step <- 1
    if(trace>1){
        cat(" - initialization of the breakpoints: ",paste(round(initialization, digits), collapse =", "),"\n",
            " - iteration:", sep="")
        if(trace<3){
            cat(" ")
        }
        
    }

    for(iIter in 1:n.iter){ ## iIter <- 1
        Hist.breakpoint[iIter,] <- iBreakpoint

        ## ** update design
        iData <- model.frame.lmbreak(list(breakpoint.var = variable, breakpoint = iBreakpoint),
                                     newdata = data)

        ## ** estimate model coefficients
        iE.lm <- stats::lm(formula, data = iData)

        ## ** update breakpoint
        iCoef <- coef(iE.lm)
        iBreakpoint <- unname(step * signUs * iCoef[termVs.labels]/iCoef[termUs.labels] + Hist.breakpoint[iIter,])

        if(all(!is.na(iBreakpoint)) && any(rowSums(abs(Hist.breakpoint[1:iIter,,drop=FALSE] - matrix(iBreakpoint, nrow = iIter, ncol = n.breakpoint, byrow = TRUE)))<tol/10)){
            ## same iteration as before
            step <- step/2
            iBreakpoint <- unname(step * signUs * iCoef[termVs.labels]/iCoef[termUs.labels] + Hist.breakpoint[iIter,])
        }
        
        ## \alpha + \beta a = -1
        ## \alpha + \beta b = 1
        ## \beta = 2/(b-a) and \alpha = - (a+b)/(b-a)
    
        ## ** display
        if(trace>2){
            cat(" ",iIter," (breakpoint/Vs/step = ",paste(paste(round(iBreakpoint, digits = digits),round(iCoef[termVs.labels], digits = digits),step, sep = "/"), collapse = ", "), ") \n             ",
                sep = "")
        }else if(trace>0){
            cat("*")
        }

        ## ** cv
        iDiff <- abs(iBreakpoint-Hist.breakpoint[iIter,])
        if(any(is.na(iBreakpoint)) || iBreakpoint[1]<bp.range[1] || iBreakpoint[n.breakpoint]>bp.range[2] || is.unsorted(iBreakpoint) ){
            ## case where one breakpoint is outside the domain
            ## or when the breakpoints are not in increasing order
            if(trace>0){
                if(any(is.na(iBreakpoint))){
                    cat(" (NA)")
                }else if(is.unsorted(iBreakpoint)){
                    cat(" (breakpoints no more ordered) ")
                }else if(iBreakpoint[1]<bp.range[1] || iBreakpoint[n.breakpoint]>bp.range[2]){
                    cat(" (breakpoints outside the range of observed values) ")
                }
            }
            escape <- TRUE
            break            

        }else if(all(iDiff<tol)){
            cv <- TRUE
            if(trace>0){
                cat(" (cv)")
            }
            escape <- TRUE
            break
        }
        
    }

    if(trace>0){
        if(escape==FALSE){
            cat(" (maximum number of iterations reached)")
        }
        cat("\n")
    }

    ## ** enforce continuity
    ## continuity: Vs terms
    test.continuity <- all(abs(iCoef[termVs.labels])<tol)
    if(enforce.continuity && !test.continuity && (cv || escape == FALSE)){
        attr(iE.lm,"continuity") <- stats::lm(stats::update(formula, paste0(".~.-",paste(paste("Vs",1:n.breakpoint,sep=""),collapse="-"))), data = iData)                
    }
    ## continuity: non-consecutive breakpoint (avoid staircase)
    if(cv && test.continuity){
        ## NOTE: data is already sorted according to the breakpoint variable
        test.staircase <- sapply(iBreakpoint, function(iBreak){which.min(abs(data[[variable]]-iBreak))})
        test.continuity2 <- all(diff(test.staircase)>1)
    }else{
        test.continuity2 <- FALSE
    }

    ## ** export
    df.opt <- data.frame(initialization = NA,
                         n.iter = iIter,
                         cv = cv,
                         tol = tol,
                         pattern = pattern,
                         continuity = ifelse(is.na(test.continuity),FALSE,test.continuity),
                         continuity2 = test.continuity2,
                         diff = NA,
                         R2 = summary(iE.lm)$r.squared
                         )
    df.opt$initialization <- list(initialization)
    df.opt$diff <- list(iDiff)
    return(list(model = iE.lm,
                breakpoint = data.frame(value = iBreakpoint, Us = termUs.labels, Vs = termVs.labels, sign = signUs),
                opt = df.opt))
}

## * .pattern2UsVs
.pattern2UsVs <- function(pattern, terms){

    ## ** read input
    n.breakpoint <- length(pattern)-1

    ## ** generate formula with Us andVs
    index.1 <- which(pattern!="0")
    addU.terms <- paste0("Us",index.1-1) ## Us0 is X
    if(any(pattern[-1]=="0")){
        index.0 <- setdiff(which(pattern=="0"),0)
        for(i0 in 1:length(index.0)){ ## i0 <- 1
            iIndex <- which(index.1<index.0[i0])
            addU.terms[iIndex] <- paste0("I(",addU.terms[iIndex],"-",paste0("Us",index.0[i0]-1),")")
            index.1[iIndex] <- Inf
        }        
    }
    addV.terms <- paste0("Vs",1:n.breakpoint)
    formula <- stats::update(terms, paste0(".~.+",paste(addU.terms, collapse="+"),"+",paste(addV.terms, collapse="+")))

    terms.formula <- stats::terms(formula)
    index.gamma <- utils::tail(attr(terms.formula,"term.labels"), n.breakpoint)
    terms.beta <- utils::tail(attr(terms.formula,"term.labels"), n.breakpoint + sum(pattern == "1"))[1:sum(pattern == "1")] ## 0.1.0, 0.1.1, 1.0.1, 1.1.0, 1.1.1
    patten.beta <- paste(pattern[1:n.breakpoint],pattern[2:(n.breakpoint+1)],c(pattern,1)[3:(n.breakpoint+2)], sep=".")
    index.beta <- terms.beta[patten.beta[1] %in% c("1.1.1","1.1.0","1.0.1")+ cumsum(patten.beta != "1.0.1")]
    sign.beta <- c(1,-1)[1+(patten.beta=="1.0.1")]

    return(list(formula = formula,
                termUs.labels = index.beta,
                termVs.labels = index.gamma,
                signUs = sign.beta
                ))

}

##----------------------------------------------------------------------
### lmbreak.R ends here
