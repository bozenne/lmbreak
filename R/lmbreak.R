### lmbreak.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Apr  5 2024 (15:33) 
## Version: 
## Last-Updated: apr 19 2024 (08:50) 
##           By: Brice Ozenne
##     Update #: 1118
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
#' @name lmbreak
#'
#' @param formula a formula where the breakpoint variable appears on the right hand side
#' with an argument pattern specifying the number of breakpoints and possible constrains. See the section details.
#' @param data [data.frame] dataset
#' @param control [list] parameters to be passed to the optimizer (n.iter,tol,enforce.continuity,optimize.step).
#' See the section details of \code{\link{lmbreak.options}}.
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
##' e.lmbreak101 <- lmbreak(Y ~ bp(X, "101"), data = df1, trace = 4)
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
lmbreak <- function(formula, data, control = NULL, trace = FALSE, digits = NULL){

    ## ** normalize user input
    options <- lmbreak.options()        

    ## *** control
    name.control <- c("n.iter", "tol", "optimizer", "enforce.continuity", "optimize.step", "minR2", "init.gam","init.quantile")
    if(is.null(control)){
        control <- options[name.control]
    }else{
        if(!is.list(control)){
            stop("Argument \'control\' should be a list. \n")
        }
        if(is.null(names(control))){
            stop("Argument \'control\' should be a named list. \n")
        }
        if(any(names(control) %in% name.control == FALSE)){
            stop("Incorrect element(s) in argument control: \"",paste(names(control)[names(control) %in% name.control == FALSE], collapse = "\", \""),"\" \n",
                 "Valid names: \"",paste(setdiff(name.control, names(control)), collapse = "\", \""),"\" \n")
        }
        control[setdiff(name.control,names(control))] <- options[setdiff(name.control,names(control))]
        if("tol" %in% names(control)){
            if(length(control$tol)==1){
                control$tol <- rep(control$tol,2)
            }else if(length(control$tol)!=2){
                stop("Element \"tol\" in argument \'control\' must be have length 2. \n")
            }
        }
        if("optimizer" %in% names(control) && control$optimizer %in% c("Muggeo","nlminb","BFGS","L-BFGS-B") == FALSE){
            stop("Incorrect argument \'control\': element \"optimizer\" can only be \"Muggeo\", \"nlminb\" (grid search), or \"BFGS\" (subgradient descent). \n")
        }
    }

    if(is.null(digits)){
        digits <- -log(control$tol[1])
    }

    ## *** formula
    if(!inherits(formula,"formula")){
        stop("Argument \'formula\' should be or inherit from formula. \n")
    }
    var.response <- all.vars(stats::update(formula,".~0"))
    if(length(var.response)!=1){
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
    if(all(all.vars(stats::delete.response(terms.formula)) %in% var.bp)){ ## no covariate
        if(attr(terms.formula,"intercept")==1){
            terms.nobp <- stats::terms(stats::update(terms.formula, .~1))
        }else{
            terms.nobp <- stats::terms(stats::update(terms.formula, .~0))
        }
    }else{ ## covariates
        terms.nobp <- stats::drop.terms(terms.formula, dropx = formula.bp-1, keep.response = TRUE)
    }
    if(any(attr(terms.formula,"order")>1) && any(attr(terms.formula,"factor")[formula.bp,attr(terms.formula,"order")>1]>0)){ ## check no interaction with breakpoint
            stop("The argument \'formula\' should not contain interaction(s) with the breakpoint variable. \n")
    }

    ## *** find pattern
    term.bp <- attr(terms.formula,"variables")[[formula.bp+1]]
    if(length(term.bp)<3){
        stop("The number breakpoints and possible constrained should be specified in argument \'formula\'. \n",
             "Something like: Y ~ bp(X, pattern = \"111\") for 2 breakpoints \n",
             "                Y ~ bp(X, pattern = \"101\") for 2 breakpoints with a plateau. \n")
    }
    pattern <- eval(term.bp[[3]])
    if(is.list(pattern)){
        if(any(!sapply(pattern,is.character)) || any(unlist(lapply(pattern,nchar))<2)){
            stop("The number breakpoints and possible constrained should be specified in argument \'formula\' using character strings of length at least 2. \n",
                 "Something like: Y ~ bp(X, pattern = \"111\") for 2 breakpoints \n",
                 "                Y ~ bp(X, pattern = \"101\") for 2 breakpoints with a plateau. \n")
        }

        n.pattern <- sum(lengths(pattern))
        hierarchy.pattern <- unlist(lapply(1:length(pattern), function(iP){rep(iP,length(pattern[[iP]]))}))
        pattern <- unlist(pattern)
    }else{
        if(!is.character(pattern) || any(nchar(pattern)<2)){
            stop("The number breakpoints and possible constrained should be specified in argument \'formula\' using character strings of length at least 2. \n",
                 "Something like: Y ~ bp(X, pattern = \"111\") for 2 breakpoints \n",
                 "                Y ~ bp(X, pattern = \"101\") for 2 breakpoints with a plateau. \n")
        }
        n.pattern <- length(pattern)
        hierarchy.pattern <- rep(1,n.pattern)
    }

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

    index.NA <- which(is.na(data[[var.response]]) | is.na(data[[var.bp]]))
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
    attr(data.fit, "min.diff.bp") <- min(diff(data.fit[[var.bp]]))

    ## ** initialization

    ## *** breakpoint
    if(length(term.bp)==4){

        term.init <- cbind(as.double(eval(term.bp[[4]])))
        if(length(term.init) < nMax.breakpoint){
            stop("Insufficient number of initialization values: should be at least ",nMax.breakpoint,"\n")
        }
        breakpoint.init <- lapply(n.breakpoint, function(iN){combn2(term.init, m = iN, name = "manual")})
        name.order <- c("manual" = "manual")

    }else{

        breakpoint.prepare <- list()
        if(control$init.gam<=0 & control$init.quantile<=0){
            stop("If no initialization value are specified for the breakpoints, argument \'init.gam\' and \'init.quantile\' cannot simultaneously be FALSE. \n")
        }
        name.order <- union(c("gam","quantile")[which.max(c(control$init.gam,control$init.quantile))], c("gam","quantile"))
        if(control$init.gam!=control$init.quantile){
            name.order <- stats::setNames(name.order,name.order)
        }else{
            name.order <- stats::setNames(name.order,c("auto","auto"))
        }

        ## via gam
        if(control$init.gam>0){
            breakpoint.prepare$gam <- .init_gam(formula = terms.nobp, data = data.fit, var.bp = var.bp)
            if(length(breakpoint.prepare$gam)<nMax.breakpoint && control$init.quantile<=0){
                stop("Cannot initialize breakpoint position using gam. \n",
                     "Consider setting the argument \'init.quantile\' to TRUE or providing starting values, e.g. ~pb(X, \"111\", c(1,2)).")
            }
        }

        ## via quantiles
        ## if((control$optimizer == "Muggeo") && (control$init.quantile>0)){
        if(control$init.quantile>0){
            breakpoint.prepare$quantile <- .init_quantile(init = control$init.quantile, data = data.fit, var.bp = var.bp, n.breakpoint = nMax.breakpoint)
        }

        ## if(control$optimizer == "Muggeo"){
        breakpoint.init <- lapply(n.breakpoint, function(iN){ ## iN <- 3
            cbind(combn2(breakpoint.prepare[[name.order[1]]], m = iN, name = names(name.order)[1]),
                  combn2(breakpoint.prepare[[name.order[2]]], m = iN, name = names(name.order)[2]))
        })
        ## }else if(control$optimizer %in% c("nlminb","L-BFGS-B")){
        ##     breakpoint.init <- lapply(n.breakpoint, function(iN){ ## iN <- 1
        ##         if(length(breakpoint.prepare$gam)>=iN){
        ##             iOut <- sort(breakpoint.prepare$gam[order(abs(attr(breakpoint.prepare$gam,"diff")), decreasing = TRUE)[1:iN]])
        ##         }else{
        ##             iOut <- quantile(data.fit[[var.bp]], seq(0, 1, length.out = iN+2)[2:(iN+1)])
        ##         }
        ##         return(cbind(iOut))
        ##     })
        ## }
    }
    breakpoint.init <- stats::setNames(breakpoint.init, pattern)

    ## *** pattern
    patternUsVs <- stats::setNames(lapply(1:n.pattern, function(iP){.pattern2UsVs(pattern = vec.pattern[[iP]], terms = terms.nobp)}), pattern)

    ## *** monitor cv
    grid.fit <- data.frame(index = NA,
                           hierarchy = NA,
                           pattern = unlist(mapply(x = pattern, times = sapply(breakpoint.init,NCOL),rep, SIMPLIFY = FALSE)),                           
                           init = unlist(lapply(breakpoint.init, function(iM){1:NCOL(iM)})),
                           init.type = factor(unlist(lapply(breakpoint.init, colnames)), names(name.order)),
                           cv = NA, continuity = NA, continuity2 = NA, R2 = NA
                           )
    rownames(grid.fit) <- NULL
    pattern2hierarchy <- stats::setNames(paste0(hierarchy.pattern,".",sapply(patternUsVs, "[[", "n.param")), pattern)
    grid.fit$hierarchy <- as.numeric(factor(pattern2hierarchy, unique(pattern2hierarchy))[grid.fit$pattern])
    grid.fit <- grid.fit[order(grid.fit$hierarchy,grid.fit$init.type,grid.fit$init),,drop=FALSE]
    grid.fit$index <- 1:NROW(grid.fit)

    pattern.possibleStop <- stats::setNames(c(!duplicated(paste0(hierarchy.pattern,".",sapply(patternUsVs, "[[", "n.param")))[-1],TRUE), pattern)
    ## check for convergence before going down in hierarchy or changing of initialization method
    grid.fit$check.cv <- rev(!duplicated(rev(grid.fit$hierarchy))) | (rev(!duplicated(rev(paste0(grid.fit$init.type,grid.fit$hierarchy)))))

    ## ** fit pattern(s)
    ls.fit <- vector(mode = "list", length = NROW(grid.fit))
    for(iGrid in 1:NROW(grid.fit)){ ## iGrid <- 3
        ## select
        iPattern <- grid.fit[iGrid,"pattern"]
        iInit <- breakpoint.init[[iPattern]][,grid.fit[iGrid,"init"]]

        ## display
        if(trace>0){
            if(trace>=2){
                cat(" - (",iGrid,") pattern ",iPattern,", initialization ",paste0(round(iInit, digits[1]),collapse = ", "),if(trace>=3){" "},if(trace>=4){"\n"}, sep="")
            }else if(trace>=1){
                if((iGrid == 1 || iPattern!=grid.fit$pattern[iGrid-1])){
                    cat("\n - Pattern ",iPattern,": ",iGrid, sep = "")
                }else{
                    cat(", ",iGrid)
                }
            }else{
                cat("*")
            }
        }

        ## optimize
        if(control$optimizer=="Muggeo"){
            ls.fit[[iGrid]] <- lmbreak.fitMuggeo(formula = patternUsVs[[iPattern]]$formula, pattern = iPattern,
                                                 Us.label = patternUsVs[[iPattern]]$Us.label, Us.sign = patternUsVs[[iPattern]]$Us.sign, Vs.label = patternUsVs[[iPattern]]$Vs.label,
                                                 var.response = var.response, var.bp = var.bp, data = data.fit,
                                                 n.iter = control$n.iter, tol = control$tol, initialization = iInit, enforce.continuity = control$enforce.continuity, optimize.step = control$optimize.step,
                                                 trace = trace-1.999999, digits = digits)
        }else if(control$optimizer=="nlminb"){
            ls.fit[[iGrid]] <- lmbreak.fitNLMINB(formula = patternUsVs[[iPattern]]$formula, pattern = iPattern,
                                             Us.label = patternUsVs[[iPattern]]$Us.label, Us.sign = patternUsVs[[iPattern]]$Us.sign, Vs.label = patternUsVs[[iPattern]]$Vs.label,
                                             var.response = var.response, var.bp = var.bp, data = data.fit,
                                             n.iter = control$n.iter, tol = control$tol, initialization = iInit, enforce.continuity = control$enforce.continuity, 
                                             trace = trace-1.999999, digits = digits)
        }else if(control$optimizer %in% c("BFGS","L-BFGS-B")){
            ls.fit[[iGrid]] <- lmbreak.fitBFGS(formula = patternUsVs[[iPattern]]$formula, pattern = iPattern, transform = control$optimizer == "BFGS",
                                               Us.label = patternUsVs[[iPattern]]$Us.label, Us.sign = patternUsVs[[iPattern]]$Us.sign, Vs.label = patternUsVs[[iPattern]]$Vs.label,
                                               var.response = var.response, var.bp = var.bp, data = data.fit,
                                               n.iter = control$n.iter, tol = control$tol, initialization = iInit, enforce.continuity = control$enforce.continuity, 
                                               trace = trace-1.999999, digits = digits)
        }
        grid.fit[iGrid,c("cv","continuity","continuity2","R2")] <- as.double(ls.fit[[iGrid]]$opt[c("cv","continuity","continuity2","R2")])

        ## do not try other patterns or other initializations if a satisfactory solution is found
        if(grid.fit[iGrid,"check.cv"] && any(rowSums(grid.fit[1:iGrid,c("cv","continuity","continuity2")])==3)){
            break
        }
    }
    if(trace>0){
        cat("\n")
    }

    ## ** select best fitting model
    if(iGrid==1){
        out <- ls.fit[[1]]
    }else{
        out <- ls.fit[[which.max(2*(grid.fit[1:iGrid,"cv"] & grid.fit[1:iGrid,"R2"]>control$minR2) + grid.fit[1:iGrid,"continuity2"]+grid.fit[1:iGrid,"R2"])]]
        attr(out$opt,"all") <- do.call(rbind,lapply(ls.fit[1:iGrid],"[[","opt"))

        ls.allbreakpoint <- lapply(ls.fit[1:iGrid], function(iM){iM$breakpoint$value})
        attr(out$breakpoint,"all") <- tapply(ls.allbreakpoint,factor(lengths(ls.allbreakpoint),unique(lengths(ls.allbreakpoint))),function(iLS){do.call(rbind, iLS)})
    }

    ## ** standard error
    ## if(out$opt$cv){
    ##     vcov.bp <- stats::vcov(out$model)
    ##     beta.Us <- coef(out$model)[out$breakpoint$Us]
    ##     beta.Vs <- coef(out$model)[out$breakpoint$Vs]
    
    ##     term1 <- diag(vcov.bp)[out$breakpoint$Us] / beta.Us^2
    ##     term2 <- diag(vcov.bp)[out$breakpoint$Vs] * (beta.Vs/beta.Us^2)^2
    ##     term3 <- -2 * out$breakpoint$sign * diag(vcov.bp[out$breakpoint$Us,out$breakpoint$Vs,drop=FALSE]) * beta.Vs / beta.Us^3
    ##     out$breakpoint$se <- sqrt(term1+term2+term3)
    ## }else{
        ## out$breakpoint$se <- NA
    ## }
    
    ## ** export
    if(out$opt$cv == FALSE){
        warning("The optimizer did not converge to a stable solution. \n")
    }else if(out$opt$continuity == FALSE){
        warning("The optimizer did not converge to a continuous solution (non-0 Vs terms). \n")
    }
    out$call <- match.call()
    out$args <- list(breakpoint.var = var.bp,
                     covariate = setdiff(all.vars(formula),c(var.bp,var.response)),
                     response.var = var.response,
                     pattern = vec.pattern)
    out$data <- data
    out$index.NA <- index.NA
    class(out) <- "lmbreak"
    return(out)
}

## * .pattern2UsVs (Us,Vs)
.pattern2UsVs <- function(pattern, terms){

    ## *** read input
    n.breakpoint <- length(pattern)-1

    ## *** generate formula with Us andVs
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
    attr(index.gamma,"index") <- utils::tail(1:length(attr(terms.formula,"term.labels")), n.breakpoint)
    terms.beta <- utils::tail(attr(terms.formula,"term.labels"), n.breakpoint + sum(pattern == "1"))[1:sum(pattern == "1")] ## 0.1.0, 0.1.1, 1.0.1, 1.1.0, 1.1.1
    patten.beta <- paste(pattern[1:n.breakpoint],pattern[2:(n.breakpoint+1)],c(pattern,1)[3:(n.breakpoint+2)], sep=".")
    index.beta <- terms.beta[patten.beta[1] %in% c("1.1.1","1.1.0","1.0.1")+ cumsum(patten.beta != "1.0.1")]
    sign.beta <- c(1,-1)[1+(patten.beta=="1.0.1")]

    ## *** export
    return(list(formula = formula,
                Us.label = index.beta,
                Vs.label = index.gamma,
                Us.sign = sign.beta,
                n.param = length(attr(terms.formula,"term.labels"))
                ))

}

## * Muggeo fitter 
## ** .init_gam (initialization)
.init_gam <- function(formula, data, var.bp, n.points = 1e3){

    ## *** fit spline
    formulaS <- stats::update(formula,paste0(".~.+1+s(",var.bp,")")) ## always restore intercept (otherwise mgcv does not estimate the spline)
    e.gam <- mgcv::gam(formulaS, data = data)
    ## plot(e.gam)

    ## *** extract fitted spline
    df.spline <- data.frame(seq(attr(data, "min.var.bp"), attr(data, "max.var.bp"), length.out = n.points))
    names(df.spline) <- var.bp

    var.nobp <- all.vars(stats::delete.response(formula))
    if(length(var.nobp)==0){
        df.spline$fit <- stats::predict(e.gam, newdata = df.spline)
    }else{
        df.spline$fit <- stats::predict(e.gam, newdata = cbind(df.spline, as.list(data[1,var.nobp,drop=FALSE])))
    }
    ## *** identify change of slope
    diff.fit <- diff(df.spline$fit)
    index.change <- which(diff(sign(diff.fit))!=0)
    out <- df.spline[[1]][index.change]
    attr(out,"diff") <- diff.fit[index.change]

    ## *** export
    return(out)
}

## ** .init_quantile (initialization)
.init_quantile <- function(init, data, var.bp, n.breakpoint){

    ## *** decide on the number of quantiles
    if(init>0){
        if(n.breakpoint==1){
            n.quantile <- n.breakpoint + 4
        }else if(n.breakpoint==2){
            n.quantile <- n.breakpoint + 3
        }else if(n.breakpoint==3){
            n.quantile <- n.breakpoint + 2
        }else if(n.breakpoint>=4){
            n.quantile <- n.breakpoint + 1
        }
    }else{
        if(length(init)!=1 || init<n.breakpoint || (init %% 1 != 0)){
            stop("Argument \'init.quantile\' should be TRUE or an integer greater or equal to ",n.breakpoint,". \n")
        }else{
            n.quantile <- init
        }
    }
            
    ## *** evaluate quantiles
    ## n.quantile+2 because 2 quantiles are removed (0,1) --> 2, n.quantile+2-1=n.quantile+1
    probs.breakpoint <- seq(0,1, length.out = n.quantile+2)[2:(n.quantile+1)]
    out <- stats::quantile(data[[var.bp]], probs = probs.breakpoint)

    ## *** export
    return(out)
}



## ** .lmbreak.fitMuggeo
## fit using muggeo algo
lmbreak.fitMuggeo <- function(formula, pattern, Us.label, Us.sign, Vs.label,
                              var.response, var.bp, data,
                              n.iter, tol, optimize.step, initialization, enforce.continuity,
                              trace, digits){

    ## *** prepare
    n.breakpoint <- length(initialization)
    bp.range <- c(attr(data, "min.var.bp"),attr(data, "max.var.bp"))
    formula.noVs <-  stats::formula(stats::drop.terms(stats::terms(formula), attr(Vs.label,"index"), keep.response = TRUE))

    fit.fun <- function(step){ ## step <- 0.5
        for(iPoint in 1:n.breakpoint){ ## iPoint <- 2
            iData[[paste0("Us",iPoint)]] <- pmax(0,iData[[var.bp]] - (step * iChange[iPoint] + Hist.breakpoint[iIter,iPoint]))
        }
        iiLM <- stats::lm.fit(iData[[var.response]], x = stats::model.matrix(formula,iData))
        ## sum of absolute deviations 
        return(sum(abs(iiLM$residuals)))        
        ## return(sum(iiLM$residuals^2))
    }

    ## *** loop
    Hist.breakpoint <- matrix(NA, nrow = n.iter, ncol = n.breakpoint)        
    iBreakpoint <- initialization
    cv <- FALSE
    escape <- FALSE
    step <- 1

    for(iIter in 1:n.iter){ ## iIter <- 1

        if(n.iter>0){ ## handle special case of no iteration
            Hist.breakpoint[iIter,] <- iBreakpoint
        }

        ## **** update design
        iData <- model.frame.lmbreak(list(breakpoint.var = var.bp, breakpoint = iBreakpoint),
                                     newdata = data)
        
        ## **** estimate model coefficients
        iE.lm <- stats::lm.fit(iData[[var.response]], x = stats::model.matrix(formula,iData))

        ## **** update breakpoint
        iCoef <- iE.lm$coef
        if(n.iter==0){
            iDiff <- rep(0, length(iBreakpoint))
            break
        }
        iChange <- unname(Us.sign * iCoef[Vs.label]/iCoef[Us.label])
        if(all(!is.infinite(iChange)) && all(!is.na(iChange)) && optimize.step>0 && (optimize.step>=1 || fit.fun(1)>=fit.fun(0))){
            optim.step <- stats::optimize(f = fit.fun, lower = 0, upper = 1, tol = min(tol))
            if(optim.step$objective < fit.fun(0)){
                iStep <- optim.step$minimum
            }else{
                iStep <- 0
            }
        }else{
            iStep <- 1
        }
        iBreakpoint <- iStep * iChange + Hist.breakpoint[iIter,]
    
        ## **** display
        if(trace>=2){
            cat("    iteration ",iIter," (step=",round(iStep,digits=digits),"): breakpoint/Vs = ",paste(paste(round(iBreakpoint, digits = digits),round(iCoef[Vs.label], digits = digits), sep = "/"),
                                                                                                       collapse = ", "), ") \n",
                sep = "")
        }else if(trace>=1){
            cat("*")
        }

        ## **** cv
        iDiff <- abs(iBreakpoint-Hist.breakpoint[iIter,])
        if(any(is.na(iBreakpoint)) || iBreakpoint[1]<bp.range[1] || iBreakpoint[n.breakpoint]>bp.range[2] || is.unsorted(iBreakpoint) ){
            ## case where one breakpoint is outside the domain
            ## or when the breakpoints are not in increasing order
            if(any(is.na(iBreakpoint))){
                if(trace>=2){
                    cat("    no convergence: some breakpoint were estimated to be NAs\n")
                }else if(trace>0){
                    cat(" (no cv: NA breakpoint)")
                }
            }else if(is.unsorted(iBreakpoint)){
                if(trace>=2){
                    cat("    no convergence: breakpoint were more ordered\n")
                }else if(trace>0){
                    cat(" (no cv: unordered breakpoint)")
                }
            }else if(iBreakpoint[1]<bp.range[1] || iBreakpoint[n.breakpoint]>bp.range[2]){
                if(trace>=2){
                    cat("    no convergence: breakpoints outside the range of observed values\n")
                }else if(trace>0){
                    cat(" (no cv: breakpoint range)")
                }
            }
            escape <- TRUE
            break            

        }else if(all(iDiff<tol[1])){
            cv <- TRUE
            if(trace>=2){
                cat("    convergence\n")
            }else if(trace>0){
                cat(" (cv)")
            }
            escape <- TRUE
            break
        }
        
    }

    if(trace>0){
        if(escape==FALSE){
            if(trace>=2){
                cat("    no convergence: maximum number of iterations reached\n")
            }else if(trace>0){
                cat(" (no cv: n.iter)")
            }            
        }
        cat("\n")
    }

    ## *** enforce continuity
    ## continuity: Vs terms
    iE.lm <- stats::lm(formula, data = iData)
    if(is.na(tol[2])){
        test.continuity <- TRUE        
    }else{
        test.continuity <- all(abs(iCoef[Vs.label])<tol[2]) && (min(c(iBreakpoint,Inf),na.rm=TRUE) - bp.range[1] > tol[1]) && (bp.range[2] - max(c(iBreakpoint,-Inf),na.rm=TRUE) > tol[1])
        if(enforce.continuity && !test.continuity && (cv || escape == FALSE)){
            attr(iE.lm,"continuity") <- stats::lm(formula.noVs, data = iData)                
        }
    }

    ## continuity: non-consecutive breakpoint (avoid staircase)
    if(cv && test.continuity){
        ## NOTE: data is already sorted according to the breakpoint variable
        test.staircase <- sapply(iBreakpoint, function(iBreak){which.min(abs(data[[var.bp]]-iBreak))})
        test.continuity2 <- all(diff(test.staircase)>1)
    }else{
        test.continuity2 <- FALSE
    }
    

    ## *** export
    df.opt <- data.frame(initialization = NA,
                         n.iter = iIter,
                         cv = cv,
                         tol = NA,
                         pattern = pattern,
                         continuity = ifelse(is.na(test.continuity),FALSE,test.continuity),
                         continuity2 = test.continuity2,
                         diff = NA,
                         R2 = summary(iE.lm)$r.squared
                         )
    df.opt$tol <- list(tol)
    df.opt$initialization <- list(initialization)
    df.opt$diff <- list(iDiff)
    return(list(model = iE.lm,
                breakpoint = data.frame(value = iBreakpoint, Us = Us.label, Vs = Vs.label, sign = Us.sign),
                opt = df.opt))
}

## * grid search fitter 
## ** lmbreak.fitNLMINB
lmbreak.fitNLMINB <- function(formula, pattern, Us.label, Us.sign, Vs.label,
                          var.response, var.bp, data,
                          n.iter, tol, initialization, enforce.continuity,
                          trace, digits){

    ## *** prepare
    n.breakpoint <- length(Us.label)
    formula.noVs <-  stats::formula(stats::drop.terms(stats::terms(formula), attr(Vs.label,"index"), keep.response = TRUE))
    bp.range <- c(attr(data,"min.var.bp"),attr(data,"max.var.bp"))
    bp.min.diff <- attr(data,"min.diff.bp")
    Id <- diag(1,NROW(data))
    data$Us0 <- data[[var.bp]]

    ## *** grid search
    res.optim <- stats::nlminb(start = initialization,
                               objective = .RRS.lmbreak, transform = FALSE, formula = formula.noVs, data = data, var.bp = var.bp, var.response = var.response, 
                               lower = rep(bp.range[1] + bp.min.diff, n.breakpoint),
                               upper = rep(bp.range[2]- bp.min.diff, n.breakpoint),
                               control = list(abs.tol = tol[1], eval.max = 10*n.iter, iter.max = 7.5*n.iter))
    iBreakpoint <- res.optim$par
    cv <- res.optim$convergence==0

    ## *** check fit
    for(iPoint in 1:n.breakpoint){ ## iPoint <- 2
        data[[paste0("Vs",iPoint)]] <- as.numeric(data[[var.bp]] > iBreakpoint[iPoint])
        data[[paste0("Us",iPoint)]] <- (data[[var.bp]] - iBreakpoint[iPoint])*(data[[var.bp]] > iBreakpoint[iPoint])
    }        

    iE.lm <- stats::lm(formula, data = data)
    iCoef <- coef(iE.lm)

    if(is.na(tol[2])){
        test.continuity <- TRUE        
    }else{
        if(any(is.na(iCoef))){
            test.continuity <- FALSE
        }else{
            test.continuity <- all(abs(iCoef[Vs.label])<tol[2]) && (min(c(iBreakpoint,Inf),na.rm=TRUE) - bp.range[1] > tol[1]) && (bp.range[2] - max(c(iBreakpoint,-Inf),na.rm=TRUE) > tol[1])
        }
        if(enforce.continuity && !test.continuity && cv){
            attr(iE.lm,"continuity") <- stats::lm(formula.noVs, data = data)                
        }
    }

    ## continuity: non-consecutive breakpoint (avoid staircase)
    if(cv && test.continuity){
        ## NOTE: data is already sorted according to the breakpoint variable
        test.staircase <- sapply(iBreakpoint, function(iBreak){which.min(abs(data[[var.bp]]-iBreak))})
        test.continuity2 <- all(diff(test.staircase)>1)
    }else{
        test.continuity2 <- FALSE
    }

    ## *** export
    df.opt <- data.frame(initialization = NA,
                         n.iter = res.optim$iteration,
                         cv = cv,
                         tol = NA,
                         pattern = pattern,
                         continuity = ifelse(is.na(test.continuity),FALSE,test.continuity),
                         continuity2 = test.continuity2,
                         diff = NA,
                         R2 = summary(iE.lm)$r.squared
                         )
    df.opt$tol <- list(tol)
    df.opt$initialization <- list(initialization)
    df.opt$diff <- list(NULL)
    return(list(model = iE.lm,
                breakpoint = data.frame(value = iBreakpoint, Us = Us.label, Vs = Vs.label, sign = Us.sign),
                opt = df.opt))
}


## * BFGS with subgradient
## ** lmbreak.fitBFGS
lmbreak.fitBFGS <- function(formula, pattern, Us.label, Us.sign, Vs.label, transform,
                            var.response, var.bp, data,
                            n.iter, tol, initialization, enforce.continuity,
                            trace, digits){

    
    ## *** prepare
    n.breakpoint <- length(Us.label)
    formula.noVs <-  stats::formula(stats::drop.terms(stats::terms(formula), attr(Vs.label,"index"), keep.response = TRUE))
    bp.range <- c(attr(data,"min.var.bp"),attr(data,"max.var.bp"))
    bp.min.diff <- attr(data,"min.diff.bp")
    Id <- diag(1,NROW(data))
    data$Us0 <- data[[var.bp]]

    dX.skeleton <-lapply(1:n.breakpoint, function(iPoint){ ## iPoint <- 2
        ddata <- data
        ddata[,paste0("Us",0:n.breakpoint)] <- matrix(as.numeric(0:n.breakpoint == iPoint), byrow = TRUE, nrow = NROW(data), ncol = n.breakpoint+1)
        iX <- stats::model.matrix(formula.noVs, ddata)
        if("(Intercept)" %in% colnames(iX)){
            iX[,"(Intercept)"] <- 0
        }
        return(iX)
    })

    ## *** gradient descent
    if(transform){
        initialization.trans <- transformPsi(initialization, min = bp.range[1], max = bp.range[2], mindiff = bp.min.diff, jacobian = FALSE)
        attr(transform,"range") <- bp.range
        attr(transform,"mindiff") <- bp.min.diff

        ## SANITY CHECK
        ## .RRS.lmbreak(initialization, transform = FALSE, formula = formula.noVs, data = data, var.bp = var.bp, var.response = var.response)
        ## .RRS.lmbreak(initialization.trans, transform = transform, formula = formula.noVs, data = data, var.bp = var.bp, var.response = var.response)
        
        ## numDeriv::jacobian(.RRS.lmbreak, initialization.trans, transform = transform, formula = formula.noVs, data = data, var.bp = var.bp, var.response = var.response)
        ## .score.lmbreak(initialization.trans, transform = transform, formula = formula.noVs, data = data, var.bp = var.bp, var.response = var.response, dX.skeleton = dX.skeleton, tol = tol[1])

        res.optim <- stats::optim(par = initialization.trans, method = "BFGS",
                                  fn = function(psi){.RRS.lmbreak(psi, transform = transform, formula = formula.noVs, data = data, var.bp = var.bp, var.response = var.response)},
                                  gr = function(psi){.score.lmbreak(psi, transform = transform, formula = formula.noVs, data = data, var.bp = var.bp, dX.skeleton = dX.skeleton, var.response = var.response, tol = tol[1])},
                                  control = list(maxit = 5*n.iter))
        iBreakpoint <- backtransformPsi(res.optim$par, min = bp.range[1], max = bp.range[2], mindiff = bp.min.diff, jacobian = FALSE)
    }else{       
         
        ## SANITY CHECK
        ## numDeriv::jacobian(.RRS.lmbreak, initialization, transform = FALSE, formula = formula.noVs, data = data, var.bp = var.bp, var.response = var.response)
        ## .score.lmbreak(initialization, transform = FALSE, formula = formula.noVs, data = data, var.bp = var.bp, var.response = var.response, dX.skeleton = dX.skeleton, tol = tol[1])

        res.optim <- stats::optim(par = initialization, method = "L-BFGS-B",
                                  fn = function(psi){.RRS.lmbreak(psi, transform = FALSE, formula = formula.noVs, data = data, var.bp = var.bp, var.response = var.response)},
                                  gr = function(psi){.score.lmbreak(psi, transform = FALSE, formula = formula.noVs, data = data, var.bp = var.bp, dX.skeleton = dX.skeleton, var.response = var.response, tol = tol[1])},
                                  lower = rep(bp.range[1], n.breakpoint) + bp.min.diff/2,
                                  upper = rep(bp.range[2], n.breakpoint) - bp.min.diff/2,
                                  control = list(maxit = 5*n.iter))
        iBreakpoint <- res.optim$par
    }
    cv <- res.optim$convergence==0

    ## *** uncertainty quantification
    ## SANITY CHECK
    ## GS <- numDeriv::jacobian(.score.lmbreak, iBreakpoint, transform = FALSE, formula = formula.noVs, data = data, var.bp = var.bp, var.response = var.response, dX.skeleton = dX.skeleton, tol = tol[1])
    if(cv && all(!is.na(iBreakpoint))){
        iHessian <- .hessian.lmbreak(iBreakpoint, transform = FALSE, formula = formula.noVs, data = data, var.bp = var.bp, var.response = var.response, dX.skeleton = dX.skeleton, tol = tol[1])

        if(all(!is.na(iHessian)) && all(abs(iHessian)<1e10) && (det(iHessian)>tol[1])){
            iBreakpoint.se <- sqrt(diag(solve(iHessian)))
        }else{
            iBreakpoint.se <- rep(NA, n.breakpoint)
        }
    }else{
        iBreakpoint.se <- rep(NA, n.breakpoint)
    }
    
    ## *** check fit
    for(iPoint in 1:n.breakpoint){ ## iPoint <- 2
        data[[paste0("Vs",iPoint)]] <- as.numeric(data[[var.bp]] > iBreakpoint[iPoint])
        data[[paste0("Us",iPoint)]] <- (data[[var.bp]] - iBreakpoint[iPoint])*(data[[var.bp]] > iBreakpoint[iPoint])
    }        

    iE.lm <- stats::lm(formula, data = data)
    iCoef <- coef(iE.lm)

    ## test Vs terms
    if(is.na(tol[2])){
        test.continuity <- TRUE        
    }else{
        if(any(is.na(iCoef))){
            test.continuity <- FALSE
        }else{
            test.continuity <- all(abs(iCoef[Vs.label])<tol[2])
        }
        if(enforce.continuity && !test.continuity && cv){
            attr(iE.lm,"continuity") <- stats::lm(formula.noVs, data = data)                
        }
    }

    ## continuity: non-consecutive breakpoint (avoid staircase)
    if(cv && test.continuity){
        test.BPinRange <- (min(c(iBreakpoint,Inf),na.rm=TRUE) - bp.range[1] > tol[1]) && (bp.range[2] - max(c(iBreakpoint,-Inf),na.rm=TRUE) > tol[1])
        if(test.BPinRange==FALSE){
            test.continuity2 <- FALSE
        }else{
            ## NOTE: data is already sorted according to the breakpoint variable
            test.staircase <- sapply(iBreakpoint, function(iBreak){which(sign(data[[var.bp]]-iBreak)==1)[1]})
            test.continuity2 <- all(diff(test.staircase)>1)
        }        
    }else{
        test.continuity2 <- FALSE
    }

    ## *** export
    df.opt <- data.frame(initialization = NA,
                         n.iter = res.optim$counts[1],
                         cv = cv,
                         tol = NA,
                         pattern = pattern,
                         continuity = ifelse(is.na(test.continuity),FALSE,test.continuity),
                         continuity2 = test.continuity2,
                         diff = NA,
                         R2 = summary(iE.lm)$r.squared
                         )
    df.opt$tol <- list(tol)
    df.opt$initialization <- list(initialization)
    df.opt$diff <- list(NULL)
    return(list(model = iE.lm,
                breakpoint = data.frame(value = iBreakpoint, se = iBreakpoint.se, Us = Us.label, Vs = Vs.label, sign = Us.sign),
                opt = df.opt))
}



##----------------------------------------------------------------------
### lmbreak.R ends here
