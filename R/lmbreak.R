### lmbreak.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Apr  5 2024 (15:33) 
## Version: 
## Last-Updated: jul 18 2024 (14:38) 
##           By: Brice Ozenne
##     Update #: 1437
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
#' @param range [numeric vector of length 2] minimum and maximum outcome value at the breakpoints.
#' @param pattern [character or character vector] alternative way to specify the pattern.
#' @param start [numeric vector] starting values when estimating the breakpoints.
#' @param se [logical] should the uncertainty about the breakpoint position be quantified (EXPERIMENTAL).
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
##' grid111 <- data.frame(tester = "(0,0.5]", X = coef(e.lmbreak111, "breakpoint.range"))
##' predict(e.lmbreak111, newdata = grid111)
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
##' eCov.lmbreak101 <- lmbreak(Y ~ tester + bp(X, "101"), data = df2, start = c(1.03474573, 2.91979538))
##' summary(eCov.lmbreak101)
##' plot(eCov.lmbreak101)
##' 
##' summary(eCov.lmbreak101, continuity = FALSE)
##' plot(eCov.lmbreak101, continuity = FALSE)
##'
##' gridCov101 <- data.frame(tester = "(0,0.5]",
##'                          X = coef(eCov.lmbreak101, "breakpoint.range"))
##' predict(eCov.lmbreak101, newdata = gridCov101)
##' predict(eCov.lmbreak101, newdata = gridCov101, continuity = FALSE, extrapolate = TRUE)

## * lmbreak (code)
#' @rdname breakpoint
#' @export
lmbreak <- function(formula, data, pattern = NULL, start = NULL, range = NULL, 
                    se = FALSE,
                    control = NULL, trace = FALSE, digits = NULL){

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
    if(!is.null(pattern)){
        if(length(term.bp)>=3 && is.character(eval(term.bp[[3]]))){
            message("The pattern specified via the argument \'formula\' are ignored since argument \'pattern\' is not NULL. \n")
        }
    }else{
        if(length(term.bp)<3){
            stop("The number breakpoints and possible constrained should be specified in argument \'formula\'. \n",
                 "Something like: Y ~ bp(X, pattern = \"111\") for 2 breakpoints \n",
                 "                Y ~ bp(X, pattern = \"101\") for 2 breakpoints with a plateau. \n")
        }
        pattern <- eval(term.bp[[3]])
    }

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
    attr(data.fit, "min.bp") <- 0.5*(data.fit[[var.bp]][1]+data.fit[[var.bp]][2])
    attr(data.fit, "max.bp") <- 0.5*(data.fit[[var.bp]][NROW(data.fit)]+data.fit[[var.bp]][NROW(data.fit)-1])
    
    ## *** find range
    if(!is.null(range)){
        if(!is.numeric(range)){
            stop("Argument \'range\' should be a numeric vector. \n")
        }
        if(length(range)!=2){
            stop("Argument \'range\' should have length 2. \n")
        }
        if(any(is.na(range))){
            stop("Argument \'range\' should not contain NAs. \n")
        }
        if(range[1]>=range[2]){
            stop("The first element of argument \'range\' should be strictly superior to the second element. \n")
        }
    }else{
        range <- c(-Inf,Inf)
    }

    ## ** initialization

    ## *** breakpoint
    if(length(term.bp)==3 && is.numeric(eval(term.bp[[3]])) || length(term.bp)==4 || !is.null(start)){

        if(!is.null(start)){
            term.init <- start
            if(length(term.bp)==3 && is.numeric(eval(term.bp[[3]])) || length(term.bp)==4){
                message("The starting specified via the argument \'formula\' are ignored since argument \'start\' is not NULL. \n")
            }
        }else if(is.numeric(eval(term.bp[[3]]))){
            term.init <- cbind(as.double(eval(term.bp[[3]])))
        }else{
            term.init <- cbind(as.double(eval(term.bp[[4]])))
        }

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
            outInit.gam <- try(.init_gam(formula = terms.nobp, data = data.fit, var.bp = var.bp), silent = TRUE)
            if(!inherits(outInit.gam,"try-error")){
                breakpoint.prepare$gam <- outInit.gam[outInit.gam >= attr(data.fit, "min.bp") & outInit.gam <= attr(data.fit, "max.bp")]
                if(length(breakpoint.prepare$gam)<nMax.breakpoint && control$init.quantile<=0){
                    stop("Cannot initialize breakpoint position using gam. \n",
                         "Consider setting the argument \'init.quantile\' to TRUE or providing starting values, e.g. ~pb(X, \"111\", c(1,2)).")
                }
            }else if(control$init.quantile<=0){
                stop(outInit.gam)
            }
            
        }

        ## via quantiles
        ## if((control$optimizer == "Muggeo") && (control$init.quantile>0)){
        if(control$init.quantile>0){
            outInit.quantile <- .init_quantile(init = control$init.quantile, data = data.fit, var.bp = var.bp, n.breakpoint = nMax.breakpoint)
            breakpoint.prepare$quantile <- outInit.quantile[outInit.quantile >= attr(data.fit, "min.bp") & outInit.quantile <= attr(data.fit, "max.bp")]
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
                           init.type = factor(unlist(lapply(breakpoint.init, colnames)), unique(names(name.order))),
                           cv = NA, continuity = NA, regularity = NA, r2 = NA
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
        ## *** select pattern and initialization
        iPattern <- grid.fit[iGrid,"pattern"]
        iInit <- breakpoint.init[[iPattern]][,grid.fit[iGrid,"init"]]

        ## *** display
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

        ## *** optimize
        if((control$n.iter == 0) || control$optimizer=="Muggeo"){
            ls.fit[[iGrid]] <- optim.lmbreak_Muggeo(formula = patternUsVs[[iPattern]]$formula, formula.noVs = patternUsVs[[iPattern]]$formula.noVs, pattern = iPattern,
                                                    Us.label = patternUsVs[[iPattern]]$breakpoint.Us, Us.sign = patternUsVs[[iPattern]]$breakpoint.sign,
                                                    Vs.label = patternUsVs[[iPattern]]$Vs.label,
                                                    var.response = var.response, var.bp = var.bp, data = data.fit,
                                                    n.iter = control$n.iter, tol = control$tol, initialization = iInit, optimize.step = control$optimize.step,
                                                    trace = trace-1.999999, digits = digits)
        }else if(control$optimizer=="nlminb"){
            ls.fit[[iGrid]] <- optim.lmbreak_NLMINB(formula = patternUsVs[[iPattern]]$formula, formula.noVs = patternUsVs[[iPattern]]$formula.noVs, pattern = iPattern,
                                                    Us.label = patternUsVs[[iPattern]]$breakpoint.Us, Us.sign = patternUsVs[[iPattern]]$breakpoint.sign,
                                                    Vs.label = patternUsVs[[iPattern]]$Vs.label,
                                                    var.response = var.response, var.bp = var.bp, data = data.fit,
                                                    n.iter = control$n.iter, tol = control$tol, initialization = iInit, 
                                                    trace = trace-1.999999, digits = digits)
        }else if(control$optimizer %in% c("BFGS","L-BFGS-B")){
            ls.fit[[iGrid]] <- optim.lmbreak_BFGS(formula = patternUsVs[[iPattern]]$formula, formula.noVs = patternUsVs[[iPattern]]$formula.noVs, pattern = iPattern,
                                                  transform = control$optimizer == "BFGS",
                                                  Us.label = patternUsVs[[iPattern]]$breakpoint.Us, Us.sign = patternUsVs[[iPattern]]$breakpoint.sign,
                                                  Vs.label = patternUsVs[[iPattern]]$Vs.label,
                                                  var.response = var.response, var.bp = var.bp, data = data.fit,
                                                  n.iter = control$n.iter, tol = control$tol, initialization = iInit, 
                                                  trace = trace-1.999999, digits = digits)
        }

        ## *** check convergence
        iCheckCV <- .checkCV.lmbreak(object = ls.fit[[iGrid]], data = data.fit,
                                     formula = patternUsVs[[iPattern]]$formula, formula.noVs = patternUsVs[[iPattern]]$formula.noVs,
                                     Us.label = patternUsVs[[iPattern]]$Us.label, Vs.label = patternUsVs[[iPattern]]$Vs.label,
                                     var.bp = var.bp, tol =  control$tol, enforce.continuity = control$enforce.continuity, range = range)
        ls.fit[[iGrid]]$model <- iCheckCV$lm
        ls.fit[[iGrid]]$opt$continuity <- iCheckCV$continuity
        ls.fit[[iGrid]]$opt$regularity <- iCheckCV$regularity
        ls.fit[[iGrid]]$opt$r2 <- iCheckCV$r2
        grid.fit[iGrid,c("cv","continuity","regularity","r2")] <- unlist(ls.fit[[iGrid]]$opt[c("cv","continuity","regularity","r2")])

        ## *** summarize information about each phase
        ls.fit[[iGrid]]$phase <- .infoPhase.lmbreak(patternUsVs[[iPattern]], data = data.fit, breakpoint = ls.fit[[iGrid]]$breakpoint,
                                                    coef = ls.fit[[iGrid]]$model$coefficients, var.bp = var.bp, var.response = var.response)
        if(!is.null(attr(ls.fit[[iGrid]]$model,"continuity"))){
            attr(ls.fit[[iGrid]]$phase,"continuity") <- .infoPhase.lmbreak(patternUsVs[[iPattern]], data = data.fit, breakpoint = ls.fit[[iGrid]]$breakpoint,
                                                                           coef = attr(ls.fit[[iGrid]]$model,"continuity")$coefficients,var.bp = var.bp, var.response = var.response)
        }

        ## *** do not try other patterns or other initializations if a satisfactory solution is found
        if(grid.fit[iGrid,"check.cv"] && any(rowSums(grid.fit[1:iGrid,c("cv","continuity","regularity")])==3 & grid.fit[1:iGrid,"r2"]>control$minR2)){
            break
        }
    }
    if(trace>0){
        cat("\n")
    }

    ## ** select best fitting model
    if(iGrid==1){
        index.opt <- 1
    }else{
        index.opt <- which.max(2*(grid.fit[1:iGrid,"cv"] & grid.fit[1:iGrid,"r2"]>control$minR2) + grid.fit[1:iGrid,"regularity"]+grid.fit[1:iGrid,"r2"])
    }
    out <- ls.fit[[index.opt]]
    out$formula.pattern <- patternUsVs[[grid.fit[index.opt,"pattern"]]]$formula
    attr(out$formula.pattern,"noVs") <- patternUsVs[[grid.fit[index.opt,"pattern"]]]$formula.noVs

    attr(out$opt,"all") <- do.call(rbind,lapply(ls.fit[1:iGrid],"[[","opt"))
    ls.allbreakpoint <- lapply(ls.fit[1:iGrid], function(iM){iM$breakpoint$value})
    attr(out$breakpoint,"all") <- tapply(ls.allbreakpoint,factor(lengths(ls.allbreakpoint),unique(lengths(ls.allbreakpoint))),function(iLS){do.call(rbind, iLS)}, simplify = FALSE)

    ## ** standard error
    if(se && out$opt$cv && all(!is.na(out$breakpoint))){
        ## SANITY CHECK
        ## GS <- numDeriv::jacobian(.score.lmbreak, out$breakpoint$value, indiv = FALSE, transform = FALSE, formula = attr(out$formula.pattern,"noVs"), data = data, var.bp = var.bp, var.response = var.response, dX.skeleton = NULL, tol = options$tol[1])
        out$hessian <- .hessian.lmbreak(psi = out$breakpoint$value, transform = FALSE, formula = attr(out$formula.pattern,"noVs"),
                                        data = data, var.bp = var.bp, var.response = var.response, dX.skeleton = NULL, tol = options$tol[1])
        if(all(!is.na(out$hessian)) && (det(out$hessian)>options$tol[1])){
            out$breakpoint$se <- sqrt(diag(solve(out$hessian)))
        }else{
            out$breakpoint$se <- as.numeric(NA)
        }
    }else{
        out$breakpoint$se <- as.numeric(NA)
    }
    
    ## ** export
    if(out$opt$cv == FALSE){
        warning("The optimizer did not converge to a stable solution. \n")
    }else if(out$opt$regularity == FALSE){
        warning("The solution found by the optimizer has invalid breakpoint positions. \n")
    }else if(out$opt$continuity == FALSE){
        ## warning("The solution found by the optimizer is not continuous (non-0 Vs terms). \n")
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

## * helpers
## ** .pattern2UsVs (Us,Vs)
##' @description Re-express formual Y~bp(X, "11") into a formula with Vs and Us terms.
##' @noRd 
.pattern2UsVs <- function(pattern, terms){

    ## *** read input
    n.breakpoint <- length(pattern)-1

    ## *** generate Us terms
    Us.list <- vector(mode = "list", length = n.breakpoint+1)
    
    index.1 <- which(pattern!="0")
    addU.terms <- paste0("Us",index.1-1) ## Us0 is X
    sign.terms <- rep(1, length(addU.terms))
    if(any(pattern[-1]=="0")){
        index.0 <- setdiff(which(pattern=="0"),0)
        for(i0 in 1:length(index.0)){ ## i0 <- 3
            iIndex <- which(index.1<index.0[i0])
            if(length(iIndex)>0){
                addU.terms[iIndex] <- paste0("I(",addU.terms[iIndex],"-",paste0("Us",index.0[i0]-1),")")
                Us.list[[index.0[i0]]] <- data.frame(label = addU.terms[iIndex], sign = -1)
                index.1[iIndex] <- Inf
            }
        }        
    }
    Us.list[which(pattern!="0")] <- lapply(addU.terms, function(iLabel){data.frame(label = iLabel, sign = 1)})
    Us.terms <- stats::setNames(attr(stats::terms(stats::reformulate(addU.terms)), "term.labels"), addU.terms)

    ## *** generate Vs terms
    addV.terms <- paste0("Vs",1:n.breakpoint)
    Vs.terms <- attr(stats::terms(stats::reformulate(addV.terms)), "term.labels")

    ## *** breakpoint information
    breakpoint.pattern <- paste(pattern[1:n.breakpoint],pattern[2:(n.breakpoint+1)],c(pattern,1)[3:(n.breakpoint+2)], sep=".")
    breakpoint.Us <- Us.terms[breakpoint.pattern[1] %in% c("1.1.1","1.1.0","1.0.1") + cumsum(breakpoint.pattern != "1.0.1")]
    breakpoint.sign <- c(1,-1)[1+(breakpoint.pattern=="1.0.1")]

    ## *** formulas
    formula <- stats::update(terms, paste0(".~.+",paste(addU.terms, collapse="+"),"+",paste(addV.terms, collapse="+")))
    terms.formula <- stats::terms(formula)

    ## *** export
    return(list(formula = formula,
                formula.noVs = stats::update(terms, paste0(".~.+",paste(addU.terms, collapse="+"))),
                breakpoint.Us = breakpoint.Us,
                breakpoint.sign = breakpoint.sign,
                Vs.label = Vs.terms,
                Us.label = lapply(Us.list, function(iDf){if(!is.null(iDf)){iDf$label <- Us.terms[iDf$label]};return(iDf)}),
                n.param = length(attr(terms.formula,"term.labels"))
                ))

}

## ** .init_gam (initialization)
##' @description Find initialization values for breakpoints based on the change of sign of the first derivative of a fitted regression spline
##' @noRd 
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

## ** .init_quantile
##' @description Find initialization values for breakpoints based on the quantiles of the variable
##' @noRd 
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

## ** .checkCV.lmbreak
##' @description Check the properties of the estimated breakpoints
##' @noRd 
.checkCV.lmbreak <- function(object, formula, formula.noVs, Us.label, Vs.label, data, var.bp, tol, enforce.continuity, range){

    breakpoint <- object$breakpoint$value
    out <- list()
    
    ## *** test Vs terms
    if(is.na(tol[2]) || any(is.na(breakpoint)) || any(is.infinite(breakpoint))){
        out$continuity <- TRUE
        out$r2 <- NA
        out$lm <- NULL
    }else{
        data.UsVs <- model.frame.lmbreak(list(breakpoint.var = var.bp, breakpoint = breakpoint), newdata = data)
        out$lm <- stats::lm(formula, data = data.UsVs)
        param.UsVs <- stats::coef(out$lm)
        
        if(any(is.na(param.UsVs))){
            out$continuity <- FALSE
        }else{
            out$continuity <- all(abs(param.UsVs[Vs.label])<tol[2])
        }
        if(enforce.continuity){
            attr(out$lm,"continuity") <- stats::lm(formula.noVs, data = data.UsVs)
            out$r2 <- summary(attr(out$lm,"continuity"))$r.squared
        }else if(out$continuity){
            out$r2 <- summary(out$lm)$r.squared
        }else{
            out$r2 <- NA
        }
    }

    ## *** regularity: breakpoint non-consecutive (avoid staircase) and different from starting and ending values
    if(any(is.na(breakpoint))){
        out$regularity <- FALSE
    }else{

        if(min(breakpoint) < attr(data, "min.bp") || max(breakpoint) > attr(data, "max.bp")){
            out$regularity <- FALSE
        }else if(length(breakpoint)==1){
            out$regularity <- TRUE
        }else{
            ## NOTE: data is already sorted according to the breakpoint variable
            test.staircase <- sapply(breakpoint, function(iBreak){which(sign(data[[var.bp]]-iBreak)==1)[1]})
            out$regularity <- all(diff(test.staircase)>=1) & all(diff(breakpoint)>attr(data, "min.diff.bp"))
        }

        if(out$regularity && any(!is.infinite(range))){
            ## duration
            duration.bp <- diff(c(attr(data,"min.var.bp"),breakpoint,attr(data,"max.var.bp")))
            ## slope
            if(is.null(attr(out$lm,"continuity"))){
                wlm <- out$lm
            }else{
                wlm <- attr(out$lm,"continuity")
            }
            Us <- stats::coef(wlm)[stats::na.omit(c(ifelse("Us0" %in% names(stats::coef(wlm)),"Us0",NA),object$breakpoint$Us))]
            Us2slope <- stats::setNames(object$breakpoint$sign,object$breakpoint$Us)
            if(substring(object$opt$pattern,1,1)=="1"){
                Us2slope <- c(stats::setNames(1,names(Us)[1]),Us2slope)
                slope <- unname(cumsum(Us2slope * Us[names(Us2slope)]))
            }else{
                slope <- c(0,cumsum(Us2slope * Us[names(Us2slope)]))
            }
            y.bp <- cumsum(c(0,slope*duration.bp))
            if(any(y.bp<range[1] | y.bp>range[2])){
                out$regularity <- FALSE
            }
        }
    }

    ## *** export
    return(out)

}

## ** .infoPhase.lmbreak
##' @description Summarize information about the phase
##' @noRd
.infoPhase.lmbreak <- function(pattern, data, breakpoint, coef,
                               var.bp, var.response){

    ## *** prepare
    n.phase <- NROW(pattern$Us.label)
    bp.start <- range(data[[var.bp]][!is.na(data[[var.response]])], na.rm = TRUE)
    Vs <- stats::setNames(rep(0, n.phase-1), paste0("Vs",1:(n.phase-1)))
    if(length(intersect(names(Vs),names(coef)))>0){
        Vs[intersect(names(Vs),names(coef))] <- -coef[intersect(names(Vs),names(coef))]
    }

    ## *** compute
    vec.breakpoint <- c(bp.start[1],breakpoint$value,bp.start[2])
    vec.sign <- c(sapply(pattern$Us.label,function(iDf){ifelse(is.null(iDf),0,unique(iDf$sign))}),NA)
    ls.Us <- c(lapply(pattern$Us.label,"[[","label"), list(NULL))

    if(!is.null(coef)){
        vec.slope <- cumsum(c(sapply(pattern$Us.label, function(iDf){
            ifelse(is.null(iDf),0,sum(coef[iDf$label]*iDf$sign))
        }),NA))
        if("(Intercept)" %in% names(coef)){
            intercept <- coef["(Intercept)"]
        }else{
            intercept <- 0
        }
        vec.intercept <- intercept + cumsum(c(0,0,Vs)+c(vec.slope[c(1,1:n.phase)]*diff(c(0,vec.breakpoint))))
    }else{
        vec.intercept <- NA
        vec.slope <- NA
    }

    ## *** export
    out <- data.frame(breakpoint = vec.breakpoint,
                      intercept = vec.intercept,
                      slope = round(vec.slope,12), ## round to put to 0 slopes +/- same number that are not exactly 0 due to numerical approximation
                      Us = NA,
                      sign = vec.sign
                      )
    out$Us <- ls.Us
    return(out)
}

##----------------------------------------------------------------------
### lmbreak.R ends here
