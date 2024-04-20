### optimizer.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Apr 20 2024 (15:24) 
## Version: 
## Last-Updated: Apr 20 2024 (18:41) 
##           By: Brice Ozenne
##     Update #: 34
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * optim.lmbreak_Muggeo
optim.lmbreak_Muggeo <- function(formula, formula.noVs, pattern, Us.label, Us.sign, Vs.label,
                                 var.response, var.bp, data,
                                 n.iter, tol, optimize.step, initialization,
                                 trace, digits){

    ## *** prepare
    n.breakpoint <- length(initialization)
    bp.range <- c(attr(data, "min.var.bp"),attr(data, "max.var.bp"))

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

    ## *** export
    df.opt <- data.frame(initialization = NA,
                         optimizer = "Muggeo",
                         n.iter = iIter,
                         cv = cv,
                         continuity = NA,
                         regularity = NA,
                         tol = NA,
                         pattern = pattern,
                         diff = NA,
                         R2 = NA
                         )
    df.opt$tol <- list(tol)
    df.opt$initialization <- list(initialization)
    df.opt$diff <- list(iDiff)
    return(list(model = iE.lm,
                breakpoint = data.frame(value = iBreakpoint, Us = Us.label, Vs = Vs.label, sign = Us.sign),
                opt = df.opt))
}

## * optim.lmbreak_NLMINB (grid search fitter)
optim.lmbreak_NLMINB <- function(formula, formula.noVs, pattern, Us.label, Us.sign, Vs.label,
                                 var.response, var.bp, data,
                                 n.iter, tol, initialization, 
                                 trace, digits){

    ## *** prepare
    n.breakpoint <- length(Us.label)
    bp.range <- c(attr(data,"min.var.bp"),attr(data,"max.var.bp"))
    bp.min.diff <- attr(data,"min.diff.bp")
    Id <- diag(1,NROW(data))
    data$Us0 <- data[[var.bp]]

    ## *** grid search
    res.optim <- stats::nlminb(start = initialization,
                               objective = .RSS.lmbreak, indiv = FALSE, transform = FALSE, formula = formula.noVs, data = data, var.bp = var.bp, var.response = var.response, 
                               lower = rep(bp.range[1] + bp.min.diff, n.breakpoint),
                               upper = rep(bp.range[2]- bp.min.diff, n.breakpoint),
                               control = list(abs.tol = tol[1], eval.max = 10*n.iter, iter.max = 7.5*n.iter))
    iBreakpoint <- res.optim$par
    cv <- res.optim$convergence==0

    ## *** export
    df.opt <- data.frame(initialization = NA,
                         optimizer = "nlminb",
                         n.iter = res.optim$iteration,
                         cv = cv,
                         continuity = NA,
                         regularity = NA,
                         tol = NA,
                         pattern = pattern,
                         diff = NA,
                         R2 = NA
                         )
    df.opt$tol <- list(tol)
    df.opt$initialization <- list(initialization)
    df.opt$diff <- list(NULL)
    return(list(model = NULL,
                breakpoint = data.frame(value = iBreakpoint, Us = Us.label, Vs = Vs.label, sign = Us.sign),
                opt = df.opt))
}


## * optim.lmbreak_BFGS (BFGS with subgradient)
optim.lmbreak_BFGS <- function(formula, formula.noVs, pattern, Us.label, Us.sign, Vs.label, transform,
                               var.response, var.bp, data,
                               n.iter, tol, initialization, 
                               trace, digits){

    
    ## *** prepare
    n.breakpoint <- length(Us.label)
    bp.range <- c(attr(data,"min.bp"),attr(data,"max.bp"))
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
        ## .RSS.lmbreak(initialization, indiv = FALSE, transform = FALSE, formula = formula.noVs, data = data, var.bp = var.bp, var.response = var.response)
        ## .RSS.lmbreak(initialization.trans, indiv = FALSE, transform = transform, formula = formula.noVs, data = data, var.bp = var.bp, var.response = var.response)
        ## numDeriv::jacobian(.RSS.lmbreak, indiv = FALSE, initialization.trans, transform = transform, formula = formula.noVs, data = data, var.bp = var.bp, var.response = var.response)
        ## .score.lmbreak(initialization.trans, indiv = FALSE, transform = transform, formula = formula.noVs, data = data, var.bp = var.bp, var.response = var.response, dX.skeleton = dX.skeleton, tol = tol[1])
        res.optim <- stats::optim(par = initialization.trans, method = "BFGS",
                                  fn = function(psi){.RSS.lmbreak(psi, indiv = FALSE, transform = transform, formula = formula.noVs, data = data, var.bp = var.bp, var.response = var.response)},
                                  gr = function(psi){.score.lmbreak(psi, indiv = FALSE, transform = transform, formula = formula.noVs, data = data, var.bp = var.bp, dX.skeleton = dX.skeleton, var.response = var.response, tol = tol[1])},
                                  control = list(maxit = 5*n.iter))
        iBreakpoint <- backtransformPsi(res.optim$par, min = bp.range[1], max = bp.range[2], mindiff = bp.min.diff, jacobian = FALSE)
    }else{       
         
        ## SANITY CHECK
        ## .RSS.lmbreak(initialization, indiv = FALSE, transform = FALSE, formula = formula.noVs, data = data, var.bp = var.bp, var.response = var.response)
        ## .RSS.lmbreak(rep(bp.range[1], n.breakpoint) + bp.min.diff/2, indiv = FALSE, transform = FALSE, formula = formula.noVs, data = data, var.bp = var.bp, var.response = var.response)
        ## .RSS.lmbreak(rep(bp.range[2], n.breakpoint) - bp.min.diff/2, indiv = FALSE, transform = FALSE, formula = formula.noVs, data = data, var.bp = var.bp, var.response = var.response)
        ## numDeriv::jacobian(.RSS.lmbreak, indiv = FALSE, initialization, transform = FALSE, formula = formula.noVs, data = data, var.bp = var.bp, var.response = var.response)
        ## .score.lmbreak(initialization, indiv = FALSE, transform = FALSE, formula = formula.noVs, data = data, var.bp = var.bp, var.response = var.response, dX.skeleton = dX.skeleton, tol = tol[1])
        res.optim <- stats::optim(par = initialization, method = "L-BFGS-B",
                                  fn = function(psi){.RSS.lmbreak(psi, indiv = FALSE, transform = FALSE, formula = formula.noVs, data = data, var.bp = var.bp, var.response = var.response)},
                                  gr = function(psi){.score.lmbreak(psi, indiv = FALSE, transform = FALSE, formula = formula.noVs, data = data, var.bp = var.bp, dX.skeleton = dX.skeleton, var.response = var.response, tol = tol[1])},
                                  lower = rep(bp.range[1], n.breakpoint) + bp.min.diff/2,
                                  upper = rep(bp.range[2], n.breakpoint) - bp.min.diff/2,
                                  control = list(maxit = 5*n.iter))
        iBreakpoint <- res.optim$par
    }
    cv <- res.optim$convergence==0

    ## *** export
    df.opt <- data.frame(initialization = NA,
                         optimizer = ifelse(transform,"BFGS","L-BFGS-B"),
                         n.iter = res.optim$counts[1],
                         cv = cv,
                         continuity = NA,
                         regularity = NA,
                         tol = NA,
                         pattern = pattern,
                         diff = NA,
                         R2 = NA
                         )
    df.opt$tol <- list(tol)
    df.opt$initialization <- list(initialization)
    df.opt$diff <- list(NULL)
    return(list(model = NULL,
                breakpoint = data.frame(value = iBreakpoint, Us = Us.label, Vs = Vs.label, sign = Us.sign),
                opt = df.opt))
}



##----------------------------------------------------------------------
### optimizer.R ends here
