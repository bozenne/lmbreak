### test-manual-psilo2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr  8 2024 (20:05) 
## Version: 
## Last-Updated: Apr 20 2024 (18:52) 
##           By: Brice Ozenne
##     Update #: 61
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

if(FALSE){
    library(readxl)
    library(data.table)
    library(ggplot2)

    library(lmbreak)
}

## * Dataset 1
if(FALSE){

    ## ** load data
    dtW <- as.data.table(read_excel("source/SDI_Psilo_22.11.2023.xlsx"))
    dtL <- melt(dtW, id.vars = c("Study","PatientID"), variable.name = "time", value.name = "score")
    dtL$time.num <- as.numeric(as.character(dtL$time))

    ## ** display
    ggTraj <- ggplot(dtL, aes(x = time.num, y = score, group = PatientID, color = Study))
    ggTraj <- ggTraj + geom_line() + geom_point()
    ggTraj <- ggTraj + facet_wrap(~PatientID)
    ggTraj

    ## ** Analyse only SERT individuals
    dtL.SERT <- dtL[grep("SERT",dtL$PatientID),,drop=FALSE]
    ggTraj %+% dtL.SERT

    e.SERT <- mlmbreak(score ~ 0 + bp(time.num), pattern = c("01010","1010","101","011","11"), cluster = "PatientID", data = dtL.SERT, control = list(optimizer = "Muggeo"))
    expect_equal(logLik(e.SERT), -955.0167, tol = 1e-3)
    eDESCENTUN.SERT <- mlmbreak(score ~ 0 + bp(time.num), pattern = c("01010","1010","101","011","11"), cluster = "PatientID", data = dtL.SERT, control = list(optimizer = "BFGS"))
    expect_equal(logLik(eDESCENTUN.SERT), -842.4136, tol = 1e-3)
    
    plot(e.SERT, eDESCENTUN.SERT, labeller = label_value)

    ## manual solution
    ## e019.SERT <- lmbreak(score ~ 0 + bp(time.num), pattern = c("01010","1010","101","011","11"), data = dtL.SERT[dtL.SERT$PatientID == "SERT-Psilo-019",], control = list(optimizer = "BFGS"))
    ## plot(e019.SERT)
    ## plot(as.lmbreak(eDESCENTUN.SERT, cluster = "SERT-Psilo-019"))
 
    ## ** Analyse only LPS individuals
    dtL.LPS <- dtL[grep("^LPS",dtL$PatientID),,drop=FALSE]
    ggTraj %+% dtL.LPS

    e.LPS <- mlmbreak(score ~ 0 + bp(time.num), pattern = c("01010","10101","1010","0101","101","011","110","11"), cluster = "PatientID", data = dtL.LPS, control = list(optimizer = "Muggeo"))
    expect_equal(logLik(e.LPS), -1688.759, tol = 1e-3)
    eDESCENTUN.LPS <- mlmbreak(score ~ 0 + bp(time.num), pattern = c("01010","10101","1010","0101","101","011","110","11"), cluster = "PatientID", data = dtL.LPS, control = list(optimizer = "BFGS"))
    expect_equal(logLik(eDESCENTUN.LPS), -1499.123, tol = 1e-3)
 
    plot(e.LPS, eDESCENTUN.LPS, labeller = label_value)
    model.tables(eDESCENTUN.LPS, format = "list")

    ## manual solution
    ## e006.LPS <- lmbreak(score ~ 0 + bp(time.num, c("01010","10101","1010","0101","101","011","110","11")), data = dtL.LPS[dtL.LPS$PatientID == "LPS-006",], control = list(optimizer = "BFGS"))
    ## plot(e006.LPS)


    ## ** Analyse only L-LPS individuals
    dtL.LLPS <- dtL[grep("^L-LPS",dtL$PatientID),,drop=FALSE]
    ggTraj %+% dtL.LLPS

    e.LLPS <- mlmbreak(score ~ 0 + bp(time.num), pattern = c("01010","10101","1010","0101","101","011","110","11"), cluster = "PatientID", data = dtL.LLPS, control = list(optimizer = "Muggeo"))
    expect_equal(logLik(e.LLPS), -1684.527, tol = 1e-3)
    eDESCENTUN.LLPS <- mlmbreak(score ~ 0 + bp(time.num), pattern = c("01010","10101","1010","0101","101","011","110","11"), cluster = "PatientID", data = dtL.LLPS, control = list(optimizer = "BFGS"))
    expect_equal(logLik(eDESCENTUN.LLPS), -1429.695, tol = 1e-3)
   
    plot(e.LLPS, eDESCENTUN.LLPS, labeller = label_value)
   
    coef(eDESCENTUN.LLPS, "R2")
    model.tables(eDESCENTUN.LLPS, format = "list")

    ## ** Analyse only LPM individuals
    dtL.LPM <- dtL[grep("^LPM",dtL$PatientID),,drop=FALSE]
    ggTraj %+% dtL.LPM

    e.LPM <- mlmbreak(score ~ 0 + bp(time.num), pattern = c("01010","10101","1010","0101","101","011","110","11"), cluster = "PatientID", data = dtL.LPM, control = list(optimizer = "Muggeo"))
    expect_equal(logLik(e.LPM), -2013.906, tol = 1e-3)
    eDESCENTUN.LPM <- mlmbreak(score ~ 0 + bp(time.num), pattern = c("01010","10101","1010","0101","101","011","110","11"), cluster = "PatientID", data = dtL.LPM, control = list(optimizer = "BFGS"))
    expect_equal(logLik(eDESCENTUN.LPM), -1673.175, tol = 1e-3)

    plot(e.LPM, eDESCENTUN.LPM, labeller = label_value)
}


## * Dataset 2
if(FALSE){

    ## ** load data
    dtW <- as.data.table(read_excel("source/28_09_2018_Intensitetsratings_18_fp.xlsx"))
    dtW.red <- dtW[which(rowSums(!is.na(dtW))>0),.SD,.SDcols = c("CIMBI ID","0 minutes/adm.",seq(20,440,by=20))]
    names(dtW.red)[1:2] <- c("id","0")
    dtL <- suppressWarnings(melt(dtW.red, id.vars = c("id"), variable.name = "time", value.name = "score"))
    dtL$time.num <- as.numeric(as.character(dtL$time))

    ## ** display
    ggTraj <- ggplot(dtL, aes(x = time.num, y = score, group = id))
    ggTraj <- ggTraj + geom_line() + geom_point()
    ggTraj <- ggTraj + facet_wrap(~id)
    ggTraj

    ## ** trajectory model
    e.FP <- mlmbreak(score ~ 0 + bp(time.num), pattern = c("01010","10101","1010","0101","101","011","110","11"), cluster = "id", data = dtL, control = list(optimizer = "Muggeo"))
    expect_equal(logLik(e.FP), -239.0657, tol = 1e-3)
    eDESCENTUN.FP <- mlmbreak(score ~ 0 + bp(time.num), pattern = c("01010","10101","1010","0101","101","011","110","11"), cluster = "id", data = dtL, control = list(optimizer = "BFGS"))
    expect_equal(logLik(eDESCENTUN.FP), -188.1214, tol = 1e-3)

    plot(e.FP, eDESCENTUN.FP, labeller = label_value)
}


##----------------------------------------------------------------------
### test-manual-psilo2.R ends here
