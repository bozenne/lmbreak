### test-manual-psilo2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr  8 2024 (20:05) 
## Version: 
## Last-Updated: apr 12 2024 (15:13) 
##           By: Brice Ozenne
##     Update #: 47
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

    e.SERT <- mlmbreak(score ~ 0 + bp(time.num, c("01010","1010","101","011","11")), cluster = "PatientID", data = dtL.SERT)
    expect_equal(logLik(e.SERT), -1019.248, tol = 1e-3)
    e.SERT
    plot(e.SERT, labeller = label_value)
    model.tables(e.SERT, format = "list")

    ## ee.SERT <- mlmbreak(score ~ 0 + bp(time.num, c("01010","1010","101","110","011","11")), cluster = "PatientID", data = dtL.SERT, control = list(optimize.step = 0.5, n.iter = 20))
    
    ## ** Analyse only LPS individuals
    dtL.LPS <- dtL[grep("^LPS",dtL$PatientID),,drop=FALSE]
    ggTraj %+% dtL.LPS

    e.LPS <- mlmbreak(score ~ 0 + bp(time.num, c("01010","10101","1010","0101","101","011","110","11")), cluster = "PatientID", data = dtL.LPS)
    expect_equal(logLik(e.LPS), -1741.223, tol = 1e-3)
    summary(e.LPS)
    coef(e.LPS, "R2")
    plot(e.LPS, labeller = label_value)
    model.tables(e.LPS, format = "list")

    ## manual solution
    ## e012.LPS <- lmbreak(score ~ 0 + bp(time.num, c("01010")), data = dtL.LPS[dtL.LPS$PatientID == "LPS-012",], control = list(optimize.step = 0.5))
    e012.LPS <- lmbreak(score ~ 0 + bp(time.num, c("01010"), c(30,100,200,450)), data = dtL.LPS[dtL.LPS$PatientID == "LPS-012",], control = list(optimize.step = 0.5))
    plot(e012.LPS)
    
    ## e.GS <- segmented(lm(score ~ time.num, data = dtL.LPS[dtL.LPS$PatientID == "LPS-012",]), npsi = 4)
    

    ## ** Analyse only L-LPS individuals
    dtL.LLPS <- dtL[grep("^L-LPS",dtL$PatientID),,drop=FALSE]
    ggTraj %+% dtL.LLPS

    e.LLPS <- mlmbreak(score ~ 0 + bp(time.num, c("01010","10101","1010","0101","101","011","110","11")), cluster = "PatientID", data = dtL.LLPS)
    ## ee.LLPS <- mlmbreak(score ~ 0 + bp(time.num, c("01010","10101","1010","0101","101","011","110","11")), cluster = "PatientID", data = dtL.LLPS, control = list(optimize.step = 0.5))
    expect_equal(logLik(e.LLPS), -1694.284, tol = 1e-3)
    summary(e.LLPS)
    coef(e.LLPS, "R2")
    plot(e.LLPS, labeller = label_value)
    model.tables(e.LPS, format = "list")

    ## manual initialization
    dtL2.LLPS <- dtL.LLPS[dtL.LLPS$PatientID=="L-LPS-002",]
    e2.LLPS <- lmbreak(score ~ 0 + bp(time.num, c("110")), data = dtL2.LLPS, control = list(optimize.step = 0.5))
    e2.LLPS$opt
    plot(e2.LLPS)

    dtL5.LLPS <- dtL.LLPS[dtL.LLPS$PatientID=="L-LPS-005",]
    e5.LLPS <- lmbreak(score ~ 0 + bp(time.num, c("110")), data = dtL5.LLPS, control = list(optimize.step = 0.5))
    plot(e5.LLPS)

    dtL6.LLPS <- dtL.LLPS[dtL.LLPS$PatientID=="L-LPS-006",]
    e6.LLPS <- lmbreak(score ~ 0 + bp(time.num, c("0110")), data = dtL6.LLPS, control = list(optimize.step = 0.5))
    plot(e6.LLPS)

    dtL7.LLPS <- dtL.LLPS[dtL.LLPS$PatientID=="L-LPS-007",]
    e19.LLPS <- lmbreak(score ~ 0 + bp(time.num, c("110")), data = dtL19.LLPS, optimize.step = 0.5)
    plot(e19.LLPS)

    dtL24.LLPS <- dtL.LLPS[dtL.LLPS$PatientID=="L-LPS-024",]
    e24.LLPS <- lmbreak(score ~ 0 + bp(time.num, c("0110")), data = dtL24.LLPS, optimize.step = 0.5)
    plot(e24.LLPS)

    ## ** Analyse only LPM individuals
    dtL.LPM <- dtL[grep("^LPM",dtL$PatientID),,drop=FALSE]
    ggTraj %+% dtL.LPM

    e.LPM <- mlmbreak(score ~ 0 + bp(time.num, c("1010","0110","110","101","11")), cluster = "PatientID", data = dtL.LPM)
    expect_equal(logLik(e.LPM), -2008.148, tol = 1e-3)
    summary(e.LPM)
    coef(e.LPM, "R2")
    plot(e.LPM, labeller = label_value, scales = "free", ylim = c(0,100))
    model.tables(e.LPM, format = "list")

    dtL12.LPM <- dtL.LPM[dtL.LPM$PatientID=="LPM-012",]
    e12.LPM <- lmbreak(score ~ 0 + bp(time.num, c("110"), c(90,200)), data = dtL12.LPM, optimize.step = 1, trace = 5)
    plot(e12.LPM)
}


## * Dataset 2
if(FALSE){

    ## ** load data
    dtW <- as.data.table(read_excel("source/28_09_2018_Intensitetsratings_18_fp.xlsx"))
    dtW.red <- dtW[which(rowSums(!is.na(dtW.red))>0),.SD,.SDcols = c("CIMBI ID","0 minutes/adm.",seq(20,440,by=20))]
    names(dtW.red)[1:2] <- c("id","0")
    dtL <- melt(dtW.red, id.vars = c("id"), variable.name = "time", value.name = "score")
    dtL$time.num <- as.numeric(as.character(dtL$time))

    ## ** display
    ggTraj <- ggplot(dtL, aes(x = time.num, y = score, group = id))
    ggTraj <- ggTraj + geom_line() + geom_point()
    ggTraj <- ggTraj + facet_wrap(~id)
    ggTraj

    ## ** trajectory model
    e.mbp <- mlmbreak(score ~ 0 + bp(time.num, c("1010","0110","110","101","11")), cluster = "id", data = dtL)
    plot(e.mbp)
    model.tables(e.mbp, format = "list")
}


##----------------------------------------------------------------------
### test-manual-psilo2.R ends here
