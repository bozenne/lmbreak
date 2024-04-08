### test-manual-psilo2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr  8 2024 (20:05) 
## Version: 
## Last-Updated: apr  8 2024 (20:16) 
##           By: Brice Ozenne
##     Update #: 1
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:
library(readxl)
dtW <- as.data.table(read_excel("SDI_Psilo_22.11.2023.xlsx"))
dtL <- melt(dtW, id.vars = c("Study","PatientID"), variable.name = "time", value.name = "score")
dtL$time.num <- as.numeric(as.character(dtL$time))

library(ggplot2)
ggTraj <- ggplot(dtL, aes(x = time.num, y = score, group = PatientID, color = Study))
ggTraj <- ggTraj + geom_line() + geom_point()
ggTraj <- ggTraj + facet_wrap(~PatientID)
ggTraj

SERT-PSILO-008
SERT-PSILO-025

dtL.008red <- dtL[dtL$PatientID=="SERT-Psilo-008",][!is.na(score)]
e.lmbreak_008 <- lmbreak(score ~ bp(time.num, "1010"), data = dtL.008red)
plot(e.lmbreak_008)

dtL.017red <- dtL[dtL$PatientID=="SERT-Psilo-017",][!is.na(score)]
e.lmbreak_017 <- lmbreak(score ~ 0 + bp(time.num, "1010"), data = dtL.017red, nQuantile.init = 7, trace = TRUE)
plot(e.lmbreak_017)
e.lmbreak_017$model

dtL.025red <- dtL[dtL$PatientID=="SERT-Psilo-025",][!is.na(score)]
e.lmbreak_025 <- lmbreak(score ~ bp(time.num, "1010"), data = dtL.025red)
plot(e.lmbreak_025)
##----------------------------------------------------------------------
### test-manual-psilo2.R ends here
