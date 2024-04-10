## * load packages
library(testthat)
library(lmbreak)

library(ggplot2)
library(segmented)

## * run tests
## setwd("~/Documents/GitHub/lmbreak/tests/")
test_check("lmbreak")
