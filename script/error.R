#!/bin/Rscript

absDiff <- function(arg1,arg2) {
  return(sqrt((arg1-arg2)^2))
}

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop("Input data file name not provided\n", call. = FALSE)
}

in_original <- args[1]
in_acr <- args[2]

fdtdValOriginal <- read.table(in_original, col.names=c("xpos","ypos","physicalQuantity"))
fdtdValApprox <- read.table(in_acr, col.names=c("xpos","ypos","physicalQuantity"))

if (!all.equal(fdtdValOriginal[,c("xpos","ypos")],fdtdValApprox[,c("xpos","ypos")])) {
  stop("The data grid of the two files does not match\n", call. = FALSE)
}

maxdens <- max(fdtdValOriginal$physicalQuantity)
mindens <- min(fdtdValOriginal$physicalQuantity)
range  <- absDiff(maxdens,mindens)

summaryData <- data.frame(unclass(summary(absDiff(fdtdValApprox$physicalQuantity, fdtdValOriginal$physicalQuantity) / range * 100)))
names(summaryData) <- ""
print(summaryData, digits = 3, zero.print = "0")
cat(sprintf("Range %.3f\n", range))
