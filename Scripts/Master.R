#Master script for executing the analyses for the article "Estimating the time to statistical cure"
#As the simulations utilize the mclapply function, the scripts have to be executed on a UNIX platform.

#Figure and table directories
fig.out <- "."
tab.out <- "."

#Set generated data directory
data.out <- "."

#Load absolute paths to use for project. Comment out if absolute paths are set in the above lines
source("Scripts/directories.R")


#Create directories if not existing
dir.create(fig.out, showWarnings = F, recursive = T)
dir.create(tab.out, showWarnings = F, recursive = T)


#Load required packages
library(survival)
library(relsurv)
library(ggplot2)
library(numDeriv)
library(xtable)
library(rootSolve)
library(rstpm2)
library(parallel)
library(gridExtra)

#Load own package (can be found on github (https://github.com/LasseHjort/cuRe))
library(cuRe)

#Table format
tab.format <- "%.3f"

#Set year
ayear <- 365.24

#Figure format variable
pdf <- TRUE

#Original figure settings
mai_par <- par("mai")
mfrow_par <- par("mfrow")

#Color palette
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# #Load clinical data
# source("Scripts/LoadData2.R", encoding = "utf-8")
# 
# #Run analysis
# source("Scripts/RunAnalysis.R")
# source("Scripts/Analyze_ColonCancer.R")

#Run simulations
source("Scripts/Setup_simulations.R")

