# Main script for running analyses
library(dplyr)
library(tidyr)
library(RBesT)
library(SAMprior)
library(patchwork)
library(Rcpp)
library(forcats)
library(RcppArmadillo)



# Main analyses -----------------------------------------------------------

# Power and error rate for normal data (< 5 minutes)
source("scripts/error_rates_normal.R")

# Power and error rate for binary data (< 5 minutes)
source("scripts/error_rates_binary.R")

# Power and error rate for survival data (< 30 minutes)
source("scripts/error_rates_survival.R")

# Simulation comparison for normal data (< 3 hours)
source("scripts/sims_normal.R")

# Simulation comparison for binary data (< 3 hours)
source("scripts/sims_binary.R")

# Application
source("scripts/application.R")

# Sensitivity analyses ----------------------------------------------------

# Simulation comparison for normal data and small sample size (< 5 minutes)
source("scripts/sims_normal_smalln.R")

# Simulation comparison for normal data (< 3 hours)
source("scripts/sims_normal_recalibrated.R")

# Simulation comparison for binary data (< 3 hours)
source("scripts/sims_binary_recalibrated.R")

# Power and error rate for survival data (< 45 minutes)
source("scripts/error_rates_survival_weibull.R")