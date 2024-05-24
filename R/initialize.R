#' Project JAN56830

#' This code should be run to initialize the R model implementation. It installs
#' all required R packages, and creates input baseline data tables if they are missing.
 

# Install packages used in project ----

app.packages <- c("hesim","data.table","tidyr", "flexsurv", "ggplot2", 'readxl', 
                  "scales","survminer", "bookdown", "knitr", "kableExtra",
                  "WriteXLS", "writexl")
to.install.packages <- app.packages[!(app.packages %in% installed.packages()[,"Package"])]
if(length(to.install.packages) > 0) install.packages(to.install.packages)

if(length(ls())>0) rm(list=ls())

# Create baseline input tables in case they are missing from the data folder
library(data.table)
library(dplyr)
source('R/aux_functions.R')
source('R/TXcostsUpdate.R')

load_create_fun()

rm(list=ls())

#' ########### END  ################