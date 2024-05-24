#' ################################################
#' #      Project: JAN56830                       #
#' #      Survival Curve Estimation              #
#' #                                    
#' #The script imports data from digitalized
#' #Survival Curves and saves to a single data structure
#' #used by the master model.
#' #
#' #
#' #
#' #
#' ################################################
#Change Log
#
#Date           Programmer      Description
#20231001       OD              Created initial program
#                                   Saved as: surv_curves.R
#20231116       TJF             Modified to apply Step 1 HR from KOLs
#                               Added nonMM Mortality adjustments
#20231119       TJF             Adjusted for updated nonMM psuedo-data created on 11/19

# Load Needed Libraries ---------------------------------------------------
library("data.table")
library("survival")
library("flexsurv")

# Utility Functions ---------------------------------------------------
bind_ntimes <- function(dt, n=4){
  # Function that stacks a dataframe with itself 'n' times and creates a colum
  # giving the number a the table has been stacked.
  niter <- 1
  nrows <- nrow(dt)
  dt0 <- copy(dt)
  V4 <- rep(niter,nrows)
  while( niter < n){
    V4 <- c(V4, rep(niter +1, nrows))
    dt0 <- rbind(dt0,dt)
    niter <- niter +1
  }
  dt0[, V3:=V4]
  return(dt0)
}


# Initialize Data Objects --------------------------------
pfs.data <- vector(mode = 'list', length = 5)
os.data <-  vector(mode = 'list', length = 6)
nonMM.data <- vector(mode= 'list', length = 5)

# Read output from digitized curves  ----
## 1. P1-> R1_P2 ----
pfs.data[[1]] <- data.table(read.delim2("data_raw/KMdataIPD_Griffin_LOT1.txt", row.names = NULL,
                                         sep = "", header = T))
## 2, "P1 -> Death" ----
os.data[[1]] <- data.table(read.delim2("data_raw/KMdataIPD_OS_LOT1.txt", row.names = NULL,
                                         sep = "", header = T))

## 3. "R1_P2 -> R2_P3" ----
pfs.data[[2]] <- data.table(read.delim2("data_raw/KMdataIPD_Candor_LOT2.txt", row.names = NULL,
                                   sep = "", header = T))
## 4.  "R1_P2 -> Death" ----
os.data[[2]] <- data.table(read.delim2("data_raw/KMdataIPD_OS_LOT2.txt", row.names = NULL,
                                       sep = "", header = T))

## 5. "R2_P3 -> R3_P4" ----
pfs.data[[3]] <- data.table(read.delim2("data_raw/KMdataIPD_Eloquent_LOT3.txt", row.names = NULL,
                                   sep = "", header = T))
## 6. "R2_P3 -> Death" ----
os.data[[3]] <- data.table(read.delim2("data_raw/KMdataIPD_OS_LOT3.txt", row.names = NULL,
                                       sep = "", header = T))


## 7.  "R3_P4 -> R4_P5" ----
pfs.data[[4]] <- data.table(read.delim2("data_raw/KMdataIPD_Storm_LOT4.txt", row.names = NULL,
                                   sep = "", header = T))
## 8. "R3_P4 -> Death" ----
os.data[[4]] <- data.table(read.delim2("data_raw/KMdataIPD_OS_LOT4.txt", row.names = NULL,
                                       sep = "", header = T))

## 9. "R4_P5 -> R5" ----
pfs.data[[5]] <- data.table(read.delim2("data_raw/KMdataIPD_MajesTEC_LOT5.txt", row.names = NULL,
                                   sep = "", header = T))
## 10. "R4_P5 -> Death" ----
os.data[[5]] <- data.table(read.delim2("data_raw/KMdataIPD_OS_LOT5.txt", row.names = NULL,
                                       sep = "", header = T))
## 11 "P5 -> Death" ----
os.data[[6]] <- copy(os.data[[5]])

## Non MM-related Mortality estimation form psuedo IPD------------
nonMM.data[[1]] <- data.table(read.delim2("data_raw/KMdataIPD_NonMM_LOT1.txt", row.names = NULL,
                                        sep = "", header = T))
nonMM.data[[2]] <- data.table(read.delim2("data_raw/KMdataIPD_NonMM_LOT2.txt", row.names = NULL,
                                          sep = "", header = T))
nonMM.data[[3]] <- data.table(read.delim2("data_raw/KMdataIPD_NonMM_LOT3.txt", row.names = NULL,
                                          sep = "", header = T))
nonMM.data[[4]] <- data.table(read.delim2("data_raw/KMdataIPD_NonMM_LOT4.txt", row.names = NULL,
                                          sep = "", header = T))
nonMM.data[[5]] <- data.table(read.delim2("data_raw/KMdataIPD_NonMM_LOT5.txt", row.names = NULL,
                                          sep = "", header = T))


### Set time as numeric; event and strategy as integers---------------------
pfs.data <- lapply(pfs.data, function(x){bind_ntimes(x,4)})
pfs.data <- lapply(pfs.data, function(x){setNames(x,c("patient_id", "time", "event", "strategies"))})
pfs.data <- lapply(pfs.data, function(x){ x[, `:=`(time =as.numeric(time)/12, strategies = as.integer(strategies))] })
pfs.data <- setNames(pfs.data,paste0('line',1:length(pfs.data)))

os.data <- lapply(os.data, function(x){bind_ntimes(x,4)})
os.data <- lapply(os.data, function(x){setNames(x,c("patient_id", "time", "event", "strategies"))})
os.data <- lapply(os.data, function(x){ x[, `:=`(time =as.numeric(time)/12, strategies = as.integer(strategies))] })
os.data <- setNames(os.data,paste0('line',1:length(os.data)))

nonMM.data <- lapply(nonMM.data, function(x){bind_ntimes(x,4)})
nonMM.data <- lapply(nonMM.data, function(x){setNames(x,c("patient_id","time","event","strategies"))})
nonMM.data <- lapply(nonMM.data, function(x){x[,':=' (time=as.numeric(time)/12)]})

# Estimate Survival Functions (weibullPH function)-----------

ntransitions <- 11  #modified for adjusted curves

weilph.fit <- vector( mode = 'list',length = ntransitions)
nonMM.fit <-vector( mode='list', length=7) #output for mortality adjustments

# P1-> R1_P2
weilph.fit[[1]] <- flexsurvreg(Surv(time, event)~ I(strategies %in% c(1,2)) + I(strategies==3) ,dist = "weibullPH", 
                             data = pfs.data$line1) 
# "P1 -> Death"
weilph.fit[[2]] <- flexsurvreg(Surv(time, event)~ I(strategies %in% c(1,2)) + I(strategies==3) ,dist = "weibullPH", 
                               data = os.data$line1) 
# "P1 -> Death"
nonMM.fit[[1]] <- flexsurvreg(Surv(time, event)~ I(strategies %in% c(1,2)) + I(strategies==3) ,dist = "weibullPH", 
                               data = nonMM.data[[1]]) 


# "R1_P2 -> R2_P3"
weilph.fit[[3]] <- flexsurvreg(Surv(time, event)~ I(strategies ==1) + I(strategies == 2) + I(strategies==3) ,dist = "weibullPH", 
                               data = pfs.data$line2)
# "R1_P2 -> Death"
weilph.fit[[4]] <- flexsurvreg(Surv(time, event)~ I(strategies ==1) + I(strategies == 2) + I(strategies==3) ,dist = "weibullPH", 
                               data = os.data$line2) 
# "R1_P2 -> R2_P3"
nonMM.fit[[2]] <- flexsurvreg(Surv(time, event)~ I(strategies ==1) + I(strategies == 2) + I(strategies==3) ,dist = "weibullPH", 
                               data = nonMM.data[[2]])

# "R2_P3 -> R3_P4"
weilph.fit[[5]] <- flexsurvreg(Surv(time, event)~ I(strategies %in% c(1,3)) + I(strategies == 2) ,dist = "weibullPH", 
                               data = pfs.data$line3)
# "R2_P3 -> Death"
weilph.fit[[6]] <- flexsurvreg(Surv(time, event)~ I(strategies %in% c(1,3)) + I(strategies == 2) ,dist = "weibullPH", 
                               data = os.data$line3)
# "R2_P3 -> Death"
nonMM.fit[[3]] <- flexsurvreg(Surv(time, event)~ I(strategies %in% c(1,3)) + I(strategies == 2) ,dist = "weibullPH", 
                               data = nonMM.data[[3]])

# "R3_P4 -> R4_P5"
weilph.fit[[7]] <- flexsurvreg(Surv(time, event)~ I(strategies ==1) + I(strategies == 2) + I(strategies==3) ,dist = "weibullPH", 
                               data = pfs.data$line4)
# "R3_P4 -> Death"
weilph.fit[[8]] <- flexsurvreg(Surv(time, event)~ I(strategies ==1) + I(strategies == 2) + I(strategies==3) ,dist = "weibullPH", 
                               data = os.data$line4)
# "R3_P4 -> Death"
nonMM.fit[[4]] <- flexsurvreg(Surv(time, event)~ I(strategies ==1) + I(strategies == 2) + I(strategies==3) ,dist = "weibullPH", 
                               data = nonMM.data[[4]])
# "R4_P5 -> R5"
weilph.fit[[9]] <- flexsurvreg(Surv(time, event)~ I(strategies ==1) + I(strategies == 2) + I(strategies==3) ,dist = "weibullPH", 
                               data = pfs.data$line5)
# "R4_P5 -> Death"
weilph.fit[[10]] <- flexsurvreg(Surv(time, event)~ I(strategies ==1) + I(strategies == 2) + I(strategies==3) ,dist = "weibullPH", 
                               data = os.data$line5)
# "R4_P5 -> Death"
nonMM.fit[[5]] <- flexsurvreg(Surv(time, event)~ I(strategies ==1) + I(strategies == 2) + I(strategies==3) ,dist = "weibullPH", 
                                data = nonMM.data[[5]])
# "R5 -> Death"
weilph.fit[[11]] <- flexsurvreg(Surv(time, event)~ 1,dist = "weibullPH", 
                                data = os.data$line6[strategies == 4,])
            # RETIRED TO SIMPLY MODEL MORTALITY CURVE OPTIONS
            # for (i in 1:length(nonMM.data)) {
            #     nonMM.fit[[i]]  <- flexsurvreg(Surv(time, event)~ 1,dist = "weibullph",data =nonMM.data[[i]]) 
            # }

# "OS adjustments" ----------------
## These are based upon the percentage of MM-attributable deaths occurring up to median time of progression
##  for each line of therapy
##  this is an iterative function using median PFS times
### Key Elements: 
###   mPFS: vector of median progression free survival times using estimated weibull curves in weibph.fit
###   mOS:  vector of median overall survival times using estimated weibull curves in weibph.fit (based upon Braunlin curves)
###   nonMM.fit: (saved for use in overall model) 3 element list with:
###             1) estimated weibullph survival curve of non MM mortality (US Life Tables 2022 and Seer Incidence 2022)
###             2) estimated %of deaths attributable to MM up to median PFS time by LOT
###             3) estimated %of deaths attributable to MM up to median OS time by LOT
###       NOTE: elements 2 and 3 will be used to adjusted user input HR when applied to survival curves

mPFS<-  c(qweibullPH(.5,weilph.fit[[1]]$res[1],weilph.fit[[1]]$res[2],lower.tail=TRUE),
          qweibullPH(.5,weilph.fit[[3]]$res[1],weilph.fit[[3]]$res[2],lower.tail=TRUE),
          qweibullPH(.5,weilph.fit[[5]]$res[1],weilph.fit[[5]]$res[2],lower.tail=TRUE),
          qweibullPH(.5,weilph.fit[[7]]$res[1],weilph.fit[[7]]$res[2],lower.tail=TRUE),
          qweibullPH(.5,weilph.fit[[9]]$res[1],weilph.fit[[9]]$res[2],lower.tail=TRUE)
          )

mOS<-  c(qweibullPH(.5,weilph.fit[[2]]$res[1],weilph.fit[[2]]$res[2],lower.tail=TRUE),
          qweibullPH(.5,weilph.fit[[4]]$res[1],weilph.fit[[4]]$res[2],lower.tail=TRUE),
          qweibullPH(.5,weilph.fit[[5]]$res[1],weilph.fit[[6]]$res[2],lower.tail=TRUE),
          qweibullPH(.5,weilph.fit[[8]]$res[1],weilph.fit[[8]]$res[2],lower.tail=TRUE),
          qweibullPH(.5,weilph.fit[[10]]$res[1],weilph.fit[[10]]$res[2],lower.tail=TRUE)
)

# Pct of Deaths up to median PFS
nonMM.fit[[6]] <- c(1-((pweibullPH(mPFS[1],weilph.fit[[2]]$res[1],weilph.fit[[2]]$res[2],lower.tail=TRUE))^-1)*
                       (pweibullPH(mPFS[1],nonMM.fit[[1]]$res[1],nonMM.fit[[1]]$res[2],lower.tail=TRUE)
                        ),  #LOT-1
                    1-((pweibullPH(mPFS[2],weilph.fit[[4]]$res[1],weilph.fit[[4]]$res[2],lower.tail=TRUE))^-1)*
                      (pweibullPH(mPFS[2],nonMM.fit[[2]]$res[1],nonMM.fit[[2]]$res[2],lower.tail=TRUE)
                      ),  #LOT-2
                    1-((pweibullPH(mPFS[3],weilph.fit[[6]]$res[1],weilph.fit[[6]]$res[2],lower.tail=TRUE))^-1)*
                      (pweibullPH(mPFS[3],nonMM.fit[[3]]$res[1],nonMM.fit[[3]]$res[2],lower.tail=TRUE)
                      ),  #LOT-3
                    1-((pweibullPH(mPFS[4],weilph.fit[[8]]$res[1],weilph.fit[[8]]$res[2],lower.tail=TRUE))^-1)*
                      (pweibullPH(mPFS[4],nonMM.fit[[4]]$res[1],nonMM.fit[[4]]$res[2],lower.tail=TRUE)
                      ),  #LOT-4
                    1-((pweibullPH(mPFS[5],weilph.fit[[10]]$res[1],weilph.fit[[10]]$res[2],lower.tail=TRUE))^-1)*
                      (pweibullPH(mPFS[5],nonMM.fit[[5]]$res[1],nonMM.fit[[5]]$res[2],lower.tail=TRUE)
                      )  #LOT-5
)
  
# Pct of Deaths up to median OS
nonMM.fit[[7]] <- c(1-((pweibullPH(mOS[1],weilph.fit[[2]]$res[1],weilph.fit[[2]]$res[2],lower.tail=TRUE))^-1)*
                      (pweibullPH(mOS[1],nonMM.fit[[1]]$res[1],nonMM.fit[[1]]$res[2],lower.tail=TRUE)
                      ),  #LOT-1
                    1-((pweibullPH(mOS[2],weilph.fit[[4]]$res[1],weilph.fit[[4]]$res[2],lower.tail=TRUE))^-1)*
                      (pweibullPH(mOS[2],nonMM.fit[[2]]$res[1],nonMM.fit[[2]]$res[2],lower.tail=TRUE)
                      ),  #LOT-2
                    1-((pweibullPH(mOS[3],weilph.fit[[6]]$res[1],weilph.fit[[6]]$res[2],lower.tail=TRUE))^-1)*
                      (pweibullPH(mOS[3],nonMM.fit[[3]]$res[1],nonMM.fit[[3]]$res[2],lower.tail=TRUE)
                      ),  #LOT-3
                    1-((pweibullPH(mOS[4],weilph.fit[[8]]$res[1],weilph.fit[[8]]$res[2],lower.tail=TRUE))^-1)*
                      (pweibullPH(mOS[4],nonMM.fit[[4]]$res[1],nonMM.fit[[4]]$res[2],lower.tail=TRUE)
                      ),  #LOT-4
                    1-((pweibullPH(mOS[5],weilph.fit[[10]]$res[1],weilph.fit[[10]]$res[2],lower.tail=TRUE))^-1)*
                      (pweibullPH(mOS[5],nonMM.fit[[5]]$res[1],nonMM.fit[[5]]$res[2],lower.tail=TRUE)
                      )  #LOT-5
)

save(weilph.fit, file = "data/weilph.data.rda", compress = "bzip2")
save(nonMM.fit, file = "data/nonMM.data.rda", compress = "bzip2")

