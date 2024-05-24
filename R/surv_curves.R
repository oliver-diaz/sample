# Create survival models from output of digitized curves ----

# Function that stacks a data frame with itself 'n' times 

library("data.table")
library("survival")
library("flexsurv")

bind_ntimes <- function(dt, n=4){
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

pfs.data <- vector(mode = 'list', length = 5)
os.data <-  vector(mode = 'list', length = 6)

# KMdataIPD_OS_LOTX.txt death progressin data comes from Braunlin study

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


### Set time as numeric; event and strategy as integers

pfs.data <- lapply(pfs.data, function(x){bind_ntimes(x,4)})
pfs.data <- lapply(pfs.data, function(x){setNames(x,c("patient_id", "time", "event", "strategies"))})
pfs.data <- lapply(pfs.data, function(x){ x[, `:=`(time =as.numeric(time)/12, strategies = as.integer(strategies))] })
pfs.data <- setNames(pfs.data,paste0('line',1:length(pfs.data)))

os.data <- lapply(os.data, function(x){bind_ntimes(x,4)})
os.data <- lapply(os.data, function(x){setNames(x,c("patient_id", "time", "event", "strategies"))})
os.data <- lapply(os.data, function(x){ x[, `:=`(time =as.numeric(time)/12, strategies = as.integer(strategies))] })
os.data <- setNames(os.data,paste0('line',1:length(os.data)))


ntransitions <- 11

weilph.fit <- vector( mode = 'list',length = ntransitions)

# P1-> R1_P2
weilph.fit[[1]] <- flexsurvreg(Surv(time, event)~ I(strategies %in% c(1,2)) + I(strategies==3) ,dist = "weibullPH", 
                             data = pfs.data$line1) 
# "P1 -> Death"
weilph.fit[[2]] <- flexsurvreg(Surv(time, event)~ I(strategies %in% c(1,2)) + I(strategies==3) ,dist = "weibullPH", 
                               data = os.data$line1) 


# "R1_P2 -> R2_P3"
weilph.fit[[3]] <- flexsurvreg(Surv(time, event)~ I(strategies ==1) + I(strategies == 2) + I(strategies==3) ,dist = "weibullPH", 
                               data = pfs.data$line2)
# "R1_P2 -> Death"
weilph.fit[[4]] <- flexsurvreg(Surv(time, event)~ I(strategies ==1) + I(strategies == 2) + I(strategies==3) ,dist = "weibullPH", 
                               data = os.data$line2) 


# "R2_P3 -> R3_P4"
weilph.fit[[5]] <- flexsurvreg(Surv(time, event)~ I(strategies %in% c(1,3)) + I(strategies == 2) ,dist = "weibullPH", 
                               data = pfs.data$line3)
# "R2_P3 -> Death"
weilph.fit[[6]] <- flexsurvreg(Surv(time, event)~ I(strategies %in% c(1,3)) + I(strategies == 2) ,dist = "weibullPH", 
                               data = os.data$line3)

# "R3_P4 -> R4_P5"
weilph.fit[[7]] <- flexsurvreg(Surv(time, event)~ I(strategies ==1) + I(strategies == 2) + I(strategies==3) ,dist = "weibullPH", 
                               data = pfs.data$line4)
# "R3_P4 -> Death"
weilph.fit[[8]] <- flexsurvreg(Surv(time, event)~ I(strategies ==1) + I(strategies == 2) + I(strategies==3) ,dist = "weibullPH", 
                               data = os.data$line4)

# "R4_P5 -> R5"
weilph.fit[[9]] <- flexsurvreg(Surv(time, event)~ I(strategies ==1) + I(strategies == 2) + I(strategies==3) ,dist = "weibullPH", 
                               data = pfs.data$line5)
# "R4_P5 -> Death"
weilph.fit[[10]] <- flexsurvreg(Surv(time, event)~ I(strategies ==1) + I(strategies == 2) + I(strategies==3) ,dist = "weibullPH", 
                               data = os.data$line5)

# "R5 -> Death"
weilph.fit[[11]] <- flexsurvreg(Surv(time, event)~ 1,dist = "weibullPH", 
                                data = os.data$line6[strategies == 4,])

save(weilph.fit, file = "data/weilph.data.rda", compress = "bzip2")


#####
# Check the "logical" order of survival in terms of line of treatment
# PFS
# plot(weilph.fit[[1]], col = 'black')
# lines(weilph.fit[[3]], col = 'darkgreen')
# lines(weilph.fit[[5]], col = 'red')
# lines(weilph.fit[[7]], col = 'magenta')
# lines(weilph.fit[[9]], col = 'blue') # perhaps some adjustment to lower it down

# OS ckeck
# plot(weilph.fit[[2]], col = 'black')
# lines(weilph.fit[[4]], col = 'darkgreen')
# lines(weilph.fit[[6]], col = 'red')
# lines(weilph.fit[[8]], col = 'magenta')
# lines(weilph.fit[[10]], col = 'blue')


# library(magrittr)
# 
# n_samples <- 1000
# 
# transmod_coef_def <- define_rng({
#   for(i in 1:ntransitions){
#     assign(paste0("t",i),multi_normal_rng(mu = weilph.fit[[i]]$res.t[,1], Sigma = vcov(weilph.fit[[i]]))) 
#   }
#   list(
#     t1_shape = vec_to_dt(t1$shape),
#     t1_scale = t1[, -1,] |> (\(x){colnames(x)=c("cons",  "notst3", "st3"); return(x)})() ,
#     t2_shape = vec_to_dt(t2$shape),
#     t2_scale = t2[, -1,] |> (\(x){colnames(x)=c("cons", "notst3", "st3"); return(x)})() , 
#     t3_shape = vec_to_dt(t3$shape),
#     t3_scale = t3[, -1,] |> (\(x){colnames(x)=c("cons", "st1", "st2", "st3"); return(x)})() ,
#     t4_shape = vec_to_dt(t4$shape),
#     t4_scale = t4[, -1,] |> (\(x){colnames(x)=c("cons", "st1", "st2", "st3"); return(x)})() ,
#     t5_shape = vec_to_dt(t5$shape),
#     t5_scale = t5[, -1,] |> (\(x){colnames(x)=c("cons", "notst2", "st2"); return(x)})() ,
#     t6_shape = vec_to_dt(t6$shape),
#     t6_scale = t6[, -1,] |> (\(x){colnames(x)=c("cons", "notst2", "st2"); return(x)})() ,
#     t7_shape = vec_to_dt(t7$shape),
#     t7_scale = t7[, -1,] |> (\(x){colnames(x)=c("cons", "st1", "st2", "st3"); return(x)})() ,
#     t8_shape = vec_to_dt(t8$shape),
#     t8_scale = t8[, -1,] |> (\(x){colnames(x)=c("cons", "st1", "st2", "st3"); return(x)})() ,
#     t9_shape = vec_to_dt(t9$shape),
#     t9_scale = t9[, -1,] |> (\(x){colnames(x)=c("cons", "st1", "st2", "st3"); return(x)})() ,
#     t10_shape = vec_to_dt(t10$shape),
#     t10_scale = t10[, -1,] |> (\(x){colnames(x)=c("cons", "st1", "st2", "st3"); return(x)})() ,
#     t11_shape = vec_to_dt(t11$shape),
#     t11_scale = t11[, -1,] |> (\(x){colnames(x)=c("cons"); return(x)})()
#   )
# }, n = n_samples)
# 
# transmod_coef <- eval_rng(transmod_coef_def)
# 
# # Tables with HR ranges (per transition where needed) ----
# HR.data <- list(
#   t1.hr = data.frame(low = c(notst3 = .8, st3 = .9), up= c(notst3 = .91, st3 = 1.01)),
#   t2.hr = data.frame(low = c(notst3 = .8, st3 = .9), up= c(st12 = .91, st3 = 1.01)),
#   t3.hr = data.frame(low = c(st1 =.9, st2= .65, st3= .97), up=c(st1=.96, st2=.9, st3= 1.01)),
#   t4.hr = data.frame(low = c(st1 =.9, st2= .65, st3= .97), up=c(st1=.96, st2=.9, st3= 1.01)),
#   t5.hr = data.frame(low = c(notst2 =.9, st2= .65), up=c(notst2=.96, st2=.9)),
#   t6.hr = data.frame(low = c(notst2 =.7, st2= .65), up=c(notst2=.9, st2=.9)),
#   t7.hr = data.frame(low = c(st1 =.9, st2= .65, st3= .97), up=c(st1=.96, st2=.9, st3= 1.01)),
#   t8.hr = data.frame(low = c(st1 =.9, st2= .65, st3= .97), up=c(st1=.96, st2=.9, st3= 1.01)),
#   t9.hr = data.frame(low = c(st1 =.9, st2= .65, st3= .97), up=c(st1=.96, st2=.9, st3= 1.01)),
#   t10.hr = data.frame(low = c(st1 =.9, st2= .65, st3= .97), up=c(st1=.96, st2=.9, st3= 1.01))
# )
# 
# set.seed(111)
# 
# # Adjust scale according to (sampled) hazard rations from KOLs ----
# ## P1-> R1_P2 ----
# transmod_coef$t1_scale <- transmod_coef$t1_scale +
# data.table(const = 0, 
#            notst3 = log(runif(n_samples,HR.data$t1.hr$low[1], HR.data$t1.hr$up[1])),
#            st3 =  log(runif(n_samples,HR.data$t1.hr$low[2], HR.data$t1.hr$up[2]))
# )
# 
# ## P1-> Death----
# transmod_coef$t2_scale <- transmod_coef$t2_scale +
#   data.table(const = 0, 
#              notst3 = log(runif(n_samples,HR.data$t2.hr$low[1], HR.data$t2.hr$up[1])),
#              st3 =  log(runif(n_samples,HR.data$t2.hr$low[2], HR.data$t2.hr$up[2]))
#   )
# 
# ##  "R1_P2 -> R2_P3" ----
# transmod_coef$t3_scale <- transmod_coef$t3_scale +
#   data.table(const = 0, 
#               st1 = log(runif(n_samples,HR.data$t3.hr$low[1], HR.data$t3.hr$up[1])),
#               st2 =  log(runif(n_samples,HR.data$t3.hr$low[2], HR.data$t3.hr$up[2])),
#               st3 =  log(runif(n_samples,HR.data$t3.hr$low[3], HR.data$t3.hr$up[3]))
#   )
# ##  "R1_P2 -> Death" ----
# transmod_coef$t4_scale <- transmod_coef$t4_scale +
#   data.table(const = 0, 
#              st1 = log(runif(n_samples,HR.data$t4.hr$low[1], HR.data$t4.hr$up[1])),
#              st2 =  log(runif(n_samples,HR.data$t4.hr$low[2], HR.data$t4.hr$up[2])),
#              st3 =  log(runif(n_samples,HR.data$t4.hr$low[3], HR.data$t4.hr$up[3]))
#   )
# 
# # "R2_P3 -> R3_P4" ----
# transmod_coef$t5_scale <- transmod_coef$t5_scale +
#   data.table(const = 0, 
#              notst2 = log(runif(n_samples,HR.data$t5.hr$low[1], HR.data$t5.hr$up[1])),
#              st2 =  log(runif(n_samples,HR.data$t5.hr$low[2], HR.data$t5.hr$up[2]))
#   )
# # "R2_P3 -> Death" ----
# transmod_coef$t6_scale <- transmod_coef$t6_scale +
#   data.table(const = 0, 
#              notst2 = log(runif(n_samples,HR.data$t6.hr$low[1], HR.data$t6.hr$up[1])),
#              st2 =  log(runif(n_samples,HR.data$t6.hr$low[2], HR.data$t6.hr$up[2]))
#   )
# 
# # "R3_P4 -> R4_P5" ----
# transmod_coef$t7_scale <- transmod_coef$t7_scale +
#   data.table(const = 0, 
#              st1 = log(runif(n_samples,HR.data$t7.hr$low[1], HR.data$t7.hr$up[1])),
#              st2 =  log(runif(n_samples,HR.data$t7.hr$low[2], HR.data$t7.hr$up[2])),
#              st3 =  log(runif(n_samples,HR.data$t7.hr$low[3], HR.data$t7.hr$up[3]))
#   )
# # "R3_P4 -> Death" ----
# transmod_coef$t8_scale <- transmod_coef$t8_scale +
#   data.table(const = 0, 
#              st1 = log(runif(n_samples,HR.data$t8.hr$low[1], HR.data$t8.hr$up[1])),
#              st2 =  log(runif(n_samples,HR.data$t8.hr$low[2], HR.data$t8.hr$up[2])),
#              st3 =  log(runif(n_samples,HR.data$t8.hr$low[3], HR.data$t8.hr$up[3]))
#   )
# 
# # "R4_P5 -> R5"
# transmod_coef$t9_scale <- transmod_coef$t9_scale +
#   data.table(const = 0, 
#              st1 = log(runif(n_samples,HR.data$t9.hr$low[1], HR.data$t9.hr$up[1])),
#              st2 =  log(runif(n_samples,HR.data$t9.hr$low[2], HR.data$t9.hr$up[2])),
#              st3 =  log(runif(n_samples,HR.data$t9.hr$low[3], HR.data$t9.hr$up[3]))
#   )
# # "R4_P5 -> Death"
# transmod_coef$t10_scale <- transmod_coef$t10_scale +
#   data.table(const = 0, 
#              st1 = log(runif(n_samples,HR.data$t10.hr$low[1], HR.data$t10.hr$up[1])),
#              st2 =  log(runif(n_samples,HR.data$t10.hr$low[2], HR.data$t10.hr$up[2])),
#              st3 =  log(runif(n_samples,HR.data$t10.hr$low[3], HR.data$t10.hr$up[3]))
#   )