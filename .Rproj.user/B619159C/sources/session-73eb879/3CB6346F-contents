#' ##############################################################
#' #  Author: Oliver Diaz                                       #
#' #  Precision HEOR                                            #
#' #  Economic Model of Alternative Treatment Sequences in RRMM #
#' #                                                            #
#' #  Model Flow and Information                                #
#' #                                                            #
#' #  Change Log                                                #
#' #                                                            #
#' #  Code Review                                               #
#' #                                                            #
#' ##############################################################

library("hesim")
library("data.table")
library("survival")
library("flexsurv")
library("ggplot2")
library("survminer")
library("dplyr")

source('R/aux_functions.R')
source('R/eval_tx_cost.R') # function that evaluates tx-costs based on disease progression simulation

# Check and/or create inputs ----
## If data files with inputs (medical costs, utilities, treatment costs, survival data)
## have not been created from data files, it creates them
load_create_fun()

## Load predetermined model inputs ----
load("data/weilph.data.rda") # loads reference survival data
load("data/utility_data.rda") # loads utility inputs
load("data/medcost_data.rda") # loads medical cost inputs
load("data/cost_tx_seq.rda") # loads tables with treatments costs by state and strategy
disc <- 0.03 # Annual discount rate

# Population ----
n_patients <- 100 #?Initial Input ?#
n_samples <- 200  #?Initial Input ?#

### no population covariates are defined for this project
patients <- data.table(patient_id=1:n_patients)

# States ----
## 1. P1: On 1L therapy ----
## 2. R1_P2: relapsed on 1L; moved to 2L therapy ----
## 3. R2_P3: relapsed on 2L; moved to 3L therapy ----
## 4. R3_P4: relapsed on 3L; moved to 4L therapy ----
## 5. R4_P5: relapsed on 4L; moved to 5L therapy ----
## 6. R5: relapsed on 5L therapy ----
## 7. Death ----

tmat <- rbind( c(NA,1,NA,NA, NA, NA, 2),
               c(NA, NA, 3, NA, NA, NA , 4),
               c(NA, NA, NA, 5, NA, NA , 6),
               c(NA, NA, NA, NA, 7, NA , 8),
               c(NA, NA, NA, NA, NA, 9 , 10),
               c(NA, NA, NA, NA, NA, NA , 11),
               c(NA, NA, NA, NA, NA, NA , NA)
               )

colnames(tmat) <- rownames(tmat) <- c("P1", "R1_P2", "R2_P3", "R3_P4", "R4_P5","R5", "Death")
print(tmat)

# Strategies ----
## Seq1: (DVRd + ASCT + DR )  > DKd         > EloPd   > Sd   > Cilta-cel ----
## Seq2: (DVRd + ASCT + DR )  > Cilta-cel   > TalDara > Tec  > SKd ----
## Seq3: (DVrd + DRd)         > Tec         > EloPd   > Tal  > Kd ----

strategies <- data.table(strategy_id = c(1, 2, 3),
                         strategy_name = c("Sequence 1", "Sequence 2", "Sequence 3"))

no_death_states <- ncol(tmat)-1
states <- data.table(state_id = 1:no_death_states,
                     state_name = rownames(tmat)[1:no_death_states]) # Non-death health states
hesim_dat <- hesim_data( strategies = strategies,
                         patients = patients,
                         states = states)

# Transitions ----
## Labels for the plots and summary tables.
labs <- get_labels(hesim_dat)
labs$transition_id <- c("P1-> R1_P2" = 1, 
                        "P1 -> Death" = 2,
                        "R1_P2 -> R2_P3" = 3,
                        "R1_P2 -> Death" = 4,
                        "R2_P3 -> R3_P4" = 5,
                        "R2_P3 -> Death" = 6,
                        "R3_P4 -> R4_P5" = 7,
                        "R3_P4 -> Death" = 8,
                        "R4_P5 -> R5" = 9,
                        "R4_P5 -> Death" = 10,
                        "R5 -> Death" = 11)
print(labs)

# Data ----
## Progression parameters ----
## Survival models with Weibull (PH) distribution are created from digitized
## survival curves in the code surv_curves.R
##  NOTE: These are pre-analyzed and saved in the weilph.data.rda dataset

library(magrittr)  #?move earlier with others?#

set.seed(3455)  #?Input ?#

transmod_coef_def <- define_rng({
  for(i in 1:length(labs$transition_id)){
    assign(paste0("t",i),multi_normal_rng(mu = weilph.fit[[i]]$res.t[,1], Sigma = vcov(weilph.fit[[i]]))) 
  }
  list(
    t1_shape = vec_to_dt(t1$shape),
    t1_scale = t1[, -1,] |> (\(x){colnames(x)=c("cons",  "notst3", "st3"); return(x)})() ,
    t2_shape = vec_to_dt(t2$shape),
    t2_scale = t2[, -1,] |> (\(x){colnames(x)=c("cons", "notst3", "st3"); return(x)})() , 
    t3_shape = vec_to_dt(t3$shape),
    t3_scale = t3[, -1,] |> (\(x){colnames(x)=c("cons", "st1", "st2", "st3"); return(x)})() ,
    t4_shape = vec_to_dt(t4$shape),
    t4_scale = t4[, -1,] |> (\(x){colnames(x)=c("cons", "st1", "st2", "st3"); return(x)})() ,
    t5_shape = vec_to_dt(t5$shape),
    t5_scale = t5[, -1,] |> (\(x){colnames(x)=c("cons", "notst2", "st2"); return(x)})() ,
    t6_shape = vec_to_dt(t6$shape),
    t6_scale = t6[, -1,] |> (\(x){colnames(x)=c("cons", "notst2", "st2"); return(x)})() ,
    t7_shape = vec_to_dt(t7$shape),
    t7_scale = t7[, -1,] |> (\(x){colnames(x)=c("cons", "st1", "st2", "st3"); return(x)})() ,
    t8_shape = vec_to_dt(t8$shape),
    t8_scale = t8[, -1,] |> (\(x){colnames(x)=c("cons", "st1", "st2", "st3"); return(x)})() ,
    t9_shape = vec_to_dt(t9$shape),
    t9_scale = t9[, -1,] |> (\(x){colnames(x)=c("cons", "st1", "st2", "st3"); return(x)})() ,
    t10_shape = vec_to_dt(t10$shape),
    t10_scale = t10[, -1,] |> (\(x){colnames(x)=c("cons", "st1", "st2", "st3"); return(x)})() ,
    t11_shape = vec_to_dt(t11$shape),
    t11_scale = t11[, -1,] |> (\(x){colnames(x)=c("cons"); return(x)})()
  )
}, n = n_samples)

transmod_coef <- eval_rng(transmod_coef_def)

## Tables with HR ranges from KOLs (per transition where needed) ----
HR.data <- list(
  t1.hr = data.frame(low = c(notst3 = .8, st3 = .9), up= c(notst3 = .91, st3 = 1)),
  t2.hr = data.frame(low = c(notst3 = .9, st3 = .9), up= c(notst3 = 1, st3 = 1)),
  t3.hr = data.frame(low = c(st1 =.9, st2= .60, st3= .97), up=c(st1=.96, st2=.65, st3= 1)),
  t4.hr = data.frame(low = c(st1 =.95, st2= .8, st3= .97), up=c(st1=.96, st2=.85, st3= 1)),
  t5.hr = data.frame(low = c(notst2 =.9, st2= .8), up=c(notst2=.96, st2=.9)),
  t6.hr = data.frame(low = c(notst2 =.91, st2= .85), up=c(notst2=.96, st2=1)),
  t7.hr = data.frame(low = c(st1 =.9, st2= .85, st3= .97), up=c(st1=.96, st2=.9, st3= 1)),
  t8.hr = data.frame(low = c(st1 =.9, st2= .9, st3= .97), up=c(st1=1, st2=1, st3= 1)),
  t9.hr = data.frame(low = c(st1 =1, st2= 1.02, st3= 1.02), up=c(st1=1.02, st2=1.2, st3= 1.2)),
  t10.hr = data.frame(low = c(st1 =1, st2=1.02, st3= 1.02), up=c(st1=1.02, st2=1.2, st3= 1.2)),
  t11.hr = data.frame(low = c(st1 =2.8, st2=2.8, st3= 2.8), up=c(st1=3.1, st2=3.1, st3= 3.1))
)

### Survival adjustments with HRs----
### P1-> R1_P2
transmod_coef$t1_scale[,`:=`(
             notst3 = log(runif(n_samples,HR.data$t1.hr$low[1], HR.data$t1.hr$up[1])),
             st3 =  log(runif(n_samples,HR.data$t1.hr$low[2], HR.data$t1.hr$up[2]))
  )]

### P1-> Death
transmod_coef$t2_scale[,`:=`( 
             notst3 = log(runif(n_samples,HR.data$t2.hr$low[1], HR.data$t2.hr$up[1])),
             st3 =  log(runif(n_samples,HR.data$t2.hr$low[2], HR.data$t2.hr$up[2]))
  )]

###  R1_P2 -> R2_P3
transmod_coef$t3_scale[,`:=`( 
             st1 = log(runif(n_samples,HR.data$t3.hr$low[1], HR.data$t3.hr$up[1])),
             st2 =  log(runif(n_samples,HR.data$t3.hr$low[2], HR.data$t3.hr$up[2])),
             st3 =  log(runif(n_samples,HR.data$t3.hr$low[3], HR.data$t3.hr$up[3]))
  )]
###  R1_P2 -> Death
transmod_coef$t4_scale[,`:=`( 
             st1 = log(runif(n_samples,HR.data$t4.hr$low[1], HR.data$t4.hr$up[1])),
             st2 =  log(runif(n_samples,HR.data$t4.hr$low[2], HR.data$t4.hr$up[2])),
             st3 =  log(runif(n_samples,HR.data$t4.hr$low[3], HR.data$t4.hr$up[3]))
  )]

### R2_P3 -> R3_P4
transmod_coef$t5_scale[,`:=`(
             notst2 = log(runif(n_samples,HR.data$t5.hr$low[1], HR.data$t5.hr$up[1])),
             st2 =  log(runif(n_samples,HR.data$t5.hr$low[2], HR.data$t5.hr$up[2]))
  )]
### R2_P3 -> Death
transmod_coef$t6_scale[ , `:=`(
             notst2 = log(runif(n_samples,HR.data$t6.hr$low[1], HR.data$t6.hr$up[1])),
             st2 =  log(runif(n_samples,HR.data$t6.hr$low[2], HR.data$t6.hr$up[2]))
  )]

### R3_P4 -> R4_P5
transmod_coef$t7_scale[, `:=`( 
             st1 = log(runif(n_samples,HR.data$t7.hr$low[1], HR.data$t7.hr$up[1])),
             st2 =  log(runif(n_samples,HR.data$t7.hr$low[2], HR.data$t7.hr$up[2])),
             st3 =  log(runif(n_samples,HR.data$t7.hr$low[3], HR.data$t7.hr$up[3]))
  )]
### R3_P4 -> Death
transmod_coef$t8_scale[, `:=`( 
             st1 = log(runif(n_samples,HR.data$t8.hr$low[1], HR.data$t8.hr$up[1])),
             st2 =  log(runif(n_samples,HR.data$t8.hr$low[2], HR.data$t8.hr$up[2])),
             st3 =  log(runif(n_samples,HR.data$t8.hr$low[3], HR.data$t8.hr$up[3]))
  )]

### R4_P5 -> R5
transmod_coef$t9_scale[, `:=`(  
             st1 = log(runif(n_samples,HR.data$t9.hr$low[1], HR.data$t9.hr$up[1])),
             st2 =  log(runif(n_samples,HR.data$t9.hr$low[2], HR.data$t9.hr$up[2])),
             st3 =  log(runif(n_samples,HR.data$t9.hr$low[3], HR.data$t9.hr$up[3]))
  )]
### R4_P5 -> Death
transmod_coef$t10_scale[, `:=`(  
             st1 = log(runif(n_samples,HR.data$t10.hr$low[1], HR.data$t10.hr$up[1])),
             st2 =  log(runif(n_samples,HR.data$t10.hr$low[2], HR.data$t10.hr$up[2])),
             st3 =  log(runif(n_samples,HR.data$t10.hr$low[3], HR.data$t10.hr$up[3]))
  )]

### R5 -> Death
transmod_coef$t11_scale[, `:=`(  
  st1 = log(runif(n_samples,HR.data$t11.hr$low[1], HR.data$t11.hr$up[1])),
  st2 =  log(runif(n_samples,HR.data$t11.hr$low[2], HR.data$t11.hr$up[2])),
  st3 =  log(runif(n_samples,HR.data$t11.hr$low[3], HR.data$t11.hr$up[3]))
)]

summary(transmod_coef)

rm(weilph.fit)

## Define survival parameters for input into hesim model---- 
transmod_params <- params_surv_list(
  # 1. P1-> R1_P2
  params_surv(coefs = list(shape = transmod_coef$t1_shape,
                           scale = transmod_coef$t1_scale),
              dist = "weibullPH"),
  # 2. P1 -> Death
  params_surv(coefs = list(shape = transmod_coef$t2_shape,
                           scale = transmod_coef$t2_scale), 
              dist = "weibullPH"),
  # 3. R1_P2 -> R2_P3
  params_surv(coefs = list(shape = transmod_coef$t3_shape,
                           scale = transmod_coef$t3_scale), 
              dist = "weibullPH"),
  # 4. R1_P2 -> Death
  params_surv(coefs = list(shape = transmod_coef$t4_shape,
                           scale = transmod_coef$t4_scale), 
              dist = "weibullPH"),
  # 5. R2_P3 -> R3_P4
  params_surv(coefs = list(shape = transmod_coef$t5_shape,
                           scale = transmod_coef$t5_scale), 
              dist = "weibullPH"),
  # 6. R2_P3 -> Death
  params_surv(coefs = list(shape = transmod_coef$t6_shape,
                           scale = transmod_coef$t6_scale), 
              dist = "weibullPH"),
  # 7. R3_P4 -> R4_P5
  params_surv(coefs = list(shape = transmod_coef$t7_shape,
                           scale = transmod_coef$t7_scale), 
              dist = "weibullPH"),
  # 8 R3_P4 -> Death
  params_surv(coefs = list(shape = transmod_coef$t8_shape,
                           scale = transmod_coef$t8_scale),
              dist = "weibullPH"),
  # 9. R4_P5 -> R5
  params_surv(coefs = list(shape = transmod_coef$t9_shape,
                           scale = transmod_coef$t9_scale),
              dist = "weibullPH"),
  # 10. R4_P5 -> Death
  params_surv(coefs = list(shape = transmod_coef$t10_shape,
                           scale = transmod_coef$t10_scale),
              dist = "weibullPH"),
  # 11. R5 -> Death
  params_surv(coefs = list(shape = transmod_coef$t11_shape,
                           scale = transmod_coef$t11_scale),
              dist = "weibullPH")
)


## Create an 'input dataset' with all factors to be considered in hesim simulation ----
## NOTE: hesim cannot directly use factors and specific indicators are needed
transmod_data <- expand(hesim_dat, by = c("strategies", "patients"))
head(transmod_data)

transmod_data[, cons := 1]  ## const created by transmod
transmod_data[, st1:= ifelse(strategy_name == "Sequence 1",1,0) ]
transmod_data[, st2:= ifelse(strategy_name == "Sequence 2",1,0) ]
transmod_data[, st3:= ifelse(strategy_name == "Sequence 3",1,0) ]
transmod_data[, notst3:= ifelse(strategy_name != "Sequence 3",1,0)]
transmod_data[, notst2:= ifelse(strategy_name != "Sequence 2",1,0)]

## Transition model ----
## Creates heSim Transition Model Object
transmod <- create_IndivCtstmTrans(transmod_params, 
                                   input_data = transmod_data,
                                   trans_mat = tmat,
                                   clock = "reset")

## Utility by LOT----
# example of who values can be entered manually.
# utility_tbl <- stateval_tbl(data.table(expand.grid(state_id = states$state_id, strategy_id = strategies$strategy_id)[,2:1],
#                                        mean = rep(c(.85,.82, .80,.77, .73,.6),3),
#                                        se = rep(c(.107,0.113, 0.129, 0.096, 0.103, 0.128),3)
#                                        ),
#                             dist = "beta")

utility_tbl <- stateval_tbl(utility_data,
                            dist = "beta")
head(utility_tbl)

#### Create heSim State Object ----
utilitymod <- create_StateVals(utility_tbl, n = n_samples, hesim_data = hesim_dat)

## Non-Treatment Medical Costs by LOT ----

### Inpatient costs ----
# example of who values can be entered manually.
# medcostip_tbl <- stateval_tbl(data.table(expand.grid(state_id = states$state_id, strategy_id = strategies$strategy_id)[,2:1],
#                                        mean = rep(c(1000, 1100, 1100,1200,1200,2000),3),
#                                        se = rep(rep(c(10.95445, 14.14214),3),3)),
#                             dist = "gamma")
medcostip_tbl <- stateval_tbl(medcost_data$inpatient,
                            dist = "gamma")
# head(medcostip_tbl)

### Outpatient costs ----
medcostop_tbl <- stateval_tbl(medcost_data$outpatient,
                              dist = "gamma")
# head(medcostop_tbl)

### Other outpatient drug costs ----
medcostoth_tbl <- stateval_tbl(medcost_data$other_outpatient_drugs,
                              dist = "gamma")
# head(medcostoth_tbl)

### Hospice costs (after progrssion from line 5) ----
medcosthospice_tbl <- stateval_tbl(medcost_data$hospice,
                              dist = "gamma")
# head(medcosthospice_tbl)

#### Create heSim State Objects ----
#drugcostmod <- create_StateVals(drugcost_tbl, n = n_samples, hesim_data = hesim_dat)
medcostipmod <- create_StateVals(medcostip_tbl, n = n_samples, hesim_data = hesim_dat)
medcostopmod <- create_StateVals(medcostop_tbl, n = n_samples, hesim_data = hesim_dat)
medcostothmod <- create_StateVals(medcostoth_tbl, n = n_samples, hesim_data = hesim_dat)
medcosthospicemod <- create_StateVals(medcosthospice_tbl, n = n_samples, hesim_data = hesim_dat)

## Create Objects for heSim Simulation ----
  ## NOTE: Could move other objects under Create heSim State Objects here for clarity
costmods <- list(Inpatient = medcostipmod,
                 Outpatient = medcostopmod,
                 OtherOutRx = medcostothmod,
                 HOSPICE = medcosthospicemod)


econmod <- IndivCtstm$new(trans_model = transmod,
                          utility_model = utilitymod,
                          cost_models = costmods
)

# Simulation ----
econmod$sim_disease()
econmod$disprog_[, `:=`(transition = paste0(from,'_',to), time_spent= time_stop - time_start), ]
head(econmod$disprog_)  ### individual level results from simulation

## Some useful summary tables ----
disprog <- copy(econmod$disprog_)
LYS <- disprog[,.(time_spent=mean(time_spent)), by =.(sample, strategy_id, from)]
names(LYS)[3] <- "state_id"

disprog <- disprog[,.( progressed = .N,
  mean_time_spent = mean(time_spent), 
  median_time_spent = median(time_spent)
  ),
  by = .(sample, strategy_id, transition, from, to,final)]

### Progression summary ----
summary_dis <- disprog[,.( mean_progressed = mean(progressed),
                           mean_time_spent = mean(mean_time_spent), 
                           med_time_spent = median(mean_time_spent),
                           time_spent05 = quantile(mean_time_spent,.05), 
                           time_spent95 = quantile(mean_time_spent,.95)), 
                       .(strategy_id, from, to, transition)
                       ][order(strategy_id, transition)]

per_line <- summary_dis[, .(mean_entry = sum(mean_progressed)), by = .(strategy_id,from)]
# per_line[,prop_init := mean_entry/n_patients]

summary_dis <- left_join(summary_dis,per_line,by = c("strategy_id","from"))
summary_dis <- summary_dis[,.(strategy_id, from, mean_entry, to, mean_progressed,  
                              mean_time_spent, med_time_spent, time_spent05, time_spent95)]
summary_dis[, `:=`(proportion= ifelse(mean_entry==0,"-", round(mean_progressed/mean_entry,2)),
                   prop_init = ifelse(mean_entry==0,'-', round(mean_entry/n_patients,2 )))
            ]
names(summary_dis) <- c("Strategy", "from", "Entered",
                         "to", "Exited",  "Mean", "Median",
                         "Low 5%", "Top 5%", "Proportion", "% of initial pop")
rm(disprog)
summary_dis <- summary_dis[, `:=`(`From Line` = ifelse(from == 6, "5+", from),
                                  `To Line` = ifelse(to == 7, "D", ifelse(to == 6, "5+",to)),
                                  `Mean Duration` = round(Mean,2),
                                  `Median Duration` = round(Median,2),
                                  `CI(5%)` = paste0("(",round(`Low 5%`,2)," ,",
                                  round(`Top 5%`,2), ")")
                                  )]
summary_dis[, Transition := paste(`From Line`,'to', `To Line`)]
summary_dis <- summary_dis[,.(Strategy, `% of initial pop`  ,`From Line`, `To Line`,
                              Transition, Proportion,
                              `Mean Duration`, `Median Duration`,
                              `CI(5%)`)]

### Progression and survival per line of treatment summary ----
#### Time to progression ----
summary_dis[`To Line`!= "D", .(Strategy, `% of initial pop`,
                               Transition, `Proportion`,
                               `Mean Duration`, `Median Duration`, 
                               `CI(5%)`)] 

#### Death pre-progression ----
summary_dis[`To Line` == "D", .(Strategy, `% of initial pop`,
                               Transition, `Proportion`,
                               `Mean Duration`, `Median Duration`, 
                               `CI(5%)`)] 

rm(summary_dis)

####  Progression free survival ----

disprog <- copy(econmod$disprog_)
disprog <- disprog[,.( mean_time_spent = mean(time_spent), 
                       median_time_spent = median(time_spent)
),
by = .(sample, strategy_id, from)]

summary_dis <- disprog[,.(mean_time_spent = mean(mean_time_spent), 
                          med_time_spent = median(mean_time_spent),
                          time_spent05 = quantile(mean_time_spent,.05), 
                          time_spent95 = quantile(mean_time_spent,.95)), 
                       .(strategy_id, from)
][order(strategy_id)]
 per_line[,prop_init := mean_entry/n_patients]
summary_dis <- left_join(summary_dis, per_line[,.(strategy_id, from, prop_init)],
                         by = c('strategy_id', 'from'))
names(summary_dis) <- c("Strategy", "LOT",
                        "Mean", "Median",
                        "Low 5%", "Top 5%", "prop_init")
summary_dis[, `:=`(`LOT` = ifelse(LOT == 6, "5+", LOT),
                   `% of initial pop` = ifelse(prop_init==0,'-', round(prop_init,2)),
                   `Mean Duration` = round(Mean,2),
                   `Median Duration` = round(Median,2),
                   `CI(5%)` = paste0("(",round(`Low 5%`,2)," ,",
                                     round(`Top 5%`,2), ")") )]
summary_dis[,.(Strategy, `% of initial pop`, LOT,
               `Mean Duration`, `Median Duration`, `CI(5%)`)]

## Some useful summary graphs ----
### Time to progression at  second line ----
ttpl2.data <- econmod$disprog_[from ==2 & to == 3,] 
ttp_data <- lapply(seq(0,30, by =.5), 
                    function(x){data.table(t=x,ttpl2.data[,.(prob = mean(time_spent>x)),
                                                          by =.(sample,strategy_id)])})
ttp_data <- do.call(rbind, ttp_data)
ttp_data <- ttp_data[, .(prob_mean = mean(prob),
                         ci_l = quantile(prob, .025),
                         ci_u = quantile(prob, .975)),
                     by = .(t,strategy_id)]
setorder(ttp_data, strategy_id, t)
set_labels(ttp_data, labels = labs, new_names = c("strategy_name"))
ggplot(ttp_data, aes(x = t, y = prob_mean, col = strategy_name)) +
  geom_line() +  ylim(0,1) +
  geom_ribbon(aes(ymin = ci_l, ymax = ci_u, fill = strategy_name), alpha = .2, linetype=0) +
  xlab("Years") + ylab("Probability") +
  ggtitle("Time to progression at 2nd Line") +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))
# ggsave("output/ttp2ndline.pdf", width = 20, height = 15, units = "cm")
rm(ttpl2.data, ttp_data) # remove objects no longer in use

### Progression free survival at second line ----
pfsl2.data <- econmod$disprog_[from ==2,] 
pfs_data <- lapply(seq(0,30, by =.5), 
                   function(x){data.table(t=x,pfsl2.data[,.(prob = mean(time_spent>x)),
                                                         by =.(sample,strategy_id)])})
pfs_data <- do.call(rbind, pfs_data)
pfs_data <- pfs_data[, .(prob_mean = mean(prob),
                         ci_l = quantile(prob, .025),
                         ci_u = quantile(prob, .975)),
                     by = .(t,strategy_id)]
setorder(pfs_data, strategy_id, t)
set_labels(pfs_data, labels = labs, new_names = c("strategy_name"))
ggplot(pfs_data, aes(x = t, y = prob_mean, col = strategy_name)) +
  geom_line() +  ylim(0,1) +
  geom_ribbon(aes(ymin = ci_l, ymax = ci_u, fill = strategy_name), alpha = .2, linetype=0) +
  xlab("Years") + ylab("Probability") +
  ggtitle("Progression Free Survival at 2nd Line") +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))

ggsave("output/pfs2ndline.pdf", width = 20, height = 15, units = "cm")
rm(pfsl2.data, pfs_data) # remove objects no longer in use

### Overall Survival for those who achieve second line ----

osl2.data <- econmod$disprog_[from!=1, .(OS = sum(time_spent)), 
                              by = .(sample, strategy_id, patient_id)]
prob_data <- lapply(seq(0,30, by =.5), 
                    function(x){data.table(t=x,osl2.data[,.(prob = mean(OS>x)),
                                                         by =.(sample,strategy_id)])})
prob_data <- do.call(rbind, prob_data)
prob_data <- prob_data[, .(prob_mean = mean(prob),
                           ci_l = quantile(prob, .025),
                           ci_u = quantile(prob, .975)),
                       by = .(t,strategy_id)]
setorder(prob_data, strategy_id, t)
set_labels(prob_data, labels = labs, new_names = c("strategy_name"))
ggplot(prob_data, aes(x = t, y = prob_mean, col = strategy_name)) +
  geom_line() +  ylim(0,1) +
  geom_ribbon(aes(ymin = ci_l, ymax = ci_u, fill = strategy_name), alpha = .2, linetype=0) +
  xlab("Years") + ylab("Probability") +
  ggtitle("Overall Survival from 2nd Line") +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))

# ggsave("output/OS2ndline.pdf", width = 20, height = 15, units = "cm")
rm(osl2.data, prob_data) # remove objects that are no longer used

 ## Plot of state probabilities ----
econmod$sim_stateprobs(t = seq(0, 20, .5)) 
summarize_stprobs <- function(stateprobs){
  x <- stateprobs[, .(prob_mean = mean(prob), 
                      low_prob_2.5 = quantile(prob, .025),
                      up_prob_2.5 = quantile(prob, .975)),
                  by = c("strategy_id", "state_id", "t")]
  set_labels(x, labels = labs, new_names = c("strategy_name", "state_name"))
}

stprobs <- summarize_stprobs(econmod$stateprobs_)
ggplot(stprobs, aes(x = t, y = prob_mean, col = strategy_name)) +
  geom_line() +  
  geom_ribbon(aes(ymin = low_prob_2.5, ymax = up_prob_2.5, fill = strategy_name), alpha = 0.2, linetype=0) +
  facet_wrap(~state_name,ncol = 2) + 
  xlab("Years") + ylab("Probability in health state") +
  # scale_color_discrete(name = "Strategy") +
  theme(legend.position = "bottom", 
        legend.title = element_blank())

# ggsave("output/stprobs.pdf", width = 15, height = 25, units = "cm")
rm(stprobs)

# Economics analysis ----
v <- 1/(1+disc) # discount factor

# QALYs model ----
econmod$sim_qalys(dr = disc)
econmod$qalys_

## Summary of QALYs avarages per line of treatent ----
qaly_chart <- chart.LOT(econmod$qalys_,'QALYs')
qaly_chart$plot + ggtitle("Average QALYs") +
  theme(plot.title = element_text(hjust = 0.5) ) 
ggsave("output/qalys.pdf", width = 20, height = 15, units = "cm")

# Medical Costs model ----
econmod$sim_costs(dr = disc)
head(econmod$costs_)

## Treatment cost estimation ----
disprog <- copy(econmod$disprog_)

ptm <- proc.time()
disprog[,TX.cost:=cost_eval(strategy_id, from, time_start, time_stop, disc_fact = v)]
proc.time() - ptm

txcosts <- disprog[,.(costs =mean(TX.cost)), 
                   by = .(sample,  grp_id, strategy_id, from)]
txcosts <- txcosts[, `:=`(category= "Drug", dr = disc)][, .(sample, strategy_id, 
                                                            grp_id, from, dr, category, costs),]
names(txcosts)[names(txcosts) == "from"] <- "state_id"
econmod$costs_ <- rbind(econmod$costs_, txcosts)
rm(txcosts)

## Chart for drug costs ----
tx_chart <- chart.LOT(df = econmod$costs_, categ = 'Drug')
tx_chart$plot + ggtitle("Average Treatment Costs") +
  theme(plot.title = element_text(hjust = 0.5) ) +
  scale_y_continuous(labels = scales::dollar_format(scale = 1e-6, suffix = "M"))
# ggsave("output/txcosts.pdf", width = 20, height = 15, units = "cm")

## Chart for inpatient costs ----
ip_chart <- chart.LOT(df = econmod$costs_, categ = 'Inpatient')
ip_chart$plot + ggtitle("Average Inpatient Costs") +
  theme(plot.title = element_text(hjust = 0.5) ) +
  scale_y_continuous(labels = scales::dollar_format(scale = 1e-3, suffix = "K"))
# ggsave("output/ipcosts.pdf", width = 20, height = 15, units = "cm")

## Chart for outpatient costs ----
op_chart <- chart.LOT(df = econmod$costs_, categ = 'Outpatient')
op_chart$plot + ggtitle("Average Outpatient Costs") +
  theme(plot.title = element_text(hjust = 0.5) ) + 
  scale_y_continuous(labels = scales::dollar_format(scale = 1e-3, suffix = "K"))
# ggsave("output/opcosts.pdf", width = 20, height = 15, units = "cm")

## Chart for total medical costs  ----
medcost_chart <- chart.LOT(econmod$costs_,'Total Medical Expenses')
medcost_chart$plot  + ggtitle("Average Medical Expenses") +
  theme(plot.title = element_text(hjust = 0.5) ) +
  scale_y_continuous(labels = scales::dollar_format(scale = 1e-3, suffix = "K"))
# ggsave("output/medcosts.pdf", width = 20, height = 15, units = "cm")

chart_data <- list(
  QALYS = qaly_chart$chart_data,
  TX = tx_chart$chart_data,
  IP = ip_chart$chart_data,
  OP = op_chart$chart_data,
  MedCosts = medcost_chart$chart_data
)
# writexl::write_xlsx(chart_data, 'output/cost_chart_data.xlsx')
# WriteXLS::WriteXLS(chart_data, ExcelFileName = 'output/cost_chart_data.xlsx')

rm(qaly_chart, tx_chart, ip_chart, op_chart, medcost_chart, chart_data)

## Summary results ----
ce_sim <- econmod$summarize()
format(summary(ce_sim, labels = labs))

# writexl::write_xlsx(summary(ce_sim, labels = labs), path = 'output/summary_costs.xlsx')
# WriteXLS::WriteXLS(summary(ce_sim, labels = labs), ExcelFileName = 'output/summary_costs.xlsx')

cea_pw_out <- cea_pw(ce_sim, comparator = 3, dr_qalys = disc, dr_costs = disc) 
