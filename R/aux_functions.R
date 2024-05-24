#' ############################################
#' #      Project: JAN56830                   #
#' #      Auxiliary functions                 #
#' #                                          #
#' ############################################

#  Convert a vector to a data.table with a single column for the intercept ----

vec_to_dt <- function(v, n = NULL){
  if (length(v) == 1) v <- rep(v, n_samples)
  dt <- data.table(v)   
  colnames(dt) <- "cons"
  return(dt)
}

# Loads data inputs  and creates data files when they are missing ----
load_create_fun <- function(){
  
  if(!(file.exists("data/cost_tx_seq.rda")&
       file.exists("data/base_wac.rda")&
       file.exists("data/baselineTX_data.rda"))){
    source("R/tx_cost_data.R") # creates baseline treatment costs data
  }
  if(!file.exists("data/utility_data.rda")){
    source("R/utility_data.R") # creates baseline utilities input data
  }
  if(!file.exists("data/medcost_data.rda")){
    source("R/medcost_data.R") # create baseline medical costs input data
  }
  if(!(file.exists("data/weilph.data.rda") &
       file.exists("data/nonMM.data.rda"))){
    source("R/surv_curves_alt.R") # create reference survival data
  }
  if(!file.exists("data/kol.HR.rda")){
    source("R/HR_data.R") # create KOL HR data
  }
}


# Add empty line at every x rows in a table ----
linesep<-function(x,y=character()){
  if(!length(x))
    return(y)
  linesep(x[-length(x)], c(rep('',x[length(x)]-1),'\\addlinespace',y))  
}


# Chart of averge cost per line of treatment ----
## For QALYs or LYs, use df = econmod$qalys_ 
## For Inpatient, Outpatient, Drug, OtherOutRx, HOSPICE, use df = econmod$qalys_

chart.LOT <- function(df = econmod$costs_, categ = 'Drug' ){
  if(categ == 'QALYs'){
    df_summary <- df[,.(mean = mean(qalys)), 
                     by =.(strategy_id, state_id, dr)]
    y_label <- 'Mean QALYs'
    disc <- min(df_summary$dr)
  } else if(categ=='LYs'){
    df_summary <- df[,.(mean = mean(lys)), 
                     by =.(strategy_id, state_id, dr)]
    y_label <- 'LYs'
    disc <- min(df_summary$dr)
  } else if(categ =='Other Medical Costs'){
    newdf = econmod$costs_[category != 'Drug', .(costs = sum(costs)),
                           by =.(sample,strategy_id,state_id, grp_id,dr)]
    df_summary <- newdf[,.(mean = mean(costs)), 
                        by =.(strategy_id, state_id, dr)]
    y_label <- "Cost (in USD)"
    disc <- min(df_summary$dr)
  } else if(categ == 'Total Costs'){
    df_summary <- df[ ,.(total = sum(costs)), 
                      by = .(sample, strategy_id,grp_id, state_id, dr)]
    df_summary <- df_summary[, .(mean = mean(total)),
                             by =.(strategy_id, state_id, dr)]
    y_label <- "Cost (in USD)"
    disc <- min(df_summary$dr)
  } else {
    df_summary <- df[category == categ,
                     .(mean = mean(costs)), 
                     by =.(strategy_id, state_id, dr)]
    y_label <- "Cost (in USD)"
    disc <- min(df_summary$dr)
  }
  # labs0 <- list("state_id" = c("L1"=1, "L2"=2,"L3"=3,"L4"=4, "L5"=5, "L5+"= 6))
  # set_labels(df_summary, labels = labs0, new_names = "Line_TX")
  set_labels(df_summary, labels = labs, new_names = c("strategy_name", "state_name"))
  chart <- ggplot(df_summary[dr == disc],
                  # aes(x = strategy_name, y = mean, fill = Line_TX)) + 
                  aes(x = strategy_name, y = mean, fill = state_name)) + 
    geom_bar(stat = "identity") +
    scale_fill_discrete(name = "") +
    xlab("Strategy") + ylab(y_label) 
  return(list(chart_data = df_summary, plot = chart))
}


# Plot state probability over time ----
## This function adds a 5% CI band around the mean state probability overtime.
stprobsPlot <- function(state_nm = "P1", alpha = .2){
 p <-  ggplot(stprobs[state_name == state_nm,], aes(x = t, y = prob_mean, col = strategy_name)) +
   geom_line() +  ylim(0,1) +
   geom_ribbon(aes(ymin = low_prob_2.5, ymax = up_prob_2.5, fill = strategy_name), alpha = alpha, linetype=0) +
   xlab("Years") + ylab("Probability in health state") +
   theme(legend.position = "bottom", 
         legend.title = element_blank())
}

