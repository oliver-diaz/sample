#' ################################################
#' #      Project: JAN56830                       #
#' #      Treatment scheduling and cost           #
#' #                                              #
#' ################################################

# Create Cost tables ----

## Description:  The following code generates 2-level list of data input tables.
## These data tables reflect the 'per-protocol' dosage and costs of the drugs.
## The dimensions of this list are:
##  1) Sequences: 1-3 per the current Janssen Model
##  2) Line of Therapy (LOT): 1-7 consistent with the  states of the Janssen Model (See below and t-mat from model)# States
### 1. P1:    preprogress on 1L therapy
### 2. R1_P2: relapse on 1L therapy-preprogression, on 2L therapy 
### 3. R2_P3: relapse on 2L therapy-preprogression, on 3L therapy
### 4. R3_P4: relapse on 3L therapy-preprogression, on 4L therapy
### 5. R4_P5: relapse on 4L therapy-preprogression, on 5L therapy
### 6. R5: relapse on 5L therapy
### 7. Death

##  The data table for each LOT has 4 required elements:
##  Drugs:  The name or highlight of the drug
##  NOTE: For combination therapies, there will be one line per drug in the combination
##  Cost:   The annual cost of the drug.  This is the TOTAL COST over a 12 month period 
##.         for receiving the drug per-protocol
##  t_init: This is the time of drug initiation relative to the start of this treatment.
##  t_stop: The time of drug termination relative to the start of treatment
##

# Create baseline table of annualized costs and drug schedule ----
baselineTX_data <- data.table::data.table(readxl::read_excel("data_raw/Treatment_costs.xlsx",
                                                             sheet = 'R_TreatmentCosts',
                                                             range = 'A1:L113'))
names(baselineTX_data)[c(3,12)] <- c("lineTX", "Cost") # Cost = Annualized cost per administration
save(baselineTX_data, file = "data/baselineTX_data.rda", compress = "bzip2")

# Baseline WAC ----
base_wac <- data.table::data.table(readxl::read_excel("data_raw/Treatment_costs.xlsx",
                                                      sheet = 'R_TreatmentCosts',
                                                      range = 'P1:T15'))
save(base_wac, file = "data/base_wac.rda", compress = "bzip2")

# Separate cost and schedules by sequence strategy and LOT ----
cost_tx_seq <- cost.tx.data.fun(baselineTX_data)
save(cost_tx_seq, file = "data/cost_tx_seq.rda", compress = "bzip2")

#' ##################### END #########################

#' cost_tx_seq[[1]] <- list(
#'   line1 = data.table(Drugs= c("D",  "V", "R", "d", "ASCT", "DR"),
#'                      Cost= runif(6,200,1000),
#'                      t_init = c(0,0,3/12,3/12, 5/12,5/12),
#'                      t_stop = c(2/12, 2/12, 5/12, 5/12, 9999,9999)
# 
# 
#' # set.seed(2189) # fix seed to make results reproducible.
#' ## Sequence 1 ----
#' ### Seq1: (DVRd + ASCT + DR )  > DKd > EloPd  > Sd > Cilta-cel ----
#' 
#' cost_tx_seq[[1]] <- list(
#'   line1 = data.table(Drugs= c("D",  "V", "R", "d", "ASCT", "DR"),
#'                      Cost= runif(6,200,1000),
#'                      t_init = c(0,0,3/12,3/12, 5/12,5/12),
#'                      t_stop = c(2/12, 2/12, 5/12, 5/12, 9999,9999)
#'   ),
#'   line2  = data.table(Drugs= c("D", "Kd"),
#'                       Cost= runif(2,200,1000),
#'                       t_init = c(0,0),
#'                       t_stop = c(9999,9999)
#'   ),
#'   line3 = data.table(Drugs= c("Elo", "Pd"),
#'                      Cost= runif(2,200,500),
#'                      t_init = c(0,1/12),
#'                      t_stop = c(9999,9999)
#'   ),
#'   line4 = data.table(Drugs= c("Sd"),
#'                      Cost= runif(1,500,1000),
#'                      t_init = 0,
#'                      t_stop = 9999
#'   ),
#'   line5 = data.table(Drugs= c("Cilta", "cel"),
#'                      Cost= c(2000, 5000),
#'                      t_init = c(1, 0),
#'                      t_stop = c(9999,0)
#'   ),
#'   line6 = data.table(Drugs= "paliative",
#'                      Cost= 10,
#'                      t_init = 0,
#'                      t_stop = 9999
#'   ),
#'   line7 = data.table(Drugs= "none",
#'                      Cost= 0,
#'                      t_init = 0,
#'                      t_stop = 9999
#'   )
#' )
#' 
#' ## Sequence 2 ----
#' ### Seq2: (DVrd + ASCT + DR )  > Cilta-cel  > TalDara > Tec  > SKd ----
#' 
#' cost_tx_seq[[2]] <- list(
#'   line1 = cost_tx_seq[[1]][[1]],
#'   line2 = cost_tx_seq[[1]][[5]],
#'   line3 = data.table(Drugs= c("Tal", "Dara"),
#'                      Cost= runif(2,200,500),
#'                      t_init = c(0,1/12),
#'                      t_stop = c(9999,9999)
#'   ),
#'   line4 = data.table(Drugs= c("Tec"),
#'                      Cost= c(2000),
#'                      t_init = 0,
#'                      t_stop = 9999
#'   ),
#'   line5 = data.table(Drugs= c("S", "Kd"),
#'                      Cost= runif(2,500,700),
#'                      t_init = c(0,0),
#'                      t_stop = c(9999,9999)
#'   ),
#'   line6 = cost_tx_seq[[1]][[6]],
#'   line7 = cost_tx_seq[[1]][[7]]
#' )
#' 
#' ## Sequence 3 ----
#' ### Seq3: (DVRd + DRd) > Tec > EloPd  > Tal  > Kd ----
#' 
#' cost_tx_seq[[3]] <- list(
#'   line1 = data.table(Drugs = c("D",  "V", "R", "d1", "DR", "d2"),
#'                      Cost = c(834.9673, 444.1380, 448.3985, 895.2274,
#'                               1283.366, 448.3985),
#'                      t_init = c(0, 0, 3/12, 3/12, 5/12, 5/12),
#'                      t_stop = c(2/12, 2/12, 5/12, 5/12, 9999, 9999)
#'   ),
#'   line2 = cost_tx_seq[[2]][[4]],
#'   line3 = cost_tx_seq[[1]][[3]],
#'   line4 = data.table(Drugs= "Tal",
#'                      Cost= 243.1151,
#'                      t_init = 0,
#'                      t_stop = 9999
#'   ),
#'   line5 = data.table(Drugs= c("Kd"),
#'                      Cost= 585.2421,
#'                      t_init = 0,
#'                      t_stop = 9999
#'   ),
#'   line6 = cost_tx_seq[[1]][[6]],
#'   line7 = cost_tx_seq[[1]][[7]]
#' )
#' 
#' save(cost_tx_seq, file = "data/cost_tx_seq.rda", compress = "bzip2")