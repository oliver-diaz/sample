#' ############################################
#' #      Project: JAN56830                   #
#' #      Hazard ratio tables (default)       #
#' #                                          #
#' ############################################


# Create default hazard ratio (HR) ranges for each tansition in the model ----
# Tables with HR ranges for each transition (t1,...,t10)
## t1: "P1-> R1_P2", 
## t2: "P1 -> Death",
## t3: "R1_P2 -> R2_P3",
## t4: "R1_P2 -> Death", 
## t5: "R2_P3 -> R3_P4", 
## t6: "R2_P3 -> Death",
## t7: "R3_P4 -> R4_P5",
## t8: "R3_P4 -> Death", 
## t9: "R4_P5 -> R5", 
## t10: "R4_P5 -> Death", 
## t11: "R5 -> Death"

# Sequence strategies 
## st1 = Sequence 1
## st2 = Sequence 2
## st3 = Sequence 3
## notst3 = Either sequence 1 or 2
## notst2 = Either sequence 1 or 3

kol.HR <- list(
  t1.hr = data.frame(low = c(notst3 = .21, st3 = 1), up= c(notst3 = .95, st3 = 1)), 
  t3.hr = data.frame(low = c(st1 = 1, st2= .66, st3= 1), up=c(st1= 1, st2=1.09, st3= 1)),
  t5.hr = data.frame(low = c(notst2 =.96, st2= .4), up=c(notst2= 1, st2=.79)),
  t7.hr = data.frame(low = c(st1 =1, st2= .05, st3= .16), up=c(st1=1, st2=.11, st3= 1)),
  t9.hr = data.frame(low = c(st1 =.35, st2= 1.51, st3= 1.04), up=c(st1=.53, st2=2.21, st3= 1.42))
)

save(kol.HR, file = "data/kol.HR.rda", compress = "bzip2")

