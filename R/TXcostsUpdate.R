#' ################################################################
#' #      Project: JAN56830                                       #
#' #      Create treatment cost and schedule data                 #
#' #      Update annualized costs from user provided WAC costs.   #
#' #                                                              #
#' ################################################################


######## Create TX cost and scheduling table

# Calculate annualized TX cost ----
annual.cost <- function(wac,h,method_admin,AvgDose_mg, bsa.weight_flag,
                        weight = 75, bsa = 1.7){
  ifelse(h==0, (method_admin == 'IV')*( (bsa.weight_flag == 1)*weight +  
                                          (bsa.weight_flag == 0)*bsa ) + (method_admin != 'IV'),
         ( (method_admin == 'IV')*( (bsa.weight_flag == 1)*weight +  
                                      (bsa.weight_flag == 0)*bsa ) + (method_admin != 'IV') )/h
  ) * AvgDose_mg* wac
  
}

# Create updated treatment cost and schedule table ----
#' It uses prices set by user
tx.tbl.fun <- function(wac, df.TX){
 TXupdate <- left_join(df.TX, wac, by = c("Drug_Code" = "Abbr"))
 names(TXupdate)[which(names(TXupdate) == 'WAC')] <- "newWAC"
 TXupdate[, Cost := annual.cost(wac = newWAC, h,method_admin,AvgDose_mg, bsa.weight_flag)]
 return(TXupdate[,.(LOT, lineTX, Drug_Code, Drug_Cost_Code,
                    t_init, t_stop, Cost)])
}

# Creates a list of separate treatments costs and schedule by sequence and LOT ----
cost.tx.data.fun <- function(drug.df){
  cost_tx_seq <- vector('list', length = 3) 
  
  ## Sequence 1 ----
  ### Seq1: (DVRd + ASCT + DR )  > DKd > EloPd  > Sd > Cilta-cel ----
  cost_tx_seq[[1]] <- list(
    line1 = drug.df[LOT ==1 & lineTX == "DVRd + ASCT + DR", .(Drug_Cost_Code, Cost, t_init, t_stop)],
    line2 = drug.df[LOT == 2 & lineTX == "DKd", .(Drug_Cost_Code, Cost, t_init, t_stop)],
    line3 = drug.df[LOT == 3 & lineTX == "EloPd", .(Drug_Cost_Code, Cost, t_init, t_stop)],
    line4 = drug.df[LOT == 4 & lineTX == "Sd", .(Drug_Cost_Code, Cost, t_init, t_stop)],
    line5 = drug.df[LOT == 5 & lineTX == "Cilta-cel", .(Drug_Cost_Code, Cost, t_init, t_stop)],
    line6 = data.table(Drug_Cost_Code = "6_Paliative", Cost= 0,
                       t_init = 0, t_stop = 99),
    line7 = data.table(Drug_Cost_Code= "none", Cost= 0,
                       t_init = 0,
                       t_stop = 99)
  )
  ## Sequence 2 ----
  ### Seq2: (DVrd + ASCT + DR )  > Cilta-cel  > TalDara > Tec  > SKd ----
  cost_tx_seq[[2]] <- list(
    line1 = drug.df[LOT ==1 & lineTX == "DVRd + ASCT + DR", .(Drug_Cost_Code, Cost, t_init, t_stop)],
    line2 = drug.df[LOT == 2 & lineTX == "Cilta-cel", .(Drug_Cost_Code, Cost, t_init, t_stop)],
    line3 = drug.df[LOT == 3 & lineTX == "TalD", .(Drug_Cost_Code, Cost, t_init, t_stop)],
    line4 = drug.df[LOT == 4 & lineTX == "Tec", .(Drug_Cost_Code, Cost, t_init, t_stop)],
    line5 = drug.df[LOT == 5 & lineTX == "SKd", .(Drug_Cost_Code, Cost, t_init, t_stop)],
    line6 = data.table(Drug_Cost_Code = "6_Paliative", Cost= 0,
                       t_init = 0, t_stop = 99),
    line7 = data.table(Drug_Cost_Code= "none", Cost= 0,
                       t_init = 0,
                       t_stop = 99)
  )
  ## Sequence 3 ----
  ### Seq3: (DVRd + DRd) > Tec > EloPd  > Tal  > Kd ----
  cost_tx_seq[[3]] <- list(
    line1 = drug.df[LOT ==1 & lineTX == "DVRd + DRd", .(Drug_Cost_Code, Cost, t_init, t_stop)],
    line2 = drug.df[LOT ==2 & lineTX == "Tec", .(Drug_Cost_Code, Cost, t_init, t_stop)],
    line3 = drug.df[LOT ==3 & lineTX == "EloPd", .(Drug_Cost_Code, Cost, t_init, t_stop)],
    line4 = drug.df[LOT ==4 & lineTX == "Tal", .(Drug_Cost_Code, Cost, t_init, t_stop)],
    line5 = drug.df[LOT ==5 & lineTX == "Kd", .(Drug_Cost_Code, Cost, t_init, t_stop)],
    line6 = data.table(Drug_Cost_Code = "6_Paliative", Cost= 0,
                       t_init = 0, t_stop = 99),
    line7 = data.table(Drug_Cost_Code= "none", Cost= 0,
                       t_init = 0,
                       t_stop = 99)
  )
  return(cost_tx_seq)
}

