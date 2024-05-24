#' ############################################
#' #      Project: JAN56830                   #
#' #      Medical cost data                   #
#' #                                          #
#' ############################################

# Create baseline medical costs tables ----

## Inputs are read from file "Non_TreatmentCosts.xlsx" in the directory 
## "data_raw", sheet 'CompilationWithYearVariation'

## Hospital costs ----
inpat.df <- data.table(readxl::read_excel(
  "data_raw/Non_TreatmentCosts.xlsx", sheet = 'CompilationWithYearVariation',
  range = "a2:e20"
) )
inpat.df <- inpat.df[,.(strategy_id, state_id, mean, se)]

## Outpatient costs ----
outpat.df <- data.table(readxl::read_excel(
  "data_raw/Non_TreatmentCosts.xlsx", sheet = 'CompilationWithYearVariation',
  range = "a23:e41"
) )
outpat.df <- outpat.df[,.(strategy_id, state_id, mean, se)]

## Other outpatient Rx costs ----
oth_out_rx.df <- data.table(readxl::read_excel(
  "data_raw/Non_TreatmentCosts.xlsx", sheet = 'CompilationWithYearVariation',
  range = "a44:e62"
) )
oth_out_rx.df <- oth_out_rx.df[,.(strategy_id, state_id, mean, se)]

## Hospice costs
hospice.df <- data.table(readxl::read_excel(
  "data_raw/Non_TreatmentCosts.xlsx", sheet = 'CompilationWithYearVariation',
  range = "a65:e83"
) )
hospice.df <- hospice.df[,.(strategy_id, state_id, mean, se)]

medcost_data <- list(
  inpatient = inpat.df,
  outpatient = outpat.df,
  other_outpatient_drugs = oth_out_rx.df,
  hospice = hospice.df
)

save(medcost_data, file = "data/medcost_data.rda", compress = "bzip2")
