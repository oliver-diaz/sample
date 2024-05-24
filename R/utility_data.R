# Create QALY input tables ----

utility_data <- data.table(readxl::read_excel("data_raw/Utility_Scores.xlsx",
                                            sheet = 'Summary_v2'))
utility_data <- utility_data[state_id !=7,.(strategy_id, state_id, mean, se)]

save(utility_data, file = "data/utility_data.rda", compress = "bzip2")
