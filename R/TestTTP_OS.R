LOT <- 2

### TTP
ttpl2.data <- econmod$disprog_[from ==LOT & to!="7",] 
ttp_data <- lapply(seq(0,30, by =.2), 
                   function(x){data.table(t=x,ttpl2.data[,.(prob = mean(time_spent > x)),
                                                         by =.(sample,strategy_id)])})
ttp_data <- do.call(rbind, ttp_data)
ttp_data <- ttp_data[, .(prob_mean = mean(prob),
                         ci_l = quantile(prob, .025),
                         ci_u = quantile(prob, .975)),
                     by = .(t,strategy_id)]
setorder(ttp_data, strategy_id, t)
set_labels(ttp_data, labels = labs, new_names = c("strategy_name"))

### PFS
pfsl2.data <- econmod$disprog_[from ==LOT,] 
pfs_data <- lapply(seq(0,30, by =.2), 
                   function(x){data.table(t=x,pfsl2.data[,.(prob = mean(time_spent > x)),
                                                         by =.(sample,strategy_id)])})
pfs_data <- do.call(rbind, pfs_data)
pfs_data <- pfs_data[, .(prob_mean = mean(prob),
                         ci_l = quantile(prob, .025),
                         ci_u = quantile(prob, .975)),
                     by = .(t,strategy_id)]
setorder(pfs_data, strategy_id, t)
set_labels(pfs_data, labels = labs, new_names = c("strategy_name"))

### OS
osl2.data <- econmod$disprog_[from >=  LOT, .(OS = sum(time_spent)), 
                              by = .(sample, strategy_id, patient_id)]
prob_data <- lapply(seq(0,30, by =.2), 
                    function(x){data.table(t=x,osl2.data[,.(prob = mean(OS>x)),
                                                         by =.(sample,strategy_id)])})
prob_data <- do.call(rbind, prob_data)
prob_data <- prob_data[, .(prob_mean = mean(prob),
                           ci_l = quantile(prob, .025),
                           ci_u = quantile(prob, .975)),
                       by = .(t,strategy_id)]
setorder(prob_data, strategy_id, t)
set_labels(prob_data, labels = labs, new_names = c("strategy_name"))

#### Test PFS - OS
time_ <- unique(pfs_data$t)
strategy <- 1

plot(time_, pfs_data[strategy_id == strategy,  prob_mean],
     type ="l", col='red',
     xlab = "Time (years)", 
     ylab = "Probability",
     main = paste0('Survival Curves starting at LOT ', LOT, ': Sequence ', strategy ))
lines(time_, ttp_data[strategy_id == strategy,  prob_mean],
     type ="l", col='magenta')
lines(time_, prob_data[strategy_id == strategy, prob_mean],
      type ="l", col='darkred')
legend(25, 1, legend = c("PFS", "TTP", "OS"), fill =c("red","magenta","darkred"))

strategy <- 2
plot(time_, pfs_data[strategy_id == strategy, prob_mean],
     type ="l", col='green',
     xlab = "Time (years)", 
     ylab = "Probability",
     main = paste0('Survival Curves starting at LOT ', LOT, ': Sequence ', strategy ))
lines(time_, ttp_data[strategy_id ==strategy, prob_mean],
     type ="l", col='olivedrab')
lines(time_, prob_data[strategy_id == strategy, prob_mean],
     type ="l", col='darkgreen')
legend(25, 1, legend = c("PFS", "TTP", "OS"), fill =c("green","olivedrab","darkgreen"))

strategy <- 3
plot(time_, pfs_data[strategy_id == strategy, prob_mean],
     type ="l", col='skyblue',
     xlab = "Time (years)", 
     ylab = "Probability",
     main = paste0('Survival Curves starting at LOT ', LOT, ': Sequence ', strategy ))
lines(time_, ttp_data[strategy_id == strategy, prob_mean],
      type ="l", col='blue')
lines(time_, prob_data[strategy_id == strategy, prob_mean],
      type ="l", col='darkblue')
legend(25, 1, legend = c("PFS", "TTP", "OS"), fill =c("skyblue","blue","darkblue"))
