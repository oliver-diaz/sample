# Change Log 
# Date       Programmer     Description
#  20230828     OD            Created Code
#  20230829     tjf           Documentation and Review

####

# Treatment cost evaluation ----
## func.tx.cost: Calculates the discounted costs of a treatment  ----
### Input parameters:
### T_start: Within model time of entry into the state
### T_stop: Within model time of progression/transition out of state
### tau_0:  Time of drug initiation after entering state
### tau_1:  Time of drug cessation after entering state
### Cost:   annualized cost of the treatment
### v:      discount factor where discount rate = (1-v)/v

func.tx.cost <- function(T_start,T_stop ,tau_0, tau_1,Cost, v){
  if(v!=1){
    return(ifelse(Cost==0,0,
                  ifelse(tau_0 != tau_1, Cost * (v^(pmin(T_start+ tau_0 ,T_stop)) - v^(pmin(T_start+ tau_1 ,T_stop)))/(-log(v)),
                         (T_start+ tau_0 < T_stop)* Cost * v^(T_start+tau_0))))
  } else{
    return(ifelse(Cost==0,0,
                  ifelse(tau_0 != tau_1, 
                         Cost * (pmin(T_start+ tau_1 ,T_stop) - pmin(T_start+ tau_0 ,T_stop)),
                         Cost*(T_start+ tau_0 < T_stop))))
  }
}

#### A wrapper function to vectorize and stage entry and exit times from states 
func.tx.cost <- Vectorize(func.tx.cost, vectorize.args = c("T_start", "T_stop"))

## cost_eval: calculate PV (discounted) drug costs using per-protocol dosing schedules and modeled transitions ----
### strategy: treatment strategy id specified in model
### lineTX: line of therapy ("from") indicated in model
### T_start: Time of entry into current state simulated in model
### T_stop:  Time of exit out of current state simulated by model
### disc_fact: discount factor where discount rate = (1-v)/v
### It used treatment cost tables that contain the drugs used in any specific strategy and line 
### of treatment. These tables are contained in the list cost_tx_seq, which is created or loaded at
### the begining of the model_structure code.

cost_eval <- function( strategy, lineTX,T_start, T_stop, 
                       disc_fact = 1/(1.03), tbl = cost_tx_seq){
  sum(tbl[[strategy]][[lineTX]][,
                                TX_cost:= func.tx.cost(T_start,
                                                       T_stop, 
                                                       t_init, 
                                                       t_stop, 
                                                       Cost, 
                                                       disc_fact)]$TX_cost)
}

#### A wrapper function to vectorize arguments used by cost_eval function 
cost_eval <- Vectorize(cost_eval, 
                       vectorize.args = c("strategy", "lineTX", "T_start", "T_stop"))
