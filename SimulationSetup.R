# simulation script for ADS thesis
# original script by Hanne Oberman

########################
### SETUP SIMULATION ###
########################

# packages
library(dplyr)
library(mvtnorm)
library(mice)
library(miceadds)
library(reticulate)

# functions
miceadds::source.all("./functions")

# randomness
set.seed(1)

# parameters
n_sim <- 1000
n_obs <- 200
betas <- c(-0.5, -0.1, 0.1, 0.5)
mis_mech = c("MCAR", "MAR")
mis_prop = c(0.1, 0.25, 0.5)


######################
### IMPUTE DATASET ###
######################

# impute data with IterativeImputer
complete_py <- function(incomplete, m) {
  #impute amputed dataset with IterativeImputer functions fit() and transform()
  it_imp_fit <- it_imp$fit_transform(incomplete)
  #convert to dataframe
  py_it_imp <- as.data.frame(it_imp_fit)
  #reassign column names
  colnames(py_it_imp) <- c("Y", "X1", "X2", "X3", "X4")
  py_it_imp <- cbind(.imp = m, py_it_imp)
  return(py_it_imp)
}

apply_PYTHON <- function(amp){
  #convert cells with missing values from na to NaN
  amp_nan <- amp$amp
  amp_nan[is.na(amp_nan)]<- NaN
  
  # fit regression on each imputation
  long <- purrr::map_dfr(1:5, ~complete_py(amp_nan, m = .x)) %>% 
    rbind(cbind(.imp = 0, amp$amp), .)
  imp <- as.mids(long)
  est <- with(imp, lm(Y ~ X1 + X2 + X3 + X4)) %>%
    # pool results
    mice::pool() %>%
    # clean results
    broom::tidy(conf.int = TRUE) %>%
    # select estimates
    select(term, estimate, conf.low, conf.high) %>%
    # add method name
    cbind(method = "IT_IMP", mech = amp$mech, prop = amp$prop, .)
  # output
  return(est)
}


##############################################
### COMBINE ALL FUNCTIONS INTO ONE FUCTION ###
##############################################

simulate_once <- function(n_obs, betas, mis_mech, mis_prop) {
  # generate incomplete data
  amps <- create_data(
    sample_size = n_obs,
    effects = betas,
    mechanisms = mis_mech,
    proportions = mis_prop
  )
  # estimate regression coefficients
  ests <- apply_methods(amps)
  # output
  return(ests)
}


######################
### RUN SIMULATION ###
######################

# repeat the simulation function n_sim times
results_raw <- replicate(
  n_sim, 
  simulate_once(n_obs, betas, mis_mech, mis_prop),
  simplify = FALSE
)
# save raw results
saveRDS(results_raw, "./Results/raw.RDS")


########################
### EVALUATE RESULTS ###
########################

# function(s) to evaluate the estimates
evaluate_est <- function(results) {
  performance <- purrr::map_dfr(results, ~{
    mutate(.x,
           bias = truth - estimate,
           cov = conf.low <= truth & conf.high >= truth,
           ciw = conf.high - conf.low,
           .keep = "unused"
    )
  })
  return(performance)
}

# calculate bias, coverage rate and CI width
performance <- evaluate_est(results_raw)

# simulation results across all conditions
res_all_cond <- performance %>% 
  group_by(method) %>% 
  summarise(across(c(bias, cov, ciw), mean))
res_all_cond

# simulation results split by condition
res_split_cond <- performance %>% 
  group_by(method, mech, prop) %>% 
  summarise(across(c(bias, cov, ciw), mean))
res_split_cond

# simulation results split by condition and regression coefficient
res_cond_reg <- performance %>% 
  group_by(method, mech, prop, term) %>% 
  summarise(across(c(bias, cov, ciw), mean))
res_cond_reg