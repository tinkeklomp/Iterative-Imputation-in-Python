# functions to apply missing data methods

# complete case analysis
apply_CCA <- function(amp) {
  # list-wise deletion
  est <- na.omit(amp$amp) %>% 
    # fit regression 
    lm(Y ~ X1 + X2 + X3 + X4, .) %>% 
    # clean results
    broom::tidy(conf.int = TRUE) %>% 
    # choose estimates
    select(term, estimate, conf.low, conf.high) %>% 
    # add method name and missingness
    cbind(method = "CCA", mech = amp$mech, prop = amp$prop, .) 
  # output
  return(est)
}

# MICE imputation
apply_MICE <- function(amp) {
  # imputation with MICE
  imp <- mice::mice(amp$amp, method = "norm", printFlag = FALSE)
  # fit regression on each imputation
  est <- with(imp, lm(Y ~ X1 + X2 + X3 + X4)) %>% 
    # pool results
    mice::pool() %>% 
    # clean results
    broom::tidy(conf.int = TRUE) %>% 
    # select estimates
    select(term, estimate, conf.low, conf.high) %>% 
    # add method name
    cbind(method = "MICE", mech = amp$mech, prop = amp$prop, .)
  # output
  return(est)
}

# Python imputation
# impute data with python
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

# combine into one function
apply_methods <- function(amps, betas) {
  # apply CCA to each incomplete dataset
  CCA <- purrr::map_dfr(amps, ~{apply_CCA(.)})
  # impute with MICE and estimate effects
  MICE <-  purrr::map_dfr(amps, ~{apply_MICE(.)})
  # impute with Python and estimate effects
  IT_IMP <- purrr::map_dfr(amps, ~{apply_PYTHON(.)})
  # combine estimates 
  ests <- rbind(CCA, MICE, IT_IMP) %>% 
    cbind(truth = c(0, betas))
  # output
  return(ests)
}
