# function(s) to evaluate the estimates
### [YOUR FUNCTION HERE] ###

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
