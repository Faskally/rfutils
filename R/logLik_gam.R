logLik_gam <- function(model) {
  
  stopifnot(
    class(model)[1] %in% c("gam", "bam"),
    model$method == "ML"
  )
  
  log_likelihood <- - model$gcv.ubre
  attributes(log_likelihood) <- NULL
  
  
  # fixed df
  # rank gives total number of 'covariates'
  # df_smooth gives those columns that are not fixed effects
  
  if (length(model$smooth) == 0) {
    df_smooth <- 0
  } else {
    df_smooth <- sapply(
      model$smooth,
      function(x) if (x$fixed) 0 else sum(x$rank)
    )
  }
  
  df_fixed <- model$rank - sum(df_smooth)
  
  
  # random df
  # smoothing paramters and scale estimate
  
  df_random <- length(model$sp) + model$scale.estimated
  
  df <- df_fixed + df_random
  
  attributes(log_likelihood) <- list(df = df)
  
  class(log_likelihood) <- "logLik"
  
  log_likelihood
}