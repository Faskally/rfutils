#' logLik_gam
#'
#' Log likelihood function for use in model selection of mixed models using gam or bam
#' Here’s a function for extracting the log likelihood and (harder) the degrees of freedom 
#' out of a gam or bam model.  The degrees of freedom are the number of fixed effects plus 
#' the number of variances (not the edf as is usually reported by gam) and assumes the 
#' models are fitted by ML.
#' 
#' The function has been tested with a range of smoothers, but as good check would involve 
#' fitting models in gam and also in gamm4 and check the df are the same. 
#'
#' @param model detail
#'
#' @references
#'
#' @details
#'
#' @export


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