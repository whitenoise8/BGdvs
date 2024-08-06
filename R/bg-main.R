#' Inference for Bernoulli-Gaussian time varying parameter model with dynamic variable selection (! number of covariates should be p > 2). Several settings are available and both Variational Bayes and MCMC can be used.
#'
#' @param y vector of observations (response variable)
#' @param X design matrix (\code{n x p})
#' @param hyper list of hyperparameters (see paper for details)
#' @param meanField type of mean-field factorization (see paper for details). Only "full" is available for now. "group", "joint" will be available soon. If group, group indicators should be given in the \code{hyper} list. If \code{meanField = 'null'}, MCMC is run (slow) and number of posterior samples \code{ndraws} and burn-in \code{nburn} should be provided in \code{hyper}.
#' @param options list containing two fields. \code{sv=0} (no stochastic volatility) \code{sv=1} (stochastic volatility, default) and \code{smooth=0} (no smoothing inclusion probabilities, default) \code{sv=1} (smoothing inclusion probabilities) (see the paper for details). Currently not all the options are supported by all the methods.
#' @param initialize default "null", that is \beta=0. If "OLS", OLS estimates are used to initialize \beta
#' @param Tol_Par tolerance in the relative variation of parameters to detect convergence in Variational Bayes
#' @param Trace if 1, print progress of the algorithm
#' 
#'
#' @return  Variational parameters of the approximated posterior distributions according to the selected method (see the paper for details, the names of the variables match). If MCMC is used, posterior samples after burnin are returned.
#' @export
bg_dvs = function(y,X,hyper,
                  meanField='full',
                  options=list(sv=1,smooth=0),
                  initialize='null',
                  Tol_Par = 0.01,Trace = 0) {
  
  if (meanField == 'full') {
    mod = BGTVP_independent(y,X,hyper,
                            options = options,
                            initialize = initialize,
                            Tol_Par = Tol_Par,Trace = Trace)
  }
  
  # if (meanField == 'group') {
  #   group_index = hyper$group_index
  #   mod = BGTVP_group_dependent(y,X,group_index,hyper,
  #                           options = options,
  #                           initialize = initialize,
  #                           Tol_Par = Tol_Par,Trace = Trace)
  # }
  # 
  # if (meanField == 'joint') {
  #   mod = BGTVP_dependent(y,X,hyper,
  #                           options = options,
  #                           initialize = initialize,
  #                           Tol_Par = Tol_Par,Trace = Trace)
  # }
  
  if (meanField == 'none') {
    ndraws = hyper$ndraws
    nburn = hyper$nburn
    
    mod = NULL
    
    if ((options$sv==1)&(options$smooth==0)) mod = BGTVPSV_MCMC(y,X,hyper,ndraws,nburn,initialize=initialize,Trace = Trace)
    
    if ((options$sv==0)&(options$smooth==0)) mod = BGTVP_MCMC(y,X,hyper,ndraws,nburn,initialize=initialize,Trace = Trace)
    
    if (is.null(mod)) print("Options setting not available!")
  }
  
  mod
}