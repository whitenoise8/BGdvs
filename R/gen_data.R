#' Generate data from the model and the latent coefficients (time-varying parameters) with dynamic sparsity
#'
#' @param n length of the series
#' @param p_in number of coefficients always included in the model
#' @param p_dyn number of coefficients dynamically included in the model
#' @param p_low number of coefficients included in the model only for a short period (low signal)
#' @param p_null number of coefficients always zero
#' @param lam average length of non-zero periods for \code{p_dyn} coefficients (should have same length!)
#' @param delta average length of non-zero period in the low signal case \code{p_low} coefficients (should have same length!)
#' @param thresh threshold for a parameter to be non-zero#' 
#' @param eta2 variance of the random walk process for \beta
#' @param s2 variance of the error term
#' @param nu2 if \code{NULL} (default) homoskedastic errors are considered. If a value is assigned, then stochastic volatility is included and \code{nu2} is the variance of the random walk process 
#' 
#'
#' @return  \code{y}: simulated observations
#' @return  \code{X}: simulated design matrix \code{n x p} (independent standard normally distributed columns) 
#' @return  \code{beta}: the simulated latent coefficients \beta (a matrix \code{p x n})
#' @return  \code{gamma}: the simulated latent indicators \gamma (a matrix \code{p x n})
#' @return  \code{s2}: variance of the errors
#' @export
generateData = function(n,p_in,p_dyn,p_low,p_null,lam,delta,thresh=0.5,eta2=0.1,s2=0.25,nu2=NULL) {
  
  TVPparams = getTVPparams(n,p_in,p_dyn,p_low,p_null,eta2,lam,delta,thresh)
  p = p_in+sum(p_dyn)+p_low+p_null
  
  beta = TVPparams$beta
  gamma = TVPparams$gamma
  
  y = rep(0,n)
  X = matrix(rnorm(n*p),n,p)
  
  if (is.null(nu2)) {
    for (t in 1:n) {
      y[t] = rnorm(1,sum(X[t,]*gamma[t,]*beta[t,]),sqrt(s2))
    }
    s2out = s2
  }
  
  if (!is.null(nu2)) {
    h = rep(1,n+1)
    m = log(s2)
    h[1] = rnorm(1,m,0.1)
    
    for (t in 1:n) {
      h[t+1] = rnorm(1,m+0.99*(h[t]-m),sqrt(nu2))
      y[t] = rnorm(1,sum(X[t,]*gamma[t,]*beta[t,]),sqrt(exp(h[t+1])))
    }
    s2out = exp(h[-1])
  }
  
  out = list(y=y,X=X,beta=t(beta),gamma=t(gamma),s2=s2out)
  
  out
}


