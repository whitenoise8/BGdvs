#--------  Time varying regression with dependent sparsity ----------#

BGTVPSV_MCMC = function(y,X,hyper,ndraws,nburn,initialize='null',Trace = 0){
  
  require(stochvol)
  require(pgdraw)
  
  priors_sv = specify_priors(
    mu = sv_normal(0,1e-4),
    phi = sv_constant(0.999),
    sigma2 = sv_gamma(shape = 0.01, rate = 0.01))
  
  rmnorm = function(A, b) {
    R = chol(A)
    z = rnorm(length(b))
    y = drop(forwardsolve(t(R), b, upper.tri = FALSE, transpose = FALSE))
    x = drop(backsolve(R, y + z, upper.tri = TRUE, transpose = FALSE))
    return(x)
  }
  
  #: Code
  n = length(y)
  p = ncol(X)
  
  X2 = X^2
  
  B0 = matrix(0,n+1,p,byrow=TRUE)
  
  if (initialize == "OLS") {
    B0 = solve(t(X)%*%X+diag(1,p))%*%t(X)%*%y
    B0 = matrix(B0,n+1,p,byrow=TRUE)
  }
  
  Anu = hyper$As
  Bnu = hyper$Bs
  Ae = hyper$Ae
  Be = hyper$Be
  Ax = hyper$Ax
  Bx = hyper$Bx
  k0 = hyper$k0
  
  s2_sim = matrix(2,n+1,ndraws)
  nu2_sim = 1
  
  beta_sim = array(0,dim=c(n+1,p,ndraws))
  beta_sim[,,1] = B0
  gamma_sim = array(1,dim=c(n,p,ndraws))
  beta_gamma_sim = array(0,dim=c(n,p,ndraws))
  om_sim = array(0,dim=c(n+1,p,ndraws))
  
  z_sim = matrix(0.1,n,p)
  eta2_sim = rep(1,p)
  xi2_sim = rep(1,p)
  
  Q = getQ(n+1,k0)
  dQ = diag(Q)
  
  bound = 100
  for (r in 2:ndraws) {
    
    beta_sim[,,r] = beta_sim[,,r-1]
    gamma_sim[,,r] = gamma_sim[,,r-1]
    
    for (j in 1:p) {
      
      Om_q_beta = diag(c(0,1/s2_sim[-1,r-1]*gamma_sim[,j,r]*X2[,j])) + 1/eta2_sim[j]*Q
      mu_q_eps_j = y-rowSums(X[,-j]*gamma_sim[,-j,r]*beta_sim[-1,-j,r])
      lam_q_beta = c(0,1/s2_sim[-1,r-1]*gamma_sim[,j,r]*X[,j]*mu_q_eps_j) 
      
      beta_sim[,j,r] = rmnorm(Om_q_beta,lam_q_beta)
      mu_q_betasq = beta_sim[,j,r]^2
      
      omega_q = om_sim[-1,j,r-1]-0.5*(1/s2_sim[-1,r-1]*mu_q_betasq[-1]*X2[,j]-2*1/s2_sim[-1,r-1]*beta_sim[-1,j,r]*X[,j]*mu_q_eps_j)
      omega_q[omega_q > bound] = bound
      omega_q[omega_q < (-bound)] = -bound
      
      gamma_sim[,j,r] = rbinom(n,1,exp(omega_q)/(1+exp(omega_q)))
      
      A_q_eta = Ae + 0.5*(n+1)
      B_q_eta = Be + 0.5*(sum(beta_sim[,j,r]^2*dQ)-2*sum(beta_sim[1:n,j,r]*beta_sim[2:(n+1),j,r])) 
      eta2_sim[j] = 1/rgamma(1,A_q_eta,B_q_eta)
      
      Om_q_om = diag(c(0,z_sim[,j])) + 1/xi2_sim[j]*Q
      lam_q_om = c(0,gamma_sim[,j,r]-1/2)
      om_sim[,j,r] = rmnorm(Om_q_om,lam_q_om)          
      
      z_sim[,j] = pgdraw(1,om_sim[-1,j,r])
      
      A_q_xi = Ax + 0.5*(n+1)
      B_q_xi = Bx + 0.5*(sum(om_sim[,j,r]^2*dQ)-2*sum(om_sim[1:n,j,r]*om_sim[2:(n+1),j,r]))
      xi2_sim[j] = 1/rgamma(1,A_q_xi,B_q_xi)
    }
    beta_gamma_sim[,,r] = beta_sim[-1,,r]*gamma_sim[,,r]
    
    res = y - rowSums(beta_gamma_sim[,,r]*X)
    st_param = list(mu=0,phi=0.999,sigma=sqrt(nu2_sim),latent0=log(s2_sim[1,r-1]))
    
    sv_mod = svsample_fast_cpp(res,priorspec=priors_sv,startpara=st_param,startlatent=log(s2_sim[-1,r-1]))
    s2_sim[-1,r] = exp(as(sv_mod$latent,"numeric"))
    s2_sim[1,r] = exp(as(sv_mod$latent0,"numeric"))
    nu2_sim = as(sv_mod$para[,"sigma"]^2,"numeric")
    
    if (r %in% seq(1000,ndraws,by=1000)) {
      if (Trace == 1) print(paste0("Sampling: ",r," of ",ndraws))
    }
    
  }
  
  beta_sim = beta_sim[-1,,(nburn+1):ndraws]
  gamma_sim = gamma_sim[,,(nburn+1):ndraws]
  beta_gamma_sim = beta_gamma_sim[,,(nburn+1):ndraws]
  s2_sim = s2_sim[-1,(nburn+1):ndraws]
  
  outList = list(b = beta_sim,
                 gamma = gamma_sim,
                 beta = beta_gamma_sim,
                 s2 = s2_sim)
  
  outList
}

BGTVP_MCMC = function(y,X,hyper,ndraws,nburn,initialize='null',Trace = 0){
  
  require(stochvol)
  require(pgdraw)
  
  rmnorm = function(A, b) {
    R = chol(A)
    z = rnorm(length(b))
    y = drop(forwardsolve(t(R), b, upper.tri = FALSE, transpose = FALSE))
    x = drop(backsolve(R, y + z, upper.tri = TRUE, transpose = FALSE))
    return(x)
  }
  
  #: Code
  n = length(y)
  p = ncol(X)
  
  X2 = X^2
  
  B0 = matrix(0,n+1,p,byrow=TRUE)
  
  if (initialize == "OLS") {
    B0 = solve(t(X)%*%X+diag(1,p))%*%t(X)%*%y
    B0 = matrix(B0,n+1,p,byrow=TRUE)
  }
  
  As = hyper$As
  Bs = hyper$Bs
  Ae = hyper$Ae
  Be = hyper$Be
  Ax = hyper$Ax
  Bx = hyper$Bx
  k0 = hyper$k0
  
  s2_sim = rep(1,ndraws)
  beta_sim = array(0,dim=c(n+1,p,ndraws))
  beta_sim[,,1] = B0
  gamma_sim = array(1,dim=c(n,p,ndraws))
  beta_gamma_sim = array(1,dim=c(n,p,ndraws))
  om_sim = array(0,dim=c(n+1,p,ndraws))
  
  z_sim = matrix(0.1,n,p)
  eta2_sim = rep(1,p)
  xi2_sim = rep(1,p)
  
  Q = getQ(n+1,k0)
  dQ = diag(Q)
  
  bound = 100
  for (r in 2:ndraws) {
    
    beta_sim[,,r] = beta_sim[,,r-1]
    gamma_sim[,,r] = gamma_sim[,,r-1]
    
    for (j in 1:p) {
      
      Om_q_beta = diag(c(0,1/s2_sim[r-1]*gamma_sim[,j,r]*X2[,j])) + 1/eta2_sim[j]*Q
      mu_q_eps_j = y-rowSums(X[,-j]*gamma_sim[,-j,r]*beta_sim[-1,-j,r])
      lam_q_beta = c(0,1/s2_sim[r-1]*gamma_sim[,j,r]*X[,j]*mu_q_eps_j) 
      
      beta_sim[,j,r] = rmnorm(Om_q_beta,lam_q_beta)
      mu_q_betasq = beta_sim[,j,r]^2
      
      omega_q = om_sim[-1,j,r-1]-0.5*(1/s2_sim[r-1]*mu_q_betasq[-1]*X2[,j]-2*1/s2_sim[r-1]*beta_sim[-1,j,r]*X[,j]*mu_q_eps_j)
      omega_q[omega_q > bound] = bound
      omega_q[omega_q < (-bound)] = -bound
      
      gamma_sim[,j,r] = rbinom(n,1,exp(omega_q)/(1+exp(omega_q)))

      A_q_eta = Ae + 0.5*(n+1)
      B_q_eta = Be + 0.5*(sum(beta_sim[,j,r]^2*dQ)-2*sum(beta_sim[1:n,j,r]*beta_sim[2:(n+1),j,r])) 
      eta2_sim[j] = 1/rgamma(1,A_q_eta,B_q_eta)
      
      Om_q_om = diag(c(0,z_sim[,j])) + 1/xi2_sim[j]*Q
      lam_q_om = c(0,gamma_sim[,j,r]-1/2)
      om_sim[,j,r] = rmnorm(Om_q_om,lam_q_om)          
      
      z_sim[,j] = pgdraw(1,om_sim[-1,j,r])
      
      A_q_xi = Ax + 0.5*(n+1)
      B_q_xi = Bx + 0.5*(sum(om_sim[,j,r]^2*dQ)-2*sum(om_sim[1:n,j,r]*om_sim[2:(n+1),j,r]))
      xi2_sim[j] = 1/rgamma(1,A_q_xi,B_q_xi)
    }
    beta_gamma_sim[,,r] = beta_sim[-1,,r]*gamma_sim[,,r]
    
    A_q_s2 = As + 0.5*n
    B_q_s2 = Bs + 0.5*sum((y - rowSums(beta_gamma_sim[,,r]*X))^2)
    
    s2_sim[r] = 1/rgamma(1,A_q_s2,B_q_s2)

    if (r %in% seq(1000,ndraws,by=1000)) {
      if (Trace == 1) print(paste0("Sampling: ",r," of ",ndraws))
    }
    
  }
  
  beta_sim = beta_sim[,,(nburn+1):ndraws]
  gamma_sim = gamma_sim[,,(nburn+1):ndraws]
  beta_gamma_sim = beta_gamma_sim[,,(nburn+1):ndraws]
  s2_sim = s2_sim[(nburn+1):ndraws]
  
  outList = list(b = beta_sim,
                 gamma = gamma_sim,
                 beta = beta_gamma_sim,
                 s2 = s2_sim)
  
  outList
}
