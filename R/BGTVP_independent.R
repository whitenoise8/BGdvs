#--------  Time varying regression with dependent sparsity ----------#

BGTVP_independent = function(y,X,hyper,options=list(sv=1,smooth=0),initialize="null",
                 Tol_Par = 0.01,Trace = 0){
  
  opt_psi = function(par0,s,mu_q_nu2inv,Q,
                     tol = 1e-2, 
                     maxIt = 100, 
                     maxf = 10) {
    conv = 1
    it = 0
    
    n = length(s)
    p = n+1
    i = rep(1,n)
    dQ = diag(Q)
    
    fold = par0$fold
    s2old = par0$s2old
    Sold = par0$Sold
    
    p0 = fold
    
    while (conv == 1) {
      it = it + 1
      
      ups = -1/2*c(0,i) +1/2*c(0,s)*exp(-fold+1/2*s2old) -mu_q_nu2inv*(dQ*fold-c(0,fold[-p])-c(fold[-1],0))
      invSnew = -1/2*diag(as.vector(c(0,s)*exp(-fold+1/2*s2old))) -mu_q_nu2inv*Q
      Snew = -invtridiag(invSnew)
      
      fnew = fold + Snew%*%ups
      fnew[abs(fnew)>maxf] = sign(fnew[abs(fnew)>maxf])*maxf
      
      par = fnew
      
      fold = fnew
      Sold = Snew
      s2old = diag(Snew)
      
      delta = max(abs((par-p0)/p0))
      if (it > 1) if (delta < tol) conv = 0
      p0 = par
    }
    
    list(f_opt = fnew,
         sigma2_opt = diag(Snew),
         Sigma_opt = Snew,
         it = it,
         convergence = conv)
  }
  
  BGTVPSV = function(y,X,hyper,Tol_Par,Trace,initialize){
    
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
    
    mu_q_h = rep(0,n+1)
    mu_q_s2 = rep(1,n+1)
    mu_q_s2inv = rep(1/2,n+1)
    sigma2_q_h = rep(1,n+1)
    Sigma_q_h = diag(1,n+1)
    
    mu_q_nu2inv = 3
    mu_q_nu2 = 0.3
    
    mu_q_beta = B0
    sigma2_q_beta = matrix(0.8,n+1,p)
    Sigma_q_beta = array(0.8,dim=c(n+1,n+1,p))
    
    mu_q_eta2inv = rep(sqrt(n),p)
    mu_q_eta2 = rep(0.1,p)
    
    mu_q_gamma = matrix(0.5,n,p)

    mu_q_om = matrix(0,n+1,p)
    sigma2_q_om = matrix(0.5,n+1,p)
    Sigma_q_om = array(0.5,dim=c(n+1,n+1,p))
    
    mu_q_z = matrix(1e-3,n,p)
    
    mu_q_xi2inv = rep(0.1,p)
    mu_q_xi2 = rep(1,p)
    
    Q = getQ(n+1,k0)
    dQ = diag(Q)
    
    parOld = c(mu_q_nu2,mu_q_eta2,mu_q_xi2)
    
    p_in = 1:p
    conv = 0
    it = 0
    bound = 10
    
    while (conv==0) {
      it = it + 1
      
      for (j in p_in) {
        
        Om_q_beta = diag(c(0,mu_q_s2inv[-1]*mu_q_gamma[,j]*X2[,j])) + mu_q_eta2inv[j]*Q
        Sigma_q_beta[,,j] = invtridiag(Om_q_beta)
        sigma2_q_beta[,j] = diag(Sigma_q_beta[,,j])
        
        mu_q_eps_j = y-rowSums(X[,-j]*mu_q_gamma[,-j]*mu_q_beta[-1,-j])
        mu_q_beta[,j] = Sigma_q_beta[,,j]%*%c(0,mu_q_s2inv[-1]*mu_q_gamma[,j]*X[,j]*mu_q_eps_j) 
        
        mu_q_betasq = sigma2_q_beta[,j] + mu_q_beta[,j]^2
        
        A_q_eta = Ae + 0.5*(n+1)
        B_q_eta = Be + 0.5*(sum(mu_q_beta[,j]^2*dQ)-2*sum(mu_q_beta[1:n,j]*mu_q_beta[2:(n+1),j]) + 
                              sum(sigma2_q_beta[,j]*dQ)-2*sum(diag(Sigma_q_beta[-(n+1),-1,j]))) 
        
        mu_q_eta2inv[j] = as.numeric(A_q_eta/B_q_eta)
        mu_q_eta2[j] = as.numeric(B_q_eta/(A_q_eta-1))
          
          Om_q_om = diag(c(0,mu_q_z[,j])) + mu_q_xi2inv[j]*Q
          Sigma_q_om[,,j] = invtridiag(Om_q_om)
          sigma2_q_om[,j] = diag(Sigma_q_om[,,j])
          
          mu_q_om[,j] = Sigma_q_om[,,j]%*%c(0,mu_q_gamma[,j]-1/2)
          mu_q_omsq = sigma2_q_om[,j] + mu_q_om[,j]^2
          
          mu_q_z[,j] = 0.5/sqrt(mu_q_omsq[-1])*tanh(0.5*sqrt(mu_q_omsq[-1]))
          
          A_q_xi = Ax + 0.5*(n+1)
          B_q_xi = Bx + 0.5*(sum(mu_q_om[,j]^2*dQ)-2*sum(mu_q_om[1:n,j]*mu_q_om[2:(n+1),j]) + 
                               sum(sigma2_q_om[,j]*dQ)-2*sum(diag(Sigma_q_om[-(n+1),-1,j])))
          mu_q_xi2inv[j] = as.numeric(A_q_xi/B_q_xi)
          mu_q_xi2[j] = as.numeric(B_q_xi/(A_q_xi-1))
         
          omega_q = mu_q_om[-1,j]-0.5*(mu_q_s2inv[-1]*mu_q_betasq[-1]*X2[,j]-2*mu_q_s2inv[-1]*mu_q_beta[-1,j]*X[,j]*mu_q_eps_j)
          omega_q[omega_q > bound] = bound
          omega_q[omega_q < (-bound)] = -bound
          
          mu_q_gamma[,j] = exp(omega_q)/(1+exp(omega_q)) 
        
      }
      
      mu_q_s = y^2 -
        2*rowSums(mu_q_beta[-1,]*mu_q_gamma*X)*y +
        rowSums((mu_q_beta[-1,]^2+sigma2_q_beta[-1,])*mu_q_gamma*X2) +
        rowSums(mu_q_beta[-1,]*mu_q_gamma*X)^2 -
        rowSums((mu_q_beta[-1,]*mu_q_gamma*X)^2)
      
      par0 = list(fold = mu_q_h,
                  s2old = sigma2_q_h,
                  Sold = Sigma_q_h)
      opt_out = opt_psi(par0,mu_q_s,mu_q_nu2inv,Q)
      
      mu_q_h = opt_out$f_opt
      sigma2_q_h = opt_out$sigma2_opt
      Sigma_q_h = opt_out$Sigma_opt
      mu_q_s2 = exp(mu_q_h+1/2*sigma2_q_h)
      mu_q_s2inv = exp(-mu_q_h+1/2*sigma2_q_h)
      
      A_q_nu = Anu + 0.5*(n+1)
      B_q_nu = Bnu + 0.5*(sum(mu_q_h^2*dQ)-2*sum(mu_q_h[1:n]*mu_q_h[2:(n+1)]) + 
                            sum(sigma2_q_h*dQ)-2*sum(diag(Sigma_q_h[-(n+1),-1])))
      mu_q_nu2inv = as.numeric(A_q_nu/B_q_nu)
      mu_q_nu2 = as.numeric(B_q_nu/(A_q_nu-1))
      
      parNew = c(mu_q_nu2,mu_q_eta2,mu_q_xi2)
      
      if (it > 1) {
        delta = max(abs((parNew-parOld))/parOld)
        if (delta < Tol_Par) conv = 1;
        
        p_in = (1:p)[!apply(mu_q_gamma,2,function(x) all(x<0.1))]
        p_out = (1:p)[apply(mu_q_gamma,2,function(x) all(x<0.1))]
        mu_q_gamma[,p_out] = 0
        mu_q_beta[,p_out] = 0
        
        if (Trace == 1) print(paste0("iter: ",it," - error: ",delta," - active: ",length(p_in)))
      }
      
      parOld = parNew
    }
    
    mu_q_beta = t(mu_q_beta)
    mu_q_gamma = t(mu_q_gamma)
    sigma2_q_beta = t(sigma2_q_beta)
    
    outList = list(mu_q_b = mu_q_beta[,-1],
                   sigma2_q_b = sigma2_q_beta[,-1],
                   mu_q_gamma = mu_q_gamma,
                   mu_q_s2 = mu_q_s2[-1],
                   mu_q_h = mu_q_h[-1],
                   Sigma_q_h = Sigma_q_h[-1,-1],
                   mu_q_om = mu_q_om[-1,],
                   mu_q_eta2 = mu_q_eta2,
                   mu_q_xi2 = mu_q_xi2)
    
    outList
  }
  
  BGTVPsm = function(y,X,hyper,Tol_Par,Trace,initialize){
    
    psi = function(f,W,om) {
      xf = as.vector(W%*%f)
      sum((om-xf)*plogis(xf)+log(1+exp(xf)))
    }
    
    grad_psi = function(f,W,om) {
      xf = as.vector(W%*%f)
      colSums((om-xf)*plogis(xf)/(1+exp(xf))*W)
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
    
    smIt = hyper$smIt
    if (is.null(smIt)) smIt = 50
    W = hyper$W
    f_old = matrix(0,ncol(W),p)
    
    mu_q_h = rep(0,n+1)
    mu_q_s2 = rep(1,n+1)
    mu_q_s2inv = rep(0.5,n+1)
    sigma2_q_h = rep(1,n+1)
    Sigma_q_h = diag(1,n+1)
    
    mu_q_nu2inv = 3
    mu_q_nu2 = 0.3
    
    mu_q_beta = B0
    sigma2_q_beta = matrix(0.8,n+1,p)
    Sigma_q_beta = array(0.8,dim=c(n+1,n+1,p))
    
    mu_q_eta2inv = rep(sqrt(n),p)
    mu_q_eta2 = rep(0.1,p)
    
    mu_q_gamma = matrix(0.5,n,p)

    mu_q_om = matrix(0,n+1,p)
    sigma2_q_om = matrix(0.5,n+1,p)
    Sigma_q_om = array(0.5,dim=c(n+1,n+1,p))
    
    mu_q_z = matrix(1e-3,n,p)
    
    mu_q_xi2inv = rep(0.1,p)
    mu_q_xi2 = rep(1,p)
    
    Q = getQ(n+1,k0)
    dQ = diag(Q)
    
    parOld = c(mu_q_nu2,mu_q_eta2,mu_q_xi2)
    
    p_in = 1:p
    conv = 0
    it = 0
    bound = 10
    
    while (conv==0) {
      it = it + 1
      
      for (j in p_in) {
        
        Om_q_beta = diag(c(0,mu_q_s2inv[-1]*mu_q_gamma[,j]*X2[,j])) + mu_q_eta2inv[j]*Q
        Sigma_q_beta[,,j] = invtridiag(Om_q_beta)
        sigma2_q_beta[,j] = diag(Sigma_q_beta[,,j])
        
        mu_q_eps_j = y-rowSums(X[,-j]*mu_q_gamma[,-j]*mu_q_beta[-1,-j])
        mu_q_beta[,j] = Sigma_q_beta[,,j]%*%c(0,mu_q_s2inv[-1]*mu_q_gamma[,j]*X[,j]*mu_q_eps_j) 
        
        mu_q_betasq = sigma2_q_beta[,j] + mu_q_beta[,j]^2
        
        A_q_eta = Ae + 0.5*(n+1)
        B_q_eta = Be + 0.5*(sum(mu_q_beta[,j]^2*dQ)-2*sum(mu_q_beta[1:n,j]*mu_q_beta[2:(n+1),j]) + 
                              sum(sigma2_q_beta[,j]*dQ)-2*sum(diag(Sigma_q_beta[-(n+1),-1,j]))) 
        
        mu_q_eta2inv[j] = as.numeric(A_q_eta/B_q_eta)
        mu_q_eta2[j] = as.numeric(B_q_eta/(A_q_eta-1))
        
          Om_q_om = diag(c(0,mu_q_z[,j])) + mu_q_xi2inv[j]*Q
          Sigma_q_om[,,j] = invtridiag(Om_q_om)
          sigma2_q_om[,j] = diag(Sigma_q_om[,,j])
          
          mu_q_om[,j] = Sigma_q_om[,,j]%*%c(0,mu_q_gamma[,j]-1/2)
          mu_q_omsq = sigma2_q_om[,j] + mu_q_om[,j]^2
          
          mu_q_z[,j] = 0.5/sqrt(mu_q_omsq[-1])*tanh(0.5*sqrt(mu_q_omsq[-1]))
          
          A_q_xi = Ax + 0.5*(n+1)
          B_q_xi = Bx + 0.5*(sum(mu_q_om[,j]^2*dQ)-2*sum(mu_q_om[1:n,j]*mu_q_om[2:(n+1),j]) + 
                               sum(sigma2_q_om[,j]*dQ)-2*sum(diag(Sigma_q_om[-(n+1),-1,j])))
          mu_q_xi2inv[j] = as.numeric(A_q_xi/B_q_xi)
          mu_q_xi2[j] = as.numeric(B_q_xi/(A_q_xi-1))
          
          omega_q = mu_q_om[-1,j]-0.5*(mu_q_s2inv[-1]*mu_q_betasq[-1]*X2[,j]-2*mu_q_s2inv[-1]*mu_q_beta[-1,j]*X[,j]*mu_q_eps_j)
          omega_q[omega_q > bound] = bound
          omega_q[omega_q < (-bound)] = -bound
          
          mu_q_gamma[,j] = plogis(omega_q) 
          
          if (it > smIt) {
            f_new = optim(f_old[,j],function(f) -psi(f,W,omega_q),method="BFGS",
                          gr=function(f) -grad_psi(f,W,omega_q))$par
            mu_q_gamma[,j] = plogis(W%*%f_new) 
            f_old[,j] = f_new
          }
        
      }
      
      mu_q_s = y^2 -
        2*rowSums(mu_q_beta[-1,]*mu_q_gamma*X)*y +
        rowSums((mu_q_beta[-1,]^2+sigma2_q_beta[-1,])*mu_q_gamma*X2) +
        rowSums(mu_q_beta[-1,]*mu_q_gamma*X)^2 -
        rowSums((mu_q_beta[-1,]*mu_q_gamma*X)^2)
      
      par0 = list(fold = mu_q_h,
                  s2old = sigma2_q_h,
                  Sold = Sigma_q_h)
      opt_out = opt_psi(par0,mu_q_s,mu_q_nu2inv,Q)
      
      mu_q_h = opt_out$f_opt
      sigma2_q_h = opt_out$sigma2_opt
      Sigma_q_h = opt_out$Sigma_opt
      mu_q_s2 = exp(mu_q_h+1/2*sigma2_q_h)
      mu_q_s2inv = exp(-mu_q_h+1/2*sigma2_q_h)
      
      A_q_nu = Anu + 0.5*(n+1)
      B_q_nu = Bnu + 0.5*(sum(mu_q_h^2*dQ)-2*sum(mu_q_h[1:n]*mu_q_h[2:(n+1)]) + 
                            sum(sigma2_q_h*dQ)-2*sum(diag(Sigma_q_h[-(n+1),-1])))
      mu_q_nu2inv = as.numeric(A_q_nu/B_q_nu)
      mu_q_nu2 = as.numeric(B_q_nu/(A_q_nu-1))
      
      parNew = c(mu_q_nu2,mu_q_eta2,mu_q_xi2)
      
      if (it > 1) {
        delta = max(abs((parNew-parOld))/parOld)
        if (delta < Tol_Par) conv = 1;
        
        p_in = (1:p)[!apply(mu_q_gamma,2,function(x) all(x<0.1))]
        p_out = (1:p)[apply(mu_q_gamma,2,function(x) all(x<0.1))]
        mu_q_gamma[,p_out] = 0
        mu_q_beta[,p_out] = 0
        
        if (Trace == 1) print(paste0("iter: ",it," - error: ",delta," - active: ",length(p_in)))
      }
      
      parOld = parNew
    }
    
    mu_q_beta = t(mu_q_beta)
    mu_q_gamma = t(mu_q_gamma)
    sigma2_q_beta = t(sigma2_q_beta)
    
    outList = list(mu_q_b = mu_q_beta[,-1],
                   sigma2_q_b = sigma2_q_beta[,-1],
                   mu_q_gamma = mu_q_gamma,
                   mu_q_s2 = mu_q_s2[-1],
                   mu_q_h = mu_q_h[-1],
                   Sigma_q_h = Sigma_q_h[-1,-1],
                   mu_q_om = mu_q_om[-1,],
                   mu_q_eta2 = mu_q_eta2,
                   mu_q_xi2 = mu_q_xi2)
    
    outList
  }
  
  BGTVPhomo = function(y,X,hyper,Tol_Par,Trace,initialize){
    
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
    
    mu_q_s2inv = 1/2
    mu_q_s2 = 2
    
    mu_q_beta = B0
    sigma2_q_beta = matrix(0.8,n+1,p)
    Sigma_q_beta = array(0.8,dim=c(n+1,n+1,p))
    
    mu_q_eta2inv = rep(sqrt(n),p)
    mu_q_eta2 = rep(0.1,p)
    
    mu_q_gamma = matrix(0.5,n,p)

    mu_q_om = matrix(0,n+1,p)
    sigma2_q_om = matrix(0.8,n+1,p)
    Sigma_q_om = array(0.8,dim=c(n+1,n+1,p))
    
    mu_q_z = matrix(1e-3,n,p)
    
    mu_q_xi2inv = rep(0.1,p)
    mu_q_xi2 = rep(1,p)
    
    Q = getQ(n+1,k0)
    dQ = diag(Q)
    
    parOld = c(mu_q_s2,mu_q_eta2,mu_q_xi2)
    
    p_in = 1:p
    conv = 0
    it = 0
    bound = 10
    
    while (conv==0) {
      it = it + 1
      
      for (j in p_in) {
        
        Om_q_beta = mu_q_s2inv*diag(c(0,mu_q_gamma[,j])*c(0,X2[,j])) + mu_q_eta2inv[j]*Q
        Sigma_q_beta[,,j] = invtridiag(Om_q_beta)
        sigma2_q_beta[,j] = diag(Sigma_q_beta[,,j])
        
        mu_q_eps_j = y-rowSums(X[,-j]*mu_q_gamma[,-j]*mu_q_beta[-1,-j])
        mu_q_beta[,j] = Sigma_q_beta[,,j]%*%(mu_q_s2inv*c(0,mu_q_gamma[,j]*X[,j]*mu_q_eps_j))  
        
        mu_q_betasq = sigma2_q_beta[,j] + mu_q_beta[,j]^2
        
        A_q_eta = Ae + 0.5*(n+1)
        B_q_eta = Be + 0.5*(sum(mu_q_beta[,j]^2*dQ)-2*sum(mu_q_beta[1:n,j]*mu_q_beta[2:(n+1),j]) + 
                              sum(sigma2_q_beta[,j]*dQ)-2*sum(diag(Sigma_q_beta[-(n+1),-1,j]))) 
        
        mu_q_eta2inv[j] = as.numeric(A_q_eta/B_q_eta)
        mu_q_eta2[j] = as.numeric(B_q_eta/(A_q_eta-1))
          
          Om_q_om = diag(c(0,mu_q_z[,j])) + mu_q_xi2inv[j]*Q
          Sigma_q_om[,,j] = invtridiag(Om_q_om)
          sigma2_q_om[,j] = diag(Sigma_q_om[,,j])
          
          mu_q_om[,j] = Sigma_q_om[,,j]%*%c(0,mu_q_gamma[,j]-1/2)
          mu_q_omsq = sigma2_q_om[,j] + mu_q_om[,j]^2
          
          mu_q_z[,j] = 0.5/sqrt(mu_q_omsq[-1])*tanh(0.5*sqrt(mu_q_omsq[-1]))
          
          A_q_xi = Ax + 0.5*(n+1)
          B_q_xi = Bx + 0.5*(sum(mu_q_om[,j]^2*dQ)-2*sum(mu_q_om[1:n,j]*mu_q_om[2:(n+1),j]) + 
                               sum(sigma2_q_om[,j]*dQ)-2*sum(diag(Sigma_q_om[-(n+1),-1,j])))
          mu_q_xi2inv[j] = as.numeric(A_q_xi/B_q_xi)
          mu_q_xi2[j] = as.numeric(B_q_xi/(A_q_xi-1))
          
          omega_q = mu_q_om[-1,j]-0.5*mu_q_s2inv*(mu_q_betasq[-1]*X2[,j]-2*mu_q_beta[-1,j]*X[,j]*mu_q_eps_j)
          omega_q[omega_q > bound] = bound
          omega_q[omega_q < (-bound)] = -bound
          
          mu_q_gamma[,j] = exp(omega_q)/(1+exp(omega_q))
          
      }
      
      A_q_s2 = As + 0.5*n
      B_q_s2 = Bs + 0.5*(sum(y^2) -
                           2*sum(rowSums(mu_q_beta[-1,]*mu_q_gamma*X)*y) +
                           sum((mu_q_beta[-1,]^2+sigma2_q_beta[-1,])*mu_q_gamma*X2) +
                           sum(rowSums(mu_q_beta[-1,]*mu_q_gamma*X)^2) -
                           sum(colSums((mu_q_beta[-1,]*mu_q_gamma*X)^2)))
      
      
      mu_q_s2inv = as.numeric(A_q_s2/B_q_s2)
      mu_q_s2 = as.numeric(B_q_s2/(A_q_s2-1))
      
      parNew = c(mu_q_s2,mu_q_eta2,mu_q_xi2)
      
      if (it > 1) {
        delta = max(abs((parNew-parOld))/parOld)
        if (delta < Tol_Par) conv = 1;
        
        p_in = (1:p)[!apply(mu_q_gamma,2,function(x) all(x<0.1))]
        p_out = (1:p)[apply(mu_q_gamma,2,function(x) all(x<0.1))]
        mu_q_gamma[,p_out] = 0
        mu_q_beta[,p_out] = 0
        
        if (Trace == 1) print(paste0("iter: ",it," - error: ",delta," - active: ",length(p_in)))
      }
      
      parOld = parNew
    }
    
    mu_q_beta = t(mu_q_beta)
    mu_q_gamma = t(mu_q_gamma)
    sigma2_q_beta = t(sigma2_q_beta)
    
    outList = list(mu_q_b = mu_q_beta[,-1],
                   sigma2_q_b = sigma2_q_beta[,-1],
                   mu_q_gamma = mu_q_gamma,
                   mu_q_s2 = mu_q_s2,
                   mu_q_om = mu_q_om[-1,],
                   mu_q_eta2 = mu_q_eta2,
                   mu_q_xi2 = mu_q_xi2)
    
    outList
  }
  
  BGTVPhomosm = function(y,X,hyper,Tol_Par,Trace,initialize){
    
    psi = function(f,W,om) {
      xf = as.vector(W%*%f)
      sum((om-xf)*plogis(xf)+log(1+exp(xf)))
    }
    
    grad_psi = function(f,W,om) {
      xf = as.vector(W%*%f)
      colSums((om-xf)*plogis(xf)/(1+exp(xf))*W)
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
    
    smIt = hyper$smIt
    if (is.null(smIt)) smIt = 50
    W = hyper$W
    f_old = matrix(0,ncol(W),p)
    
    mu_q_s2inv = 1
    mu_q_s2 = 1
    
    mu_q_beta = B0
    sigma2_q_beta = matrix(0.8,n+1,p)
    Sigma_q_beta = array(0.8,dim=c(n+1,n+1,p))
    
    mu_q_eta2inv = rep(sqrt(n),p)
    mu_q_eta2 = rep(0.1,p)
    
    mu_q_gamma = matrix(0.5,n,p)

    mu_q_om = matrix(0,n+1,p)
    sigma2_q_om = matrix(0.8,n+1,p)
    Sigma_q_om = array(0.8,dim=c(n+1,n+1,p))
    
    mu_q_z = matrix(1e-3,n,p)
    
    mu_q_xi2inv = rep(0.1,p)
    mu_q_xi2 = rep(1,p)
    
    Q = getQ(n+1,k0)
    dQ = diag(Q)
    
    parOld = c(mu_q_s2,mu_q_eta2,mu_q_xi2)
    
    p_in = 1:p
    conv = 0
    it = 0
    bound = 10
    
    while (conv==0) {
      it = it + 1
      
      for (j in p_in) {
        
        Om_q_beta = mu_q_s2inv*diag(c(0,mu_q_gamma[,j])*c(0,X2[,j])) + mu_q_eta2inv[j]*Q
        Sigma_q_beta[,,j] = invtridiag(Om_q_beta)
        sigma2_q_beta[,j] = diag(Sigma_q_beta[,,j])
        
        mu_q_eps_j = y-rowSums(X[,-j]*mu_q_gamma[,-j]*mu_q_beta[-1,-j])
        mu_q_beta[,j] = Sigma_q_beta[,,j]%*%(mu_q_s2inv*c(0,mu_q_gamma[,j]*X[,j]*mu_q_eps_j))  
        
        mu_q_betasq = sigma2_q_beta[,j] + mu_q_beta[,j]^2
        
        A_q_eta = Ae + 0.5*(n+1)
        B_q_eta = Be + 0.5*(sum(mu_q_beta[,j]^2*dQ)-2*sum(mu_q_beta[1:n,j]*mu_q_beta[2:(n+1),j]) + 
                              sum(sigma2_q_beta[,j]*dQ)-2*sum(diag(Sigma_q_beta[-(n+1),-1,j]))) 
        
        mu_q_eta2inv[j] = as.numeric(A_q_eta/B_q_eta)
        mu_q_eta2[j] = as.numeric(B_q_eta/(A_q_eta-1))
          
          Om_q_om = diag(c(0,mu_q_z[,j])) + mu_q_xi2inv[j]*Q
          Sigma_q_om[,,j] = invtridiag(Om_q_om)
          sigma2_q_om[,j] = diag(Sigma_q_om[,,j])
          
          mu_q_om[,j] = Sigma_q_om[,,j]%*%c(0,mu_q_gamma[,j]-1/2)
          mu_q_omsq = sigma2_q_om[,j] + mu_q_om[,j]^2
          
          mu_q_z[,j] = 0.5/sqrt(mu_q_omsq[-1])*tanh(0.5*sqrt(mu_q_omsq[-1]))
          
          A_q_xi = Ax + 0.5*(n+1)
          B_q_xi = Bx + 0.5*(sum(mu_q_om[,j]^2*dQ)-2*sum(mu_q_om[1:n,j]*mu_q_om[2:(n+1),j]) + 
                               sum(sigma2_q_om[,j]*dQ)-2*sum(diag(Sigma_q_om[-(n+1),-1,j])))
          mu_q_xi2inv[j] = as.numeric(A_q_xi/B_q_xi)
          mu_q_xi2[j] = as.numeric(B_q_xi/(A_q_xi-1))
          
          omega_q = mu_q_om[-1,j]-0.5*mu_q_s2inv*(mu_q_betasq[-1]*X2[,j]-2*mu_q_beta[-1,j]*X[,j]*mu_q_eps_j)
          omega_q[omega_q > bound] = bound
          omega_q[omega_q < (-bound)] = -bound
          
          mu_q_gamma[,j] = exp(omega_q)/(1+exp(omega_q))
          
          if (it > smIt) {
            f_new = optim(f_old[,j],function(f) -psi(f,W,omega_q),method="BFGS",
                          gr=function(f) -grad_psi(f,W,omega_q))$par
            mu_q_gamma[,j] = plogis(W%*%f_new) 
            f_old[,j] = f_new
          }
        
      }
      
      A_q_s2 = As + 0.5*n
      B_q_s2 = Bs + 0.5*(sum(y^2) -
                           2*sum(rowSums(mu_q_beta[-1,]*mu_q_gamma*X)*y) +
                           sum((mu_q_beta[-1,]^2+sigma2_q_beta[-1,])*mu_q_gamma*X2) +
                           sum(rowSums(mu_q_beta[-1,]*mu_q_gamma*X)^2) -
                           sum(colSums((mu_q_beta[-1,]*mu_q_gamma*X)^2)))
      
      
      mu_q_s2inv = as.numeric(A_q_s2/B_q_s2)
      mu_q_s2 = as.numeric(B_q_s2/(A_q_s2-1))
      
      parNew = c(mu_q_s2,mu_q_eta2,mu_q_xi2)
      
      if (it > 1) {
        delta = max(abs((parNew-parOld))/parOld)
        if (delta < Tol_Par) conv = 1;
        
        p_in = (1:p)[!apply(mu_q_gamma,2,function(x) all(x<0.1))]
        p_out = (1:p)[apply(mu_q_gamma,2,function(x) all(x<0.1))]
        mu_q_gamma[,p_out] = 0
        mu_q_beta[,p_out] = 0
        
        if (Trace == 1) print(paste0("iter: ",it," - error: ",delta," - active: ",length(p_in)))
      }
      
      parOld = parNew
    }
    
    mu_q_beta = t(mu_q_beta)
    mu_q_gamma = t(mu_q_gamma)
    sigma2_q_beta = t(sigma2_q_beta)
    
    outList = list(mu_q_b = mu_q_beta[,-1],
                   sigma2_q_b = sigma2_q_beta[,-1],
                   mu_q_gamma = mu_q_gamma,
                   mu_q_s2 = mu_q_s2,
                   mu_q_om = mu_q_om[,-1],
                   mu_q_eta2 = mu_q_eta2,
                   mu_q_xi2 = mu_q_xi2)
    
    outList
  }
  
  outMod = NULL
  
  if ((options$sv==1)&(options$smooth==0)) {
    outMod = BGTVPSV(y,X,hyper,Tol_Par,Trace,initialize)
  }
  
  if ((options$sv==0)&(options$smooth==0)) {
    outMod = BGTVPhomo(y,X,hyper,Tol_Par,Trace,initialize)
  }
  
  if ((options$sv==0)&(options$smooth==1)) {
    outMod = BGTVPhomosm(y,X,hyper,Tol_Par,Trace,initialize)
  }
  
  if ((options$sv==1)&(options$smooth==1)) {
    outMod = BGTVPsm(y,X,hyper,Tol_Par,Trace,initialize)
  }
  
  if (is.null(outMod)) print("Options setting not available!")
  
  outMod
}
