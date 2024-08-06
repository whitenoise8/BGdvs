#--------  Time varying regression with dependent sparsity ----------#

BGTVP_dependent = function(y,X,hyper,options=list(sv=1,smooth=0),initialize='null',
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
  
  BGTVPSV_joint = function(y,X,hyper,Tol_Par,Trace,initialize){
    
    #: Code
    n = length(y)
    p = ncol(X)
    
    B0 = matrix(0,n+1,p,byrow=TRUE)
    
    if (initialize == "OLS") {
      B0 = solve(t(X)%*%X+diag(1,p))%*%t(X)%*%y
      B0 = matrix(B0,n+1,p,byrow=TRUE)
    }
    
    XLs = list()
    for (t in 1:nrow(X)) XLs[[t]] = matrix(X[t,],nrow=1)
    Xbd = bdiag(XLs)
    Xtl = bdiag(matrix(0,1,p),Xbd)
    Xtl2 = t(Xtl)%*%Xtl
    
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
    mu_q_beta_vec = as(t(mu_q_beta),"vector")
    sigma2_q_beta = matrix(0.8,n+1,p)

    mu_q_eta2inv = rep(sqrt(n),p)
    mu_q_eta2 = rep(0.1,p)
    
    mu_q_gamma = matrix(0.5,n,p)
    mu_q_gamma_vec = c(rep(1,p),as(t(mu_q_gamma),"vector"))
    
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
    p_all = p_in
    p_ini = p
    dropped = NULL
    conv = 0
    it = 0
    bound = 10

    while (conv==0) {
      it = it + 1
      

      Om_q_beta = bdiag(getXX(matrix(mu_q_gamma,ncol=p),matrix(X,ncol=p),mu_q_s2inv)) + kronecker(Q,Diagonal(x=mu_q_eta2inv))
      
      if (p > 1) Sigma_elem = getQuantities(as(Om_q_beta,"matrix"),n+1,p)
      if (p == 1) {
        Sigma_full = invtridiag(as(Om_q_beta,"matrix"))
        Sigma_elem = list()
        Sigma_elem$D = array(diag(Sigma_full),dim=c(1,1,n+1))
        Sigma_elem$oD = matrix(Sigma_full[row(Sigma_full)==col(Sigma_full)+1],ncol=1)
      }
      
      Sigma_bdiag = Sigma_elem$D
      sigma2_q_beta = t(apply(Sigma_bdiag,3,diag))
      cov_q_beta = Sigma_elem$oD
      
      GX = Diagonal(x=mu_q_gamma_vec)%*%t(Xtl)
      lam = GX%*%c(0,mu_q_s2inv[-1]*y)
      mu_q_beta_vec = solve(Om_q_beta,lam)
      
      mu_q_beta = matrix(mu_q_beta_vec,n+1,p,byrow=TRUE)
      
      mu_q_betasq = sigma2_q_beta + mu_q_beta^2
      
      for (j in p_in) {
        
        A_q_eta = Ae + 0.5*(n+1)
        B_q_eta = Be + 0.5*(sum(mu_q_beta[,j]^2*dQ)-2*sum(mu_q_beta[1:n,j]*mu_q_beta[2:(n+1),j]) + 
                              sum(sigma2_q_beta[,j]*dQ)-2*sum(cov_q_beta[,j])) 
        
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
        
        for (t in 1:n) {
          cov_mat_t = matrix(Sigma_bdiag[,,t+1],p,p)
          omega_q = mu_q_om[t+1,j]-0.5*(mu_q_s2inv[t+1]*mu_q_betasq[t+1,j]*X[t,j]^2 -
                                          2*mu_q_s2inv[t+1]*mu_q_beta[t+1,j]*X[t,j]*y[t] +
                                          2*mu_q_s2inv[t+1]*sum(X[t,j]*X[t,-j]*mu_q_gamma[t,-j]*(mu_q_beta[t+1,j]*mu_q_beta[t+1,-j] + cov_mat_t[j,-j])))
          
          omega_q[omega_q > bound] = bound
          omega_q[omega_q < (-bound)] = -bound
          
          mu_q_gamma[t,j] = exp(omega_q)/(1+exp(omega_q)) 
        }
        
      }
      
      mu_q_gamma_vec = c(rep(1,p),as(t(mu_q_gamma),"vector"))
      
      mu_q_s = y^2 -
        2*rowSums(mu_q_beta[-1,]*mu_q_gamma*X)*y +
        rowSums((mu_q_beta[-1,]^2+sigma2_q_beta[-1,])*mu_q_gamma*X^2) +
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
      mu_q_s2inv = as(exp(-mu_q_h+1/2*sigma2_q_h),"vector")
      
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
        
        p_ = length(p_in)
        if (p_ < p) {
          dropped = c(dropped,p_all[p_out])
          p_all = p_all[-p_out]
          
          X = X[,p_in]
          XLs = list()
          for (t in 1:nrow(X)) XLs[[t]] = matrix(X[t,],nrow=1)
          Xbd = bdiag(XLs)
          Xtl = bdiag(matrix(0,1,p_),Xbd)
          Xtl2 = t(Xtl)%*%Xtl
          
          mu_q_beta = mu_q_beta[,p_in]
          mu_q_beta_vec = as(t(mu_q_beta),"vector")
          sigma2_q_beta = sigma2_q_beta[,p_in]
          
          mu_q_eta2inv = mu_q_eta2inv[p_in]
          mu_q_eta2 = mu_q_eta2[p_in]
          
          mu_q_gamma = mu_q_gamma[,p_in]
          mu_q_gamma_vec = c(rep(1,p_),as(t(mu_q_gamma),"vector"))
          
          mu_q_om = mu_q_om[,p_in]
          sigma2_q_om = sigma2_q_om[,p_in]
          Sigma_q_om = Sigma_q_om[,,p_in]
          
          mu_q_z = mu_q_z[,p_in]
          
          mu_q_xi2inv = mu_q_xi2[p_in]
          mu_q_xi2 = mu_q_xi2[p_in]
          
          p = p_
          p_in = 1:p
          parNew = c(mu_q_nu2,mu_q_eta2,mu_q_xi2)
        }
        if (Trace == 1) print(paste0("iter: ",it," - error: ",delta," - active: ",p))
      }
      
      parOld = parNew
    }
    
    p_in = setdiff(1:p_ini,dropped)
    
    mu_q_beta = t(mu_q_beta)
    Mu_q_beta = matrix(0,p_ini,n+1)
    Mu_q_beta[p_in,] = mu_q_beta
    
    mu_q_gamma = t(mu_q_gamma)
    Mu_q_gamma = matrix(0,p_ini,n)
    Mu_q_gamma[p_in,] = mu_q_gamma
    
    sigma2_q_beta = t(sigma2_q_beta)
    Sigma2_q_beta = matrix(0,p_ini,n+1)
    Sigma2_q_beta[p_in,] = sigma2_q_beta
    
    outList = list(mu_q_beta = Mu_q_beta[,-1],
                   sigma2_q_beta = Sigma2_q_beta[,-1],
                   mu_q_gamma = Mu_q_gamma,
                   mu_q_s2 = mu_q_s2[-1],
                   mu_q_h = mu_q_h[-1],
                   Sigma_q_h = Sigma_q_h[-1,-1],
                   mu_q_om = mu_q_om[-1,],
                   mu_q_eta2 = mu_q_eta2,
                   mu_q_xi2 = mu_q_xi2,
                   p_active = p_in)
    
    outList
  }
  
  BGTVPhomo_joint = function(y,X,hyper,Tol_Par,Trace,initialize){
    
    #: Code
    n = length(y)
    p = ncol(X)
    
    B0 = matrix(0,n+1,p,byrow=TRUE)
    
    if (initialize == "OLS") {
      B0 = solve(t(X)%*%X+diag(1,p))%*%t(X)%*%y
      B0 = matrix(B0,n+1,p,byrow=TRUE)
    }
    
    XLs = list()
    for (t in 1:nrow(X)) XLs[[t]] = matrix(X[t,],nrow=1)
    Xbd = bdiag(XLs)
    Xtl = bdiag(matrix(0,1,p),Xbd)
    Xtl2 = t(Xtl)%*%Xtl
    
    As = hyper$As
    Bs = hyper$Bs
    Ae = hyper$Ae
    Be = hyper$Be
    Ax = hyper$Ax
    Bx = hyper$Bx
    k0 = hyper$k0
    
    mu_q_s2inv = 1
    mu_q_s2 = 1
    
    mu_q_beta = B0
    mu_q_beta_vec = as(t(mu_q_beta),"vector")
    sigma2_q_beta = matrix(0.8,n+1,p)

    mu_q_eta2inv = rep(sqrt(n),p)
    mu_q_eta2 = rep(0.1,p)
    
    mu_q_gamma = matrix(0.5,n,p)
    mu_q_gamma_vec = c(rep(1,p),as(t(mu_q_gamma),"vector"))
    
    mu_q_om = matrix(0,n+1,p)
    sigma2_q_om = matrix(0.5,n+1,p)
    Sigma_q_om = array(0.5,dim=c(n+1,n+1,p))
    
    mu_q_z = matrix(1e-3,n,p)
    
    mu_q_xi2inv = rep(0.1,p)
    mu_q_xi2 = rep(1,p)
    
    Q = getQ(n+1,k0)
    dQ = diag(Q)
    
    parOld = c(mu_q_s2,mu_q_eta2,mu_q_xi2)
    
    p_in = 1:p
    p_all = p_in
    p_ini = p
    dropped = NULL
    conv = 0
    it = 0
    bound = 10
    
    while (conv==0) {
      it = it + 1
      
      G2 = mu_q_gamma_vec%*%t(mu_q_gamma_vec) + diag(mu_q_gamma_vec*(1-mu_q_gamma_vec))
      Om_q_beta = Xtl2 * kronecker(diag(mu_q_s2inv,n+1),matrix(1,p,p)) * G2 + kronecker(Q,diag(mu_q_eta2inv))
      Sigma_q_beta = solve(Om_q_beta)
      
      Sigma_bdiag = matrix(fatdiag(Sigma_q_beta,steps=n+1,size=p),p*p,n+1)
      sigma2_q_beta = matrix(diag(Sigma_q_beta),n+1,p,byrow=TRUE)
      cov_q_beta = matrix(Sigma_q_beta[row(Sigma_q_beta)-col(Sigma_q_beta)==p],n,p,byrow=TRUE)
      
      GX = diag(mu_q_gamma_vec)%*%t(Xtl)
      mu_q_beta_vec = as(Sigma_q_beta%*%GX%*%c(0,mu_q_s2inv*y),"vector")
      mu_q_beta = matrix(mu_q_beta_vec,n+1,p,byrow=TRUE)
      
      mu_q_betasq = sigma2_q_beta + mu_q_beta^2
      
      for (j in p_in) {
        
        A_q_eta = Ae + 0.5*(n+1)
        B_q_eta = Be + 0.5*(sum(mu_q_beta[,j]^2*dQ)-2*sum(mu_q_beta[1:n,j]*mu_q_beta[2:(n+1),j]) + 
                              sum(sigma2_q_beta[,j]*dQ)-2*sum(cov_q_beta[,j])) 
        
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
        
        for (t in 1:n) {
          cov_mat_t = matrix(Sigma_bdiag[,t+1],p,p)
          omega_q = mu_q_om[t+1,j]-0.5*(mu_q_s2inv*mu_q_betasq[t+1,j]*X[t,j]^2 -
                                          2*mu_q_s2inv*mu_q_beta[t+1,j]*X[t,j]*y[t] +
                                          2*mu_q_s2inv*sum(X[t,j]*X[t,-j]*mu_q_gamma[t,-j]*(mu_q_beta[t+1,j]*mu_q_beta[t+1,-j] + cov_mat_t[j,-j])))
          
          omega_q[omega_q > bound] = bound
          omega_q[omega_q < (-bound)] = -bound
          
          mu_q_gamma[t,j] = exp(omega_q)/(1+exp(omega_q)) 
        }
        
      }
      
      mu_q_gamma_vec = c(rep(1,p),as(t(mu_q_gamma),"vector"))
      
      A_q_s2 = As + 0.5*n
      B_q_s2 = Bs + 0.5*(sum(y^2) -
                           2*sum(rowSums(mu_q_beta[-1,]*mu_q_gamma*X)*y) +
                           sum((mu_q_beta[-1,]^2+sigma2_q_beta[-1,])*mu_q_gamma*X^2) +
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

        p_ = length(p_in)
        if (p_ < p) {
          dropped = c(dropped,p_all[p_out])
          p_all = p_all[-p_out]
          
          X = X[,p_in]
          XLs = list()
          for (t in 1:nrow(X)) XLs[[t]] = matrix(X[t,],nrow=1)
          Xbd = bdiag(XLs)
          Xtl = bdiag(matrix(0,1,p_),Xbd)
          Xtl2 = t(Xtl)%*%Xtl
          
          mu_q_beta = mu_q_beta[,p_in]
          mu_q_beta_vec = as(t(mu_q_beta),"vector")
          sigma2_q_beta = sigma2_q_beta[,p_in]
          
          mu_q_eta2inv = mu_q_eta2inv[p_in]
          mu_q_eta2 = mu_q_eta2[p_in]
          
          mu_q_gamma = mu_q_gamma[,p_in]
          mu_q_gamma_vec = c(rep(1,p_),as(t(mu_q_gamma),"vector"))
          
          mu_q_om = mu_q_om[,p_in]
          sigma2_q_om = sigma2_q_om[,p_in]
          Sigma_q_om = Sigma_q_om[,,p_in]
          
          mu_q_z = mu_q_z[,p_in]
          
          mu_q_xi2inv = mu_q_xi2[p_in]
          mu_q_xi2 = mu_q_xi2[p_in]
          
          p = p_
          p_in = 1:p
          parNew = c(mu_q_s2,mu_q_eta2,mu_q_xi2)
        }
        
        if (Trace == 1) print(paste0("iter: ",it," - error: ",delta," - active: ",length(p_in)))
      }
      
      parOld = parNew
    }
    
    p_in = setdiff(1:p_ini,dropped)
    
    mu_q_beta = t(mu_q_beta)
    Mu_q_beta = matrix(0,p_ini,n+1)
    Mu_q_beta[p_in,] = mu_q_beta
    
    mu_q_gamma = t(mu_q_gamma)
    Mu_q_gamma = matrix(0,p_ini,n)
    Mu_q_gamma[p_in,] = mu_q_gamma
    
    sigma2_q_beta = t(sigma2_q_beta)
    Sigma2_q_beta = matrix(0,p_ini,n+1)
    Sigma2_q_beta[p_in,] = sigma2_q_beta
    
    outList = list(mu_q_beta = Mu_q_beta[,-1],
                   sigma2_q_beta = Sigma2_q_beta[,-1],
                   mu_q_gamma = Mu_q_gamma,
                   mu_q_s2 = mu_q_s2,
                   mu_q_om = mu_q_om[-1,],
                   mu_q_eta2 = mu_q_eta2,
                   mu_q_xi2 = mu_q_xi2,
                   p_active = p_in)
    
    outList
  }
  
  
  require(Matrix)
  require(diagonals)
  
  outMod = NULL
  
  if ((options$sv==1)&(options$smooth==0)) {
    outMod = BGTVPSV_joint(y,X,hyper,Tol_Par,Trace,initialize)
  }
  
  if ((options$sv==0)&(options$smooth==0)) {
    outMod = BGTVPhomo_joint(y,X,hyper,Tol_Par,Trace,initialize)
  }
  
  if (is.null(outMod)) print("Options setting not available!")
  
  outMod
}
