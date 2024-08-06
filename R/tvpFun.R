getTVPparams = function(n,p_in,p_dyn,p_low,p_null,eta2,lam,delta,thresh=0.5) {
  
  p = p_in+sum(p_dyn)+p_low+p_null
  beta = NULL
  
  if (p_in > 0) {
    beta1 = matrix(0,n+1,p_in)
    for (j in 1:p_in) {
      c = sample(c(-1,1),1)*runif(1,2,4)
      beta1[1,j] = rnorm(1,c,0.01)
      for (t in 1:n) {
        while (abs(beta1[t+1,j]) < thresh) {
          beta1[t+1,j] = rnorm(1,c+0.98*(beta1[t,j]-c),sqrt(eta2))
        }
      }
    }
    beta = cbind(beta,beta1)
  }
  
  if (sum(p_dyn) > 0) {
    beta_ = matrix(0,n+1,sum(p_dyn))
    
    for (j in 1:sum(p_dyn)) {
      
      tint = cumsum(c(1,rpois(n,lam[j])))
      tint = tint[tint<n]
      tint = c(tint,n)
      
      signal = rbinom(1,1,0.5)
      
      for (id in 2:length(tint)-1) {
        
        if (signal == 0)  {
          beta_[tint[id]:tint[id+1],j] = 0
          signal = 1
        } else {
          c = sample(c(-1,1),1)*runif(1,thresh,1.5)
          beta_[tint[id],j] = c
          for (t in tint[id]:tint[id+1]) {
            while (abs(beta_[t+1,j]) < thresh) {
              beta_[t+1,j] = rnorm(1,c+0.98*(beta_[t,j]-c),sqrt(eta2))
            }
          }
          signal = 0
        }
      }
    }
    
    beta = cbind(beta,beta_)
  }
  
  if (p_low > 0) {
    beta1_ = matrix(0,n+1,p_low)
    
    for (j in 1:p_low) {
      
      tint = rpois(1,delta[j])
      start = sample(1:(n-tint),1)
      
      c = sample(c(-1,1),1)*runif(1,thresh,1.5)
      beta1_[start,j] = c
      for (t in start:(start+tint)) {
        while (abs(beta1_[t+1,j]) < thresh) {
          beta1_[t+1,j] = rnorm(1,c+0.98*(beta1_[t,j]-c),sqrt(eta2))
        }
      }
      
    }
    
    beta = cbind(beta,beta1_)
  }
  
  if (p_null > 0) {
    beta0 = matrix(0,n+1,p_null)
    beta = cbind(beta,beta0)
  }
  
  beta = matrix(beta,n+1,p)
  beta = beta[-1,]
  beta = matrix(beta,n,p)
  
  gamma = matrix(1,n,p)
  gamma[beta == 0] = 0
  
  out = list(beta=beta,gamma=gamma)
  
  out
}

getQ = function(n,k0) {
  m = diag(1,n)
  diag(m)[-c(1,n)] = 2
  diag(m[-nrow(m),-1]) = -1
  diag(m[-1,-ncol(m)]) = -1
  m[1,1] = 1+1/k0
  m
}

matrixplot = function(m1, ynames = NULL, xdate = NULL) {
  require(ggplot2)
  require(reshape2)
  
    p = nrow(m1)
    n = ncol(m1)
    
    if (is.null(ynames)) ynames = 1:p
    if (is.null(xdate)) xdate = 1:n
    
    m.mat = data.frame(Var1=rep(xdate,p),
                    Var2=rep(1:p,each=n),
                    Value=as.vector(t(m1)[,p:1]))
    
    if (ynames[1] != "none") {
    pl = ggplot() + 
      geom_tile(data = m.mat, aes(x=Var1, y=as.numeric(Var2), fill=Value)) + ylab('') + xlab('') +
      scale_fill_gradient2(low = "red", high = "blue", mid="white",
                           midpoint = 0) +
      geom_rect(aes(ymin=0.5,ymax=p+0.5,xmin=xdate[1],xmax=xdate[n]),col="black",fill=NA,linetype='dashed') +
      theme(panel.grid = element_blank(), panel.background = element_rect(fill='white'),
            plot.background = element_rect(color=NA), axis.title.x=element_blank(), axis.title.y=element_blank(),
            text = element_text(size=20), legend.position = "right") +
      scale_y_continuous(1:p, breaks=1:p, labels=ynames) 
    }
    
    if (ynames[1] == "none") {
      pl = ggplot() + 
        geom_tile(data = m.mat, aes(x=Var1, y=as.numeric(Var2), fill=Value)) + ylab('') + xlab('') +
        scale_fill_gradient2(low = "red", high = "blue", mid="white",
                             midpoint = 0) +
        geom_rect(aes(ymin=0.5,ymax=p+0.5,xmin=xdate[1],xmax=xdate[n]),col="black",fill=NA,linetype='dashed') +
        theme(panel.grid = element_blank(), panel.background = element_rect(fill='white'),
              plot.background = element_rect(color=NA), axis.title.x=element_blank(), axis.title.y=element_blank(),
              axis.ticks.y=element_blank(),axis.text.y=element_blank(),
              text = element_text(size=20), legend.position = "right") +
        scale_y_continuous(1:p, breaks=1:p, labels=1:p)
    }
  
  plot(pl)
}

getGroups = function(X,mincor=0.1,maxSize=20) {
  clust = corclust(X)
  cl = cvtree(clust,mincor=mincor)$cluster
  group_index = sort(cl)
  cardGroups = tapply(group_index, group_index, length)
  
  grOrd = sort.int(cardGroups,decreasing=T,index.return=T)$ix
  
  newGr = NULL
  varOrd = NULL
  for (v in 1:length(grOrd)) {
    newGr = c(newGr,rep(v,length(which(cl==grOrd[v]))))
    varOrd = c(varOrd,which(cl==grOrd[v]))
  }
  
  redGr = 0
  for (v in 1:length(unique(newGr))) {
    vid = which(newGr==v)
    lvid = length(vid)
    
    if (lvid > maxSize) {
      if (floor(lvid/maxSize) == (lvid/maxSize)) {
        idv = rep(1:floor(lvid/maxSize),each=maxSize)
      }
      if (floor(lvid/maxSize) < (lvid/maxSize)) {
        idv = c(rep(1:floor(lvid/maxSize),each=maxSize),rep(floor(lvid/maxSize)+1,lvid-length(rep(1:floor(lvid/maxSize),each=maxSize))))
      }
      idv = idv + max(redGr)
    }
    
    if (lvid <= maxSize) idv = rep(1,lvid) + max(redGr)
    
    redGr = c(redGr,idv)
  }
  redGr = redGr[-1]
  
  list(group_index=redGr,
       varOrder=varOrd)
}

getQuantities = function(M,nb,bl) {
  B = -M[1:bl,(bl+1):(2*bl)]
  dBl = matrix(fatdiag(M,nb,bl),bl*bl,nb)
  A = array(NA,dim=c(bl,bl,nb))
  for (j in 1:nb) A[,,j] = matrix(dBl[,j],bl,bl)
  
  out = getQuantitiesCpp(A, B, nb, bl)
  D = out$D
  oD = out$oD
  
  list(oD = t(oD[,-1]), D=D)
}

