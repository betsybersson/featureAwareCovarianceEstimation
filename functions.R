library(truncnorm) # for truncated normal for sampling lod
library(LaplacesDemon) # for LaplacesDemon::dmvt (multivariate Student's t) for cusp
##########################################################
### Feature Aware Covariance Estimation Model with CUSP prior on Gamma Gibbs Sampler- FACE_Gcusp_GS
## data
#
#
#
## FACE_GS code copied and edited to add cusp prior on Gamma and 
## remove Gamma dependence on T
## CUSP code adapted from https://github.com/siriolegramanti/CUSP/blob/master/cusp.R
#
#
#
# Y: n x p response matrix
# X: p x q meta covariate matrix (optional)
# k: dimension of latent factor
## classic GS parameters
# S: number of iterations for GS
# burnin: number of iterations to burnin
# thin: thinning mechanism
# save.all: logical to determine what parameters to return. If 0, save subset 
#       of params post thin/burnin. If 1, save and return all iterations 
#       of all parameters. If ICO, save and return only the sum of the inverted 
#       sampling covariance estimate (useful if large p,n)
## optional parameters
# 
##########################################################
FACE_cusp_GS = function(Y,X = NA,
                       k = round(ncol(Y)/2),
                       S = 5100,
                       burnin = 100,
                       thin = 1,
                       save.all = 0,
                       my.seed = Sys.time(),
                       #### CUSP PARAMETERS
                       alpha = round(k/2), # prior expected number of active factors
                       a.theta = 2, b.theta = 2, # prior IG parameters on nu
                       theta.inf = 0.05,
                       #### LOD PREDICTION
                       which.lod = NA,
                       row.ind.lod = NA,
                       lod = NA,
                       #### GEN. RIDGE VARIABLE SELECTION
                       which.cov.group = NA
){
  
  set.seed(my.seed)
  
  ###########################
  ## set hyper params
  p = ncol(Y)
  n = nrow(Y)
  
  ## handle some things for meta covariates
  ## incl some warnings to make sure dimensions are right
  if (is.vector(X)){
    X = matrix(X,ncol=1)
  }
  
  if (all(is.na(X))){
    # if X is NA, make it the intercept
    X = matrix(1,ncol=1,nrow = p)
  }
  
  if (nrow(X)!=p){
    stop("Y,X dimension mismatch!")
  }
  
  q = ncol(X)
  
  ## IG shape/rate parameters
  a.d = a.xi = a.tau = 2
  b.d = b.xi = b.tau = 2/5
  
  ## CUSP parameters
  # a.theta = b.theta = 2
  # alpha = 1 ####????
  # theta.inf = 1/100
  ###########################
  
  ###########################
  ## initialize values
  # sampling model
  Lambda = matrix(0,ncol = k, nrow = p)
  eta = matrix(0,ncol = k, nrow = n)
  ds = rep(1,p); D = D.inv = diag(ds)
  # factor model
  Gamma = matrix(0,nrow = q, ncol = k)
  tau2 = 1; Tau = Tau.inv = eye(k)*tau2
  # CUSP
  theta.inv = theta = rep(1,k)
  Theta = diag(1/theta); Theta.inv = diag(theta.inv)
  z = rep(2,k)
  pi = rep(1,k)/k
  nu = rep(.2,k); nu[k] = 1
  omega.function = function(nu){
    omega = array(NA,dim=k)
    omega[1] = nu[1]
    for (h in 2:k){
      omega[h]<-nu[h]*prod(1-nu[1:(h-1)]) 
    }
    return(omega)
  }
  omega = omega.function(nu)
  if (!any(is.na(which.lod))){ ### bad- writing for my fake example where each column has same # of missingness
    n.lod = length(row.ind.lod)
    y.lod = matrix(0,nrow = n.lod,ncol=p)
  } else{
    n.lod = 0
    y.lod = NA
  }
  ###########################
  
  ###########################
  ## organize shit for generalized lasso
  L = L.inv = eye(q)
  if ( !any(is.na(which.cov.group)) ){
    q.tilde = table(which.cov.group)
    n.q.tilde = length(q.tilde)
    ls = rep(1,n.q.tilde)
  } else {
    n.q.tilde = 1
    ls = 1
  }
  ###########################
  
  
  ###########################
  ## create storage 
  if (save.all == 1){
    n.save.out = S
  } else if (save.all == 0 | save.all == "ICO") {
    len.out = length(seq(from = burnin+1, to = S, by = thin))
    n.save.out = len.out
    index = 1  
  } 
  if (save.all == "ICO"){
    cov.inv.out = cov.inv.det.out = 0
  } else {
    cov.out = cov.inv.out = cov.shrinkto.inv.out = array(NA,dim = c(n.save.out,p*p))
    Lambda.out = array(NA,dim = c(n.save.out,p*k))
    eta.out = array(NA,dim=c(n.save.out, k*n))
    ds.out = array(NA,dim=c(n.save.out,p))
    Gamma.out = array(NA,dim=c(n.save.out,(q)*k))
    vars.out = array(NA,dim=c(n.save.out,1)) #for tau2
    y.lod.out = array(NA,dim=c(n.save.out,n.lod*p))
    ls.out = array(NA,dim=c(n.save.out,n.q.tilde))
  }
  
  ###########################
  
  ###########################
  ## global helpers
  XXt = X %*% t(X)
  ###########################
  
  ###########################
  ## GS
  tic = Sys.time()
  for ( s in 1:S ){
    
    ## sample y.lod
    # sample values below LOD from truncated normal full conditional
    if ( n.lod>0){
      for (i in row.ind.lod){
        
        ind_zz = which(which.lod[i,] == 1)
        
        mu_zz = c(Lambda[ind_zz,] %*% eta[i,])
        d_sqrt_zz = ds[ind_zz] |> sqrt()
        lod_zz = lod[ind_zz]
        
        temp = sapply(1:length(ind_zz), function(j)
          rtruncnorm(1,b = lod_zz[j],
                     mean = mu_zz[j],
                     sd = d_sqrt_zz[j]))
        
        Y[i,ind_zz] = y.lod[which(row.ind.lod==i),ind_zz] = temp
        
      }
    }
    
    ## sample etais
    R = t(Lambda) %*% D.inv
    S.eta.inv = qr.solve(R %*% Lambda + eye(k))
    M.eta.help = S.eta.inv %*% R
    
    eta.list = mclapply(1:n,function(i) rmvnorm(M.eta.help %*% Y[i,],S.eta.inv))
    eta = rbind.list(eta.list)
    
    ## sample djs
    for ( j in 1:p ){
      lambda.j = Lambda[j,]
      helper = (lambda.j - t(Gamma) %*% X[j,])
      S.dj = sum((Y[,j] -  eta %*% lambda.j)^2) + t(helper) %*% helper / tau2 + b.d
      ds[j] = 1/rgamma(1,(n+k+a.d)/2,c(S.dj)/2)
    }
    # update diagonal matrices
    D = diag(ds); D.inv = diag(1/ds)
    
    ## sample Lambda
    S.gamma.inv = qr.solve(Theta.inv / tau2 + t(eta) %*% eta)
    M.gamma = X %*% Gamma %*% Theta.inv / tau2 + t(Y) %*% eta
    Lambda = rmatnorm(M.gamma %*% S.gamma.inv,
                      V = S.gamma.inv, U = D)
    
    ## sample Gamma
    R = t(X) %*% D.inv / tau2
    S.gamma.inv = qr.solve(R %*% X + L.inv)
    M.gamma = R %*% Lambda
    Gamma = rmatnorm(S.gamma.inv %*% M.gamma, V = Theta, U = S.gamma.inv)
    
    
    ## sample tau2
    helper = Lambda - X %*% Gamma
    S.tau = mat_trace(t(helper) %*% D.inv %*% helper %*% Theta.inv) + b.tau
    tau2 = 1/rgamma(1, (p*k + a.tau)/2, S.tau/2)
    Tau = eye(k)*tau2; Tau.inv = eye(k)/tau2
    
    
    ### CUSP- code adapted from https://github.com/siriolegramanti/CUSP/blob/master/cusp.R
    ## sample z
    XLXt = X %*% L %*% t(X)
    lhd_spike<-rep(0,k)
    lhd_slab<-rep(0,k)
    for (h in 1:k){
      lhd_spike[h]<-exp(sum(log(dnorm(Lambda[,h], mean = 0, sd = theta.inf^(1/2), log = FALSE))))
      ## numerical issue - force to be symmetric
      tempS = (b.theta/a.theta)*(XLXt + tau2 * D)
      tempS[lower.tri(tempS)] <- t(tempS)[lower.tri(tempS)]
      lhd_slab[h]<-LaplacesDemon::dmvt(x=Lambda[,h], mu=rep(0,p), S=tempS, 
                                       df=2*a.theta)
      pi<-omega*c(rep(lhd_spike[h],h),rep(lhd_slab[h],k-h))
      if (sum(pi)==0){
        pi<-c(rep(0,k-1),1)
      } else {
        pi<-pi/sum(pi)
      }
      z[h] <- c(1:k)%*%rmultinom(n=1, size=1, prob=pi)
    }
    
    ## sample nu
    for (h in 1:(k-1)){
      nu[h] <- rbeta(1, shape1 = 1+sum(z == h), 
                     shape2 = alpha+sum(z > h))
    }
    
    ## compute omega
    omega = omega.function(nu)
    
    ## sample theta
    for (h in 1:k){
      if (z[h]<=h){
        theta.inv[h]<-theta.inf^(-1)
      } else {
        theta.inv[h] <- 
          rgamma(n = 1,
                 shape = a.theta + 0.5*p,
                 rate = b.theta + 0.5*t(Lambda[,h]) %*% solve(XLXt + tau2 * D) %*% Lambda[,h])
      }
    }
    Theta = diag(1/theta.inv); Theta.inv = diag(theta.inv)
    
    
    ## sample ls
    if (!any(is.na(which.cov.group))){
      for (j in 1:n.q.tilde){
        Gamma.tilde = Gamma[which.cov.group == j,,drop=F]
        ls[j] = 1/rgamma(1,
                         (1+k*q.tilde[j])/2,
                         (1 + mat_trace(Gamma.tilde %*% Theta.inv %*% t(Gamma.tilde))/2)
        )
        if (sum(which.cov.group == j)>1){
          diag(L[which.cov.group == j,which.cov.group == j]) = ls[j]
          diag(L.inv[which.cov.group == j,which.cov.group == j]) = 1/ls[j]
        } else if (sum(which.cov.group == j)==1){
          L[which.cov.group == j,which.cov.group == j]  = ls[j]
          L.inv[which.cov.group == j,which.cov.group == j]  = 1/ls[j]
        }
        
      }
    }
    
    
    
    ## save output
    if (save.all == 0){
      if ((s>burnin)&((s %% thin)==0)){
        
        cov.temp = Lambda %*% t(Lambda) + D
        cov.inv.temp = qr.solve(cov.temp)
        
        cov.shrinkto = D * k * tau2 + X %*% Gamma %*% t(Gamma) %*% t(X)
        
        cov.out[index,] = c(cov.temp)
        cov.inv.out[index,] = c(cov.inv.temp)
        Lambda.out[index,] = c(Lambda)
        eta.out[index,] = c(eta)
        ds.out[index,] = c(ds)
        Gamma.out[index,] = c(Gamma)
        vars.out[index,] = c(tau2) #for tau2 
        cov.shrinkto.inv.out[index,] = c(qr.solve(cov.shrinkto))
        y.lod.out[index,] = c(y.lod)
        ls.out[index,] = c(ls)
        
        index = index + 1
      }
    } else if (save.all == 1){
      
      cov.temp = Lambda %*% t(Lambda) + D
      cov.inv.temp = qr.solve(cov.temp)
      
      cov.shrinkto = D * k * tau2 + X %*% Gamma %*% t(Gamma) %*% t(X)
      
      cov.out[s,] = c(cov.temp)
      cov.inv.out[s,] = c(cov.inv.temp)
      Lambda.out[s,] = c(Lambda)
      eta.out[s,] = c(eta)
      ds.out[s,] = c(ds)
      Gamma.out[s,] = c(Gamma)
      vars.out[s,] = c(tau2) #for tau2
      cov.shrinkto.inv.out[s,] = c(qr.solve(cov.shrinkto))
      y.lod.out[index,] = c(y.lod)
      ls.out[index,] = c(ls)
      
      
    } else if (save.all == "ICO"){
      cov.temp = Lambda %*% t(Lambda) + D
      
      cov.inv.out = cov.inv.out + c(qr.solve(cov.temp))
      cov.inv.det.out = cov.inv.det.out + det(cov.temp)
    }
    
  } # end GS iteration
  ###########################
  runtime = difftime(Sys.time(),tic,units = "secs")
  
  ###########################
  ## put output in list and save to function
  if (save.all == 0 | save.all == 1){
    colnames(vars.out) = c("tau2")
    
    func.out = list("cov" = cov.out,"cov.inv" = cov.inv.out,
                    "cov.shrinkto.inv" = cov.shrinkto.inv.out,
                    "Lambda" = Lambda.out,
                    "eta" = eta.out,
                    "D" = ds.out,
                    "Gamma" = Gamma.out,
                    "vars" = vars.out,
                    "runtime" = runtime,
                    "y.lod" = y.lod.out,
                    "ls" = ls.out)
    if (save.all == 1){
      func.out = func.out #c(func.out,list("out2" = out2))
    }
  } else if (save.all == "ICO"){
    func.out = list("cov.inv.mean" = cov.inv.out/n.save.out,
                    "cov.inv.det.mean" = cov.inv.det.out/n.save.out)
  }
  
  set.seed(Sys.time())
  
  return(func.out)
  ###########################
  
} # end function

##########################################################################
# Adaptive Gibbs sampler for Gaussian factor models under the CUSP prior
# source:https://github.com/siriolegramanti/CUSP/blob/master/cusp.R
##########################################################################
cusp_factor_adapt <- function(y,my_seed,N_sampl,alpha,a_sig,b_sig,a_theta,b_theta,
                              theta_inf,start_adapt,Hmax,alpha0,alpha1,
                              save.all = 0,
                              burnin = round(N_sampl/10),
                              thin = 1){
  set.seed(my_seed)
  # statistical units
  n<-dim(y)[1]
  # number of variables
  p<-dim(y)[2]
  # uniform rvs for adaptation
  u<-runif(N_sampl)
  ########################################################################
  # VARIABLES TO UPDATE
  ########################################################################
  # number of factors
  H<-rep(NA,N_sampl)
  # number of active factors
  Hstar<-rep(NA,N_sampl)
  # factor loadings matrix
  Lambda_post<-vector("list",N_sampl)
  # factor scores matrix
  eta_post<-vector("list",N_sampl)
  # column-specific loadings precisions
  theta_inv_post<-vector("list",N_sampl)
  # stick-breaking weights
  w_post<-vector("list",N_sampl)
  # augmented data
  z_post<-vector("list",N_sampl)
  # error variances
  inv_sigma_sq_post<-matrix(NA,p,N_sampl)
  ########################################################################
  #INITIALIZATION
  ########################################################################
  H[1]<-p+1
  Hstar[1]<-p
  Lambda_post[[1]]<-matrix(rnorm(p*H[1]),p,H[1])
  eta_post[[1]]<-matrix(rnorm(n*H[1]),n,H[1])
  theta_inv_post[[1]]<-rep(1,H[1])
  w_post[[1]]<-rep(1/H[1],H[1])
  z_post[[1]]<-rep(H[1],H[1])
  inv_sigma_sq_post[,1]<-1
  ########################################################################
  # GIBBS SAMPLING #######################################################
  ########################################################################
  t0 <- Sys.time()
  for (t in 2:N_sampl){
    Lambda_post[[t]]<-matrix(NA,p,H[t-1])
    eta_post[[t]]<-matrix(NA,n,H[t-1])
    theta_inv_post[[t]]<-rep(NA,H[t-1])
    w_post[[t]]<-rep(NA,H[t-1])
    z_post[[t]]<-rep(NA,H[t-1])
    # 1) sample the j-th row of Lambda
    for (j in 1:p){
      V_j<-chol2inv(chol(diag(theta_inv_post[[t-1]],nrow=H[t-1])+inv_sigma_sq_post[j,t-1]*(t(eta_post[[t-1]])%*%eta_post[[t-1]])))
      mu_j<-V_j%*%t(eta_post[[t-1]])%*%(y[,j])*inv_sigma_sq_post[j,t-1]
      Lambda_post[[t]][j,]<-rmvnorm(mu_j,V_j) #mvrnorm(1,mu_j,V_j)
    }
    # 2) sample sigma^{-2}
    diff_y<-apply((y-(eta_post[[t-1]]%*%t(Lambda_post[[t]])))^2,2,sum)
    for (j in 1:p){
      inv_sigma_sq_post[j,t]<-rgamma(n=1,shape=a_sig+0.5*n,rate=b_sig+0.5*diff_y[j])
    }		
    # 3) sample eta
    V_eta<-chol2inv(chol(diag(H[t-1])+t(Lambda_post[[t]])%*%diag(inv_sigma_sq_post[,t])%*%Lambda_post[[t]]))
    for (i in 1:n){
      mu_eta_i<-V_eta%*%t(Lambda_post[[t]])%*%diag(inv_sigma_sq_post[,t])%*%y[i,]
      eta_post[[t]][i,]<-rmvnorm(mu_eta_i,V_eta) #mvrnorm(1,mu_eta_i,V_eta)
    }	
    # 4) sample z
    lhd_spike<-rep(0,H[t-1])
    lhd_slab<-rep(0,H[t-1])
    for (h in 1:H[t-1]){
      lhd_spike[h]<-exp(sum(log(dnorm(Lambda_post[[t]][,h], mean = 0, sd = theta_inf^(1/2), log = FALSE))))
      lhd_slab[h]<-LaplacesDemon::dmvt(x=Lambda_post[[t]][,h], mu=rep(0,p), S=(b_theta/a_theta)*diag(p), df=2*a_theta)
      prob_h<-w_post[[t-1]]*c(rep(lhd_spike[h],h),rep(lhd_slab[h],H[t-1]-h))
      if (sum(prob_h)==0){
        prob_h<-c(rep(0,H[t-1]-1),1)
      }
      else{
        prob_h<-prob_h/sum(prob_h)
      }
      z_post[[t]][h]<-c(1:H[t-1])%*%rmultinom(n=1, size=1, prob=prob_h)
    }
    # 5) sample v and update w
    v<-rep(NA,H[t-1])
    for (h in 1:(H[t-1]-1)){
      v[h]<-rbeta(1, shape1 = 1+sum(z_post[[t]]==h), shape2 = alpha+sum(z_post[[t]]>h))
    }
    v[H[t-1]]<-1
    w_post[[t]][1]<-v[1]
    for (h in 2:H[t-1]){
      w_post[[t]][h]<-v[h]*prod(1-v[1:(h-1)])  
    }
    # 6) sample theta^{-1}
    for (h in 1:H[t-1]){
      if (z_post[[t]][h]<=h){
        theta_inv_post[[t]][h]<-theta_inf^(-1)
      }
      else{
        theta_inv_post[[t]][h]<-rgamma(n=1,shape=a_theta+0.5*p,rate=b_theta+0.5*t(Lambda_post[[t]][,h])%*%Lambda_post[[t]][,h])
      }
    }
    # 7) update H[t]
    active <- which(z_post[[t]]>c(1:H[t-1]))
    Hstar[t] <- length(active)
    H[t]<-H[t-1]
    # by default, keep the same truncation, unless...
    if (t>=start_adapt & u[t]<=exp(alpha0+alpha1*t)){
      if (Hstar[t]<H[t-1]-1){
        # set truncation to Hstar[t] and subset all variables, keeping only active columns
        H[t]<-Hstar[t]+1
        eta_post[[t]]<-cbind(eta_post[[t]][,active],rnorm(n))
        theta_inv_post[[t]]<-c(theta_inv_post[[t]][active],theta_inf^(-1))
        w_post[[t]]<-c(w_post[[t]][active],1-sum(w_post[[t]][active]))
        Lambda_post[[t]]<-cbind(Lambda_post[[t]][,active],rnorm(p,mean=0,sd=sqrt(theta_inf)))
      } else if (H[t-1]<Hmax) {
        # increase truncation by 1 and extend all variables, sampling from the prior/model
        H[t]<-H[t-1]+1
        eta_post[[t]]<-cbind(eta_post[[t]],rnorm(n))
        v[H[t-1]]<-rbeta(1,shape1=1,shape2=alpha)
        v<-c(v,1)
        w_post[[t]]<-rep(NA,H[t])
        w_post[[t]][1]<-v[1]
        for (h in 2:H[t]){
          w_post[[t]][h]<-v[h]*prod(1-v[1:(h-1)])  
        }
        theta_inv_post[[t]]<-c(theta_inv_post[[t]],theta_inf^(-1))
        Lambda_post[[t]]<-cbind(Lambda_post[[t]],rnorm(p,mean=0,sd=sqrt(theta_inf)))
      }
    }
    ##print(t)
  }
  runtime <- difftime(Sys.time(),t0,units = "secs")
  
  ## get thinned cov.inv - Betsy added
  thinning <- seq(burnin+1,N_sampl,by=thin)
  
  cov.inv <- array(NA,c(length(thinning),p*p))
  
  for (i in 1:length(thinning)){
    cov.inv[i,] <- c(qr.solve(Lambda_post[[thinning[i]]]%*%t(Lambda_post[[thinning[i]]]) +
                                diag(1/inv_sigma_sq_post[,thinning[i]]))
    )
  }
  
  ## save output
  output = list("runtime"=runtime,"cov.inv" = cov.inv)
  if (save.all == 1){
    output = c(output,c("y"=y,"my_seed"=my_seed,"N_sampl"=N_sampl,
                        "alpha"=alpha,"a_sig"=a_sig,"b_sig"=b_sig,"a_theta"=a_theta,"b_theta"=b_theta,
                        "theta_inf"=theta_inf,"start_adapt"=start_adapt,"Hmax"=Hmax,
                        "alpha0"=alpha0,"alpha1"=alpha1,
                        "H"=H,"Hstar"=Hstar,"Lambda_post"=Lambda_post,"eta_post"=eta_post,
                        "theta_inv_post"=theta_inv_post,"w_post"=w_post,"z_post"=z_post,
                        "inv_sigma_sq_post"=inv_sigma_sq_post)
    )
  }
  
  return(output)
}
##########################################################
### FACE Function: get design matrix for matrix-variate data structure
# for matrix variate data Y_i (p1 x p2) <=> multivariate data y_i (p1*p2 x 1)
# p1: row dimension of matrix-variate data
# p2: column dimension of matrix-variate data
# output: X (p1*p2) x (p1+p2) design matrix indicating row/column group membership of each of the (p1*p2) variables
##########################################################
getMatDesignMat = function(p1,p2){
  
  p = p1*p2
  
  X.R = matrix(0,nrow = p, ncol = p1)
  X.C = matrix(0,nrow = p, ncol = p2)
  colnames(X.R) = paste0("R",1:p1)
  colnames(X.C) = paste0("C",1:p2)
  Y.temp = matrix(0,nrow = p1, ncol = p2)
  for ( j in 1:p1 ){
    Y.temp[j,] = 1
    X.R[,j] = c(Y.temp)
    Y.temp = Y.temp * 0
  }
  for ( j in 1:p2 ){
    Y.temp[,j] = 1
    X.C[,j] = c(Y.temp)
    Y.temp = Y.temp * 0
  }
  
  X = cbind(X.R,X.C)
  
  return(X)
}
##########################################################
### FACE Function: get design matrix for categorical vector
# CV: vector of lenght p consisting of category membership for G<=p unique categories
# output: X (p) x (G) design matrix indicating row/column group membership of each of the (p1*p2) variables
##########################################################
getCatDesignMat = function(CV){
  
  # get design matrix
  p = length(CV)
  cvs = unique(CV)
  g = length(cvs)
  X = matrix(0, nrow = p, ncol = g)
  for ( g.ind in 1:g){
    X[CV==cvs[g.ind],g.ind] = 1
  }
  colnames(X) = cvs
  return(X)
}

##########################################################
### Separable covariance MLE
## block coordinate descent algorithm
# X: p1 x p2 x n array, de-meaned
# output: Psi p1 x p1 row covariance
# output: Sig p2 x p2 column covariance
# de.meaned: logical. If de.meaned = true, use n-1 instead of n in degrees of freedom on S
##########################################################
cov.kron.mle = function(X,itmax = 100,eps = 1e-5,
                        de.meaned = F, my.seed = Sys.time()){
  
  set.seed(my.seed)
  
  p1 = dim(X)[1]
  p2 = dim(X)[2]
  N = dim(X)[3] # since removing mean
  
  if (de.meaned == T){
    N = N-1
  }
  
  # initialize sig tilde
  Sig.tilde = matrix(rowMeans(apply(X,3,cov.mle)),ncol = p2)
  Sig.tilde.inv = solve(Sig.tilde)
  
  # initialize stopping checks
  Psi.tilde = Psi.tilde.old = matrix(0,ncol = p1,nrow = p1)
  Sig.tilde.old = matrix(0,ncol = p2, nrow = p2)
  check = F
  it = 0
  tic = Sys.time()
  while((check == FALSE) & (it < itmax)){
    
    Psi.tilde = matrix(rowSums(apply(X,3,function(KK)KK %*% Sig.tilde.inv %*% t(KK))),
                       ncol = p1)/(N*p2)
    Psi.tilde.inv = solve(Psi.tilde)
    
    Sig.tilde = matrix(rowSums(apply(X,3,function(KK)t(KK) %*% Psi.tilde.inv %*% KK)),
                       ncol = p2)/(N*p1)
    Sig.tilde.inv = solve(Sig.tilde)
    
    if (all(abs(Sig.tilde.old-Sig.tilde)<eps) &
        all(abs(Psi.tilde.old-Psi.tilde)<eps)){
      check = TRUE
    }
    
    # update for next iteration in while loop
    it = it+1 
    Psi.tilde.old = Psi.tilde
    Sig.tilde.old = Sig.tilde
    
  }
  
  runtime = difftime(Sys.time(),tic,units = "secs")
  
  set.seed(Sys.time())
  
  return(list("Psi" = Psi.tilde,"Sigma" = Sig.tilde,
              "Cov" = kronecker(Sig.tilde,Psi.tilde),
              "it" = it, "runtime" = runtime))
  
}
