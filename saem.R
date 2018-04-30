
# ---------- SAEM
#X.obs:x1,x2
#X.all:x1,x2,x3,x4,x5
miss.saem <- function(X.obs,pos_var,y,maxruns=500,tol_em=1e-7,nmcmc=2,ag=1,k1=50,print_iter=TRUE, var_obs_cal=FALSE, ll_obs_cal=FALSE) {
#delete rows completely missing
  p=ncol(X.obs)#x1,x2
  if(any(apply(is.na(X.obs),1,sum)==p)){
    i_allNA=which(apply(is.na(X.obs),1,sum)==p)
    X.obs = X.obs[-i_allNA,]
    y = y[-i_allNA]
  }
  if(any((is.na(y))==TRUE)){
    i_YNA=which(is.na(y)==TRUE)
    X.obs = X.obs[-i_YNA,]
    y = y[-i_YNA]
  }
  n=length(y)

  
  rindic = as.matrix(is.na(X.obs))
  if(sum(rindic)>0){
    whichcolmissing = (1:ncol(rindic))[apply(rindic,2,sum)>0] #2 4 5
    missingcols = length(whichcolmissing) #3
  }
  if(sum(rindic)==0){missingcols=0}
  
  
  ptm <- Sys.time()
  if(missingcols>0){
    k=0
    cstop=0.1
    seqbeta = matrix(NA,nrow=ncol(X.obs)+1,ncol=(maxruns+1))
    seqbeta_avg = matrix(NA,nrow=ncol(X.obs)+1,ncol=(maxruns+1))
    
    X.mean = X.obs
    for(i in 1:ncol(X.mean)){
      X.mean[is.na(X.mean[,i]), i] <- mean(X.mean[,i], na.rm = TRUE)
    }
    X.sim <- X.mean
    
    mu = apply(X.mean,2,mean)
    Sigma = var(X.mean)*(n-1)/n
    beta= rep(0,p+1)
    beta[c(1,pos_var+1)]= glm(y~ X.mean[,pos_var],family=binomial(link='logit'))$coef
    
    while ((cstop>tol_em)*(k<maxruns)|(k<20)){
      k = k+1
      beta.old = beta
      
      if(k <k1){gamma <- 1}else{gamma <- 1/(k-(k1-1))^ag}
      
      S.inv <- solve(Sigma)
      
      for (i in (1:n)) {
        jna <- which(is.na(X.obs[i,]))
        njna <- length(jna)
        if (njna>0) {
          xi <- X.sim[i,]
          Oi <- solve(S.inv[jna,jna])
          mi <- mu[jna] 
          lobs <- beta[1]
          if (njna<p) {
            jobs <- setdiff(1:p,jna)
            mi <- mi - (xi[jobs] - mu[jobs])%*%S.inv[jobs,jna]%*%Oi
            lobs <- lobs + sum(xi[jobs]*beta[jobs+1])
          }
          
          cobs <- exp(lobs)
          if(cobs==0){cobs=.Machine$double.xmin}
          if(cobs==Inf){cobs=.Machine$double.xmax}
          
          xina <- xi[jna]
          betana <- beta[jna+1]
          for (m in (1:nmcmc)) {
            xina.c <- mi + rnorm(njna)%*%chol(Oi)
            
            if (y[i]==1) 
              alpha <- (1+exp(-sum(xina*betana))/cobs)/(1+exp(-sum(xina.c*betana))/cobs)
            else
              alpha <- (1+exp(sum(xina*betana))*cobs)/(1+exp(sum(xina.c*betana))*cobs)
            if (runif(1) < alpha){
              xina <- xina.c
            }
          }
          X.sim[i,jna] <- xina
        }
      }
      beta_new= rep(0,p+1)
      beta_new[c(1,pos_var+1)]= glm(y~ X.sim[,pos_var],family=binomial(link='logit'))$coef
      # data.sim <- data.frame(y=y,X.sim)
      # model.sim <- glm(y ~.,family=binomial(link='logit'),data=data.sim)
      beta <- (1-gamma)*beta + gamma*beta_new
      cstop = sum((beta-beta.old)^2)
      
      mu <- (1-gamma)*mu + gamma*colMeans(X.sim)
      Sigma <- (1-gamma)*Sigma + gamma*cov(X.sim)
      
      seqbeta[,k]=beta.old
      
      if(k==1){
        seqbeta_avg[,k]=beta.old
        }else{
        seqbeta_avg[,k]= 1/k*rowSums(seqbeta[,1:k])
        }
      
      if(print_iter==TRUE){
        cat(sprintf('iteration = %i ', k))
        cat(sprintf('beta ='),beta,'\n')
        cat(sprintf('Distance from estimateur from last iter ='),cstop,'\n')
      }
    }
    var_obs = ll = ll1= ll2=std_obs =NULL
    if(var_obs_cal==TRUE){
      var_obs = louis_lr_saem(beta,mu,Sigma,y,X.obs,pos_var,rindic,whichcolmissing,mc.size=1000)
      std_obs <- sqrt(diag(var_obs))
    }
    if(ll_obs_cal==TRUE){
      ll_obs = likelihood_saem(beta,mu,Sigma,y,X.obs,rindic,whichcolmissing,mc.size=1000)
      ll = ll_obs$ll
      ll1 = ll_obs$ll1
      ll2 = ll_obs$ll2
    }
  }
  if(missingcols==0){
    X.obs = matrix(X.obs,nrow=n)
    data.complete <- data.frame(y=y,X.obs)
    model.complete <- glm(y ~.,family=binomial(link='logit'),data=data.complete)
    mu = apply(X.obs,2,mean)
    Sigma = var(X.obs)*(n-1)/n
    beta <- model.complete$coefficients
    var_obs = ll = ll1 =ll2= std_obs =seqbeta_avg= seqbeta=NULL
    if(var_obs_cal==TRUE){
      P <- predict(model.complete, type = "response")
      W <- diag(P*(1-P))
      X <- model.matrix(model.complete)
      
      var_obs <- solve(t(X)%*%W%*%X)
      std_obs <- sqrt(diag(var_obs))
    }
    if(ll_obs_cal==TRUE){
      ll_obs = likelihood_saem(beta,mu,Sigma,y,X.obs,rindic,whichcolmissing,mc.size=1000)
      ll = ll_obs$ll
      ll1 = ll_obs$ll1
      ll2 = ll_obs$ll2
    }
  }
  time_run=Sys.time() - ptm
  return(list(mu=mu, sig2=Sigma, beta=beta,time_run=time_run,seqbeta=seqbeta,seqbeta_avg=seqbeta_avg,ll=ll,ll1=ll1,ll2=ll2,var_obs=var_obs,std_obs=std_obs,X.sim=X.sim))
}


louis_lr_saem = function(beta,mu,Sigma,Y,X.obs,pos_var,rindic,whichcolXmissing,mc.size){
  n=dim(X.obs)[1]
  #p=dim(X.obs)[2]
  p = length(pos_var)
  beta = beta[c(1,pos_var+1)]
  mu = mu[pos_var]
  Sigma = Sigma[pos_var,pos_var]
  X.obs = X.obs[,pos_var]
  rindic = rindic[,pos_var]
  
  X.mean = X.obs
  for(i in 1:ncol(X.mean)){X.mean[is.na(X.mean[,i]), i] <- mean(X.mean[,i], na.rm = TRUE)}
  X.sim <- X.mean
  G = D = I_obs = matrix(0,ncol = p+1,nrow = p+1)
  Delta = matrix (0,ncol=1,nrow=p+1)
  S.inv <- solve(Sigma)
  
  for (i in (1:n)) {
    jna <- which(is.na(X.obs[i,]))
    njna <- length(jna)
    if (njna==0) {
      x = matrix(c(1,X.sim[i,]))
      exp_b=exp(beta%*%x)[1]
      d2l = -x%*%t(x)*(exp_b/(1+exp_b)^2)
      I_obs = I_obs - d2l}
    if (njna>0) {
      xi <- X.sim[i,]
      Oi <- solve(S.inv[jna,jna])
      mi <- mu[jna] 
      lobs <- beta[1]
      if (njna<p) {
        jobs <- setdiff(1:p,jna)
        mi <- mi - (xi[jobs] - mu[jobs])%*%S.inv[jobs,jna]%*%Oi
        lobs <- lobs + sum(xi[jobs]*beta[jobs+1])
      }
      cobs <- exp(lobs)
      xina <- xi[jna]
      betana <- beta[jna+1]
      for (m in (1:mc.size)) {
        xina.c <- mi + rnorm(njna)%*%chol(Oi)
        if (Y[i]==1) 
          alpha <- (1+exp(-sum(xina*betana))/cobs)/(1+exp(-sum(xina.c*betana))/cobs)
        else
          alpha <- (1+exp(sum(xina*betana))*cobs)/(1+exp(sum(xina.c*betana))*cobs)
        if (runif(1) < alpha){
          xina <- xina.c
        }
        X.sim[i,jna] <- xina
        x = matrix(c(1,X.sim[i,]))
        exp_b=exp(beta%*%x)[1]
        dl = x*(Y[i]-exp_b/(1+exp_b))
        d2l = -x%*%t(x)*(exp_b/(1+exp_b)^2)
        D = D + 1/m * ( d2l - D) 
        G = G + 1/m * ( dl %*% t(dl) - G)
        Delta = Delta + 1/m * ( dl - Delta)
      }
      I_obs = I_obs - (D+G-Delta%*%t(Delta))
    }
  }
  V_obs = solve(I_obs)
  return(V_obs)
}


likelihood_saem = function(beta,mu,Sigma,Y,X.obs,rindic,whichcolXmissing,mc.size){
  n=dim(X.obs)[1]
  p=dim(X.obs)[2]
  lh=0
  lh1=0
  lh2=0
  for(i in 1:n){
    y=Y[i]
    x=X.obs[i,]
    #if(sum(rindic[i,])==0){lh=lh+log_reg(y=y,x=c(1,x),beta,iflog=FALSE)+dmvnorm(x, mean = mu, sigma = Sigma, log = FALSE)}
    if(sum(rindic[i,])==0){
      lr = log_reg(y=y,x=c(1,x),beta,iflog=TRUE)
      nm = dmvnorm(x, mean = mu, sigma = Sigma, log = TRUE)
      lh=lh+lr+nm
      lh1=lh1+lr
      lh2=lh2+nm
      #if(nm <= -100){cat(sprintf('i='),i,sprintf('logLik(lr)='),lr,sprintf('logLik(nm)='),nm,'\n')}
    }
    else{
      miss_col = which(rindic[i,]==TRUE)
      x2 = x[-miss_col]
      mu1 = mu[miss_col]
      mu2 = mu[-miss_col]
      sigma11 = Sigma[miss_col,miss_col]
      sigma12 = Sigma[miss_col,-miss_col]
      sigma22 = Sigma[-miss_col,-miss_col]
      sigma21 = Sigma[-miss_col,miss_col]
      mu_cond = mu1+sigma12 %*% solve(sigma22)%*%(x2-mu2)
      sigma_cond = sigma11 - sigma12 %*% solve(sigma22) %*% sigma21
      x1_all=mvrnorm(n = mc.size, mu_cond, sigma_cond, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
      lh_mis1=0
      for(j in 1:mc.size){
        x[miss_col] =x1= x1_all[j,]
        lh_mis1 = lh_mis1 + log_reg(y=y,x=c(1,x),beta,iflog=FALSE)
      }
      lr = log(lh_mis1/mc.size)
      lh1 = lh1 + lr
      if(length(x2)>1){nm = dmvnorm(x2, mean = mu2, sigma = sigma22, log = TRUE)}
      if(length(x2)==1){nm = dnorm(x2, mean=mu2, sd=sqrt(sigma22), log = TRUE)}
      lh2 = lh2 + nm
      lh = lh + lr + nm
      #if(nm <= -100){cat(sprintf('i='),i,sprintf('logLik(lr)='),lr,sprintf('logLik(nm)='),nm,'\n')}
    }
  }
  return(list(ll=lh,ll1=lh1,ll2=lh2))
}


log_reg <- function(y,x,beta,iflog=TRUE){
  res <- y*(x%*%beta) - log(1+exp(x%*%beta))
  if(iflog==TRUE)
    return(res)
  else
    return(exp(res))
}
