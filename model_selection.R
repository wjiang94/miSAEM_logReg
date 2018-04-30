library(MASS)
library(mvtnorm)
source("saem.R")

N <- 1000  # number of subjects
p <- 5     # number of explanatory variables
mu.star <- 1:p  #rep(0,p)  # mean of the explanatory variables
sd <- 1:p # rep(1,p) # standard deviations



# ------- generating missing data

p.miss <- 0.10 


combinations = function(n){
  comb = NULL
  if (n<15) {
    for( i in 1:n) comb = rbind(cbind(1,comb),cbind(0,comb))
    return(comb)
  }
  else {error("this value will probably block your computer, try on your own risk")}
}

subsets=combinations(p)

nb.simu = 100

subsets=combinations(p)
#subsets= subsets[rowSums(subsets)==4,]

ll=ll1 = AIC =AIC1 = BIC= BIC1 = matrix(0,nrow= nb.simu,ncol = nrow(subsets)-1)

AIC_min =BIC_min=AIC1_min =BIC1_min= matrix(1e+5,nrow= nb.simu,ncol = p)
j_AIC = j_BIC =j_ll=j_AIC1 = j_BIC1 =j_ll1 = matrix(0,nrow= nb.simu,ncol = p)
ll_max = ll1_max=matrix(-1e+5,nrow= nb.simu,ncol = p)

AIC_all_min =BIC_all_min=AIC1_all_min =BIC1_all_min= rep(1e+5,nb.simu)
ll_all_max = ll1_all_max=rep(-1e+5,nb.simu)
j_all_AIC = j_all_BIC =j_all_ll=j_all_AIC1 = j_all_BIC1 =j_all_ll1 = rep(0,nb.simu)

# C <- matrix(c(   # correlation matrix
#   1,   0.8, 0,   0,   0,
#   0.8, 1,   0,   0,   0,
#   0,   0,   1,   0.3, 0.6,
#   0,   0,   0.3, 1,   0.7,
#   0,   0,   0.6, 0.7, 1
#  ), nrow=p)
C= diag(5)

Sigma.star <- diag(sd)%*%C%*%diag(sd) # variance-covariance matrix of the explanatory variables

#beta.star <- c(0.5, -0.3, 1, 0, -0.6) # coefficients
#beta0.star <- -0.2 
#beta.star <- c(0.5, 0, 1, 0, -0.6) # coefficients
#beta0.star <- -0.2 
#
beta.star <- c(0.5, 0, 1, 0, -0.6) # coefficients
beta0.star <- -0.2  # intercept

beta.true = c(beta0.star,beta.star)


for(nb in 1:nb.simu){
  set.seed(nb)
  cat('simu ',nb,'\n')
  # ----- complete data simulation
  patterns = runif(N*p)<p.miss
  X.complete <- matrix(rnorm(N*p), nrow=N)%*%chol(Sigma.star) + matrix(rep(mu.star,N), nrow=N, byrow = TRUE)
  p1 <- 1/(1+exp(-X.complete%*%beta.star-beta0.star))
  y <- as.numeric(runif(N)<p1)
  
  
  
  X.obs <- X.complete
  X.obs[patterns] <- NA
  for (j in 1:(nrow(subsets)-1)){
    nb.var = sum(subsets[j,])
    cat('subset ',subsets[j,],'\n')  
    variables = subsets[j,]
    pos_var=which(variables==1)
    nb.x = sum(variables)
    nb.para = (nb.x + 1) + p + p*p 
    list.saem.subset=miss.saem(X.obs,pos_var,y,maxruns=1000,tol_em=1e-7,nmcmc=2,ag=1,k1=5,print_iter=FALSE,ll_obs_cal=TRUE)
    ll[nb,j] = list.saem.subset$ll
    ll1[nb,j] = list.saem.subset$ll1
    AIC[nb,j] = -2*ll[nb,j]+ 2*nb.para
    AIC1[nb,j] = -2*ll1[nb,j]+ 2*(nb.x+1)
    BIC[nb,j] = -2*ll[nb,j]+ nb.para * log(N)
    BIC1[nb,j] = -2*ll1[nb,j]+ (nb.x+1) * log(N)
    
    
    
    if(AIC[nb,j]<=AIC_min[nb,nb.x]){
      AIC_min[nb,nb.x]= AIC[nb,j]
      j_AIC[nb,nb.x] = j
    }
    if(BIC[nb,j]<=BIC_min[nb,nb.x]){
      BIC_min[nb,nb.x]= BIC[nb,j]
      j_BIC[nb,nb.x] = j
    }
    if(ll[nb,j]>=ll_max[nb,nb.x]){
      ll_max[nb,nb.x]= ll[nb,j]
      j_ll[nb,nb.x] = j
    }
    if(AIC1[nb,j]<=AIC1_min[nb,nb.x]){
      AIC1_min[nb,nb.x]= AIC1[nb,j]
      j_AIC1[nb,nb.x] = j
    }
    if(BIC1[nb,j]<=BIC1_min[nb,nb.x]){
      BIC1_min[nb,nb.x]= BIC1[nb,j]
      j_BIC1[nb,nb.x] = j
    }
    if(ll1[nb,j]>=ll1_max[nb,nb.x]){
      ll1_max[nb,nb.x]= ll1[nb,j]
      j_ll1[nb,nb.x] = j
    }
    
    
    if(AIC[nb,j]<=AIC_all_min[nb]){
      AIC_all_min[nb]= AIC[nb,j]
      j_all_AIC[nb] = j
    }
    if(BIC[nb,j]<=BIC_all_min[nb]){
      BIC_all_min[nb]= BIC[nb,j]
      j_all_BIC[nb] = j
    }
    if(ll[nb,j]>=ll_all_max[nb]){
      ll_all_max[nb]= ll[nb,j]
      j_all_ll[nb] = j
    }
    if(AIC1[nb,j]<=AIC1_all_min[nb]){
      AIC1_all_min[nb]= AIC1[nb,j]
      j_all_AIC1[nb] = j
    }
    if(BIC1[nb,j]<=BIC1_all_min[nb]){
      BIC1_all_min[nb]= BIC1[nb,j]
      j_all_BIC1[nb] = j
    }
    if(ll1[nb,j]>=ll1_all_max[nb]){
      ll1_all_max[nb]= ll1[nb,j]
      j_all_ll1[nb] = j
    }
  }
}



par(mfrow=c(1,1))
plot(ll[1,],xlab="model",ylab="observed log-likelihood")
#for(nb in 2:10 ){lines(ll[nb,])}
abline(v = 11, col = "red", lty = 2)

par(mfrow=c(1,1))
plot(AIC[1,],xlab="model",ylab="AIC")
#for(nb in 2:10 ){lines(AIC[nb,])}
abline(v = 11, col = "red", lty = 2)

par(mfrow=c(1,1))
plot(BIC[1,],xlab="model",ylab="BIC")
#for(nb in 2:10 ){lines(BIC[nb,])}
abline(v = 11, col = "red", lty = 2)


par(mfrow=c(3,1))
plot(AIC_min[1,])
#for (i in 1:10){lines(AIC_min[i+1,])}
abline(v = 4, col = "red", lty = 2)

plot(BIC_min[1,])
#for (i in 1:10){lines(BIC_min[i+1,])}
abline(v = 4, col = "red", lty = 2)

plot(ll_max[1,])
#for (i in 1:10){lines(ll_max[i+1,])}
abline(v = 4, col = "red", lty = 2)

