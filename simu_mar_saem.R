library(MASS)
library(mvtnorm)
library(mice)
source("mcem_louis.R")
source("saem.R")

N <- 1000  # number of subjects
#N <- 10000  # number of subjects
p <- 5     # number of explanatory variables

mu.star <- 1:p  #rep(0,p)  # mean of the explanatory variables
sd <- 1:p # rep(1,p) # standard deviations

# C <- matrix(c(   # correlation matrix
#   1,   0.8, 0,   0,   0,
#   0.8, 1,   0,   0,   0,
#   0,   0,   1,   0.3, 0.6,
#   0,   0,   0.3, 1,   0.7,
#   0,   0,   0.6, 0.7, 1
# ), nrow=p)

C = diag(p)
Sigma.star <- diag(sd)%*%C%*%diag(sd) # variance-covariance matrix of the explanatory variables

beta.star <- c(0.5, -0.3, 1, 0, -0.6) # coefficients
beta0.star <- -0.2  # intercept


nbsim = 100
#e.obs = e.predx = e.predxy = e.saem = f.predx= f.predxy = f.saem = NULL
EST.saem = EST.comp=EST.cc=matrix(0, nbsim,length(beta.star)+1)
BIAS.saem = rep(0, nbsim)
TIME.saem = rep(0, nbsim)
STD.saem = STD.comp= STD.cc= matrix(0, nbsim,length(beta.star)+1)
LENGTH.saem=LENGTH.comp = matrix(0, nbsim,length(beta.star)+1)


#100 simus
beta.true = c(beta0.star,beta.star)
count.saem =count.comp= rep(0,p+1)

for (NB in 1:nbsim){
  #cat(sprintf('simulation = '), NB,'\n')
  
  # ----- complete data simulation
  X.complete <- matrix(rnorm(N*p), nrow=N)%*%chol(Sigma.star) + matrix(rep(mu.star,N), nrow=N, byrow = TRUE)
  p1 <- 1/(1+exp(-X.complete%*%beta.star-beta0.star))
  y <- as.numeric(runif(N)<p1)
  
  data.complete <- data.frame(y=y,X.complete)
  model.complete <- glm(y ~.,family=binomial(link='logit'),data=data.complete)
  beta0.complete <- model.complete$coefficients[1]
  beta.complete <- model.complete$coefficients[2:(p+1)]
  bias.complete  = sum((c(beta0.complete ,beta.complete)-c(beta0.star,beta.star))^2)
  P <- predict(model.complete, type = "response")
  W <- diag(P*(1-P))
  X <- model.matrix(model.complete)
  V_complete <- solve(t(X)%*%W%*%X)
  std.complete <- sqrt(diag(V_complete))
  
  # ------- generating missing data
  X.obs <- X.complete
  for(i in c(2,4,5)){
    z <- cbind(y,X.complete[,c(1,3)])%*%matrix(sample(-5:5, 3, replace=T),ncol=1)        # linear combination 
    pr <- 1/(1+exp(-z))         # pass through an inv-logit function
    r <- rbinom(N,1,pr)      # bernoulli response variable
    X.obs[r==0,i]<-NA
  }
  cat('percentage of NA: ', mean(is.na(X.obs[,2])),mean(is.na(X.obs[,4])),mean(is.na(X.obs[,5])),'\n')
  
  
  # ------- estimation ignoring the missing data
  data.obs <- data.frame(y=y,X.obs)
  model.obs <- glm(y ~.,family=binomial(link='logit'),data=data.obs)
  beta0.cc <- model.obs$coefficients[1]
  beta.cc <- model.obs$coefficients[2:(p+1)]
  bias.cc = sum((c(beta0.cc,beta.cc)-c(beta0.star,beta.star))^2)
  
  
  #list.est = mcem_ar(X= X.obs,Y= y , maxruns=1000,tol_em=1e-5)
  #ptm <- proc.time()
  list.saem=miss.saem(X.obs,y,maxruns=500,tol_em=1e-5,nmcmc=2,ag=1,k1=50, print_iter=FALSE,var_obs_cal=TRUE, ll_obs_cal=FALSE)
  #time_mcem_ar=proc.time() - ptm
  #time_mcem_ar=proc.time() - ptm
  beta.saem = list.saem$beta
  bias.saem = sum((beta.saem -c(beta0.star,beta.star))^2)
  std.saem = list.saem$std_obs
  
  EST.comp[NB,] = c(beta0.complete,beta.complete)
  EST.cc[NB,] = c(beta0.cc,beta.cc)
  EST.saem[NB,] = beta.saem
  
  STD.comp[NB,] = std.complete
  STD.saem[NB,] = std.saem
  
  ci.comp_ceil =  c(beta0.complete,beta.complete) + 1.96*std.complete
  ci.comp_ground =  c(beta0.complete,beta.complete) - 1.96*std.complete
  ci.saem_ceil = beta.saem + 1.96*std.saem
  ci.saem_ground = beta.saem - 1.96*std.saem
  
  LENGTH.comp[NB,] = ci.comp_ceil - ci.comp_ground
  LENGTH.saem[NB,] = ci.saem_ceil - ci.saem_ground
  for(i in 1:(p+1)){
    if( ci.comp_ground[i] <=beta.true[i] & ci.comp_ceil[i]>beta.true[i]){
      count.comp[i]<-count.comp[i]+1
    }
    if( ci.saem_ground[i] <=beta.true[i] & ci.saem_ceil[i]>beta.true[i]){
      count.saem[i]<-count.saem[i]+1
    }
  }
  cat('simulation =',NB,', count(complete) =',count.comp[2],', count(saem) =',count.saem[2], '\n')
}

par(mfrow=c(2,3))
boxplot(EST.comp[,1]-beta0.star,EST.cc[,1]-beta0.star,EST.saem[,1]-beta0.star,main='bias for beta0',names=c("no NA","CC","saem"))
abline(h = 0, col = "red", lty = 2)
boxplot(EST.comp[,2]-beta.star[1],EST.cc[,2]-beta.star[1],EST.saem[,2]-beta.star[1],main='bias for beta1',names=c("no NA","CC","saem"))
abline(h = 0, col = "red", lty = 2)
boxplot(EST.comp[,3]-beta.star[2],EST.cc[,3]-beta.star[2],EST.saem[,3]-beta.star[2],main='bias for beta2',names=c("no NA","CC","saem"))
abline(h = 0, col = "red", lty = 2)
boxplot(EST.comp[,4]-beta.star[3],EST.cc[,4]-beta.star[3],EST.saem[,4]-beta.star[3],main='bias for beta3',names=c("no NA","CC","saem"))
abline(h = 0, col = "red", lty = 2)
boxplot(EST.comp[,5]-beta.star[4],EST.cc[,5]-beta.star[4],EST.saem[,5]-beta.star[4],main='bias for beta4',names=c("no NA","CC","saem"))
abline(h = 0, col = "red", lty = 2)
boxplot(EST.comp[,6]-beta.star[5],EST.cc[,6]-beta.star[5],EST.saem[,6]-beta.star[5],main='bias for beta5',names=c("no NA","CC","saem"))
abline(h = 0, col = "red", lty = 2)


par(mfrow=c(2,3))
boxplot(LENGTH.comp[,1],LENGTH.saem[,1],main='length of CI of beta0',names=c("no NA",'SAEM'))
boxplot(LENGTH.comp[,2],LENGTH.saem[,2],main='length of CI of beta1',names=c("no NA",'SAEM'))
boxplot(LENGTH.comp[,3],LENGTH.saem[,3],main='length of CI of beta2',names=c("no NA",'SAEM'))
boxplot(LENGTH.comp[,4],LENGTH.saem[,4],main='length of CI of beta3',names=c("no NA",'SAEM'))
boxplot(LENGTH.comp[,5],LENGTH.saem[,5],main='length of CI of beta4',names=c("no NA",'SAEM'))
boxplot(LENGTH.comp[,6],LENGTH.saem[,6],main='length of CI of beta5',names=c("no NA",'SAEM'))
