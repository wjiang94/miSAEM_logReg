library(MASS)
library(mvtnorm)
library(mice)
source("saem.R")

N <- 1000 # number of subjects
#N <- 10000 # number of subjects

p <- 5     # number of explanatory variables

mu.star <- 1:p  #rep(0,p)  # mean of the explanatory variables
sd <- 1:p # rep(1,p) # standard deviations

# with correlation
C <- matrix(c(   # correlation matrix
  1,   0.8, 0,   0,   0,
  0.8, 1,   0,   0,   0,
  0,   0,   1,   0.3, 0.6,
  0,   0,   0.3, 1,   0.7,
  0,   0,   0.6, 0.7, 1
), nrow=p)
## without correlation
# C = diag(p)
Sigma.star <- diag(sd)%*%C%*%diag(sd) # variance-covariance matrix of the explanatory variables

beta.star <- c(0.5, -0.3, 1, 0, -0.6) # coefficients
beta0.star <- -0.2  # intercept


nbsim = 100
#e.obs = e.predx = e.predxy = e.saem = f.predx= f.predxy = f.saem = NULL
EST.saem = EST.comp=EST.cc=EST.mice =matrix(0, nbsim,length(beta.star)+1)
TIME.saem = TIME.mice = rep(0, nbsim)
STD.saem = STD.comp=STD.mice = matrix(0, nbsim,length(beta.star)+1)
LENGTH.saem=LENGTH.comp = LENGTH.mice= matrix(0, nbsim,length(beta.star)+1)
count.saem = count.comp =  count.mice = rep(0,p+1)

#100 simus
beta.true = c(beta0.star,beta.star)
p.miss <- 0.10 
patterns = runif(N*p)<p.miss

for (NB in 1:nbsim){
  #cat(sprintf('simulation = '), NB,'\n')
  set.seed(NB)
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
  #MCAR
  X.obs <- X.complete
  X.obs[patterns] <- NA
  
  ## MAR
  # X.obs <- X.complete
  # for(i in c(2,4,5)){
  #   z <- cbind(y,X.complete[,c(1,3)])%*%matrix(sample(-5:5, 3, replace=T),ncol=1)        # linear combination 
  #   pr <- 1/(1+exp(-z))         # pass through an inv-logit function
  #   r <- rbinom(N,1,pr)      # bernoulli response variable
  #   X.obs[r==0,i]<-NA
  # }
  # cat('percentage of NA: ', mean(is.na(X.obs[,2])),mean(is.na(X.obs[,4])),mean(is.na(X.obs[,5])),'\n')
  # 
  
  # ------- estimation ignoring the missing data
  data.obs <- data.frame(y=y,X.obs)
  model.obs <- glm(y ~.,family=binomial(link='logit'),data=data.obs)
  beta0.cc <- model.obs$coefficients[1]
  beta.cc <- model.obs$coefficients[2:(p+1)]
  bias.cc = sum((c(beta0.cc,beta.cc)-c(beta0.star,beta.star))^2)
  
  
  #-multiple imputation
  ptm <- Sys.time()
  DATA.ch= cbind.data.frame(y,X.obs)
  imp.ch <- mice(DATA.ch,print = FALSE)
  fit.ch <- glm.mids(y~., data=imp.ch,family = binomial)
  #summary(fit.ch)
  beta.mice=summary(pool(fit.ch, method = "rubin1987"))[,1]
  std.mice=summary(pool(fit.ch, method = "rubin1987"))[,2]
  time.mice=Sys.time() - ptm
  
  
  
  #list.est = mcem_ar(X= X.obs,Y= y , maxruns=1000,tol_em=1e-5)
  #ptm <- proc.time()
  list.saem=miss.saem(X.obs=X.obs,y=y,maxruns=500,tol_em=1e-7,nmcmc=2,ag=1,k1=5, print_iter=FALSE,var_obs_cal=TRUE, ll_obs_cal=FALSE)
  #time_mcem_ar=proc.time() - ptm
  #time_mcem_ar=proc.time() - ptm
  beta.saem = list.saem$beta
  std.saem = list.saem$std_obs
  
  EST.comp[NB,] = c(beta0.complete,beta.complete)
  EST.cc[NB,] = c(beta0.cc,beta.cc)
  EST.saem[NB,] = beta.saem
  EST.mice[NB,] = beta.mice
  
  STD.comp[NB,] = std.complete
  STD.saem[NB,] = std.saem
  STD.mice[NB,] = std.mice
  
  TIME.saem[NB] = list.saem$time_run
  TIME.mice[NB] = time.mice
  
  ci.comp_ceil =  c(beta0.complete,beta.complete) + 1.96*std.complete
  ci.comp_ground =  c(beta0.complete,beta.complete) - 1.96*std.complete
  ci.saem_ceil = beta.saem + 1.96*std.saem
  ci.saem_ground = beta.saem - 1.96*std.saem
  ci.mice_ceil = beta.mice + 1.96*std.mice
  ci.mice_ground = beta.mice - 1.96*std.mice
  
  LENGTH.comp[NB,] = ci.comp_ceil - ci.comp_ground
  LENGTH.saem[NB,] = ci.saem_ceil - ci.saem_ground
  LENGTH.mice[NB,] = ci.mice_ceil - ci.mice_ground
  for(i in 1:(p+1)){
    if( ci.comp_ground[i] <=beta.true[i] & ci.comp_ceil[i]>beta.true[i]){
      count.comp[i]<-count.comp[i]+1
    }
    if( ci.saem_ground[i] <=beta.true[i] & ci.saem_ceil[i]>beta.true[i]){
      count.saem[i]<-count.saem[i]+1
    }
    if( ci.mice_ground[i] <=beta.true[i] & ci.mice_ceil[i]>beta.true[i]){
      count.mice[i]<-count.mice[i]+1
    }
  }
  cat('simulation =',NB,', count(complete) =',count.comp[2],', count(saem) =',count.saem[2], '\n')
}

par(mfrow=c(2,3))
boxplot(EST.comp[,1]-beta0.star,EST.cc[,1]-beta0.star,EST.mice[,1]-beta0.star,EST.saem[,1]-beta0.star,main='bias for beta0',names=c("no NA","CC","mice","saem"))
abline(h = 0, col = "red", lty = 2)
boxplot(EST.comp[,2]-beta.star[1],EST.cc[,2]-beta.star[1],EST.mice[,2]-beta.star[1],EST.saem[,2]-beta.star[1],main='bias for beta1',names=c("no NA","CC","mice","saem"))
abline(h = 0, col = "red", lty = 2)
boxplot(EST.comp[,3]-beta.star[2],EST.cc[,3]-beta.star[2],EST.mice[,3]-beta.star[2],EST.saem[,3]-beta.star[2],main='bias for beta2',names=c("no NA","CC","mice","saem"))
abline(h = 0, col = "red", lty = 2)
boxplot(EST.comp[,4]-beta.star[3],EST.cc[,4]-beta.star[3],EST.mice[,4]-beta.star[3],EST.saem[,4]-beta.star[3],main='bias for beta3',names=c("no NA","CC","mice","saem"))
abline(h = 0, col = "red", lty = 2)
boxplot(EST.comp[,5]-beta.star[4],EST.cc[,5]-beta.star[4],EST.mice[,5]-beta.star[4],EST.saem[,5]-beta.star[4],main='bias for beta4',names=c("no NA","CC","mice","saem"))
abline(h = 0, col = "red", lty = 2)
boxplot(EST.comp[,6]-beta.star[5],EST.cc[,6]-beta.star[5],EST.mice[,6]-beta.star[5],EST.saem[,6]-beta.star[5],main='bias for beta5',names=c("no NA","CC","mice","saem"))
abline(h = 0, col = "red", lty = 2)


par(mfrow=c(2,3))
boxplot(LENGTH.comp[,1],LENGTH.mice[,1],LENGTH.saem[,1],main='length of CI of beta0',names=c("no NA","mice",'SAEM'))
boxplot(LENGTH.comp[,2],LENGTH.mice[,2],LENGTH.saem[,2],main='length of CI of beta1',names=c("no NA","mice",'SAEM'))
boxplot(LENGTH.comp[,3],LENGTH.mice[,3],LENGTH.saem[,3],main='length of CI of beta2',names=c("no NA","mice",'SAEM'))
boxplot(LENGTH.comp[,4],LENGTH.mice[,4],LENGTH.saem[,4],main='length of CI of beta3',names=c("no NA","mice",'SAEM'))
boxplot(LENGTH.comp[,5],LENGTH.mice[,5],LENGTH.saem[,5],main='length of CI of beta4',names=c("no NA","mice",'SAEM'))
boxplot(LENGTH.comp[,6],LENGTH.mice[,6],LENGTH.saem[,6],main='length of CI of beta5',names=c("no NA","mice",'SAEM'))

par(mfrow=c(1,1))
boxplot(TIME.mice[1:25],TIME.saem[1:25],main='execution time',names=c('mice','SAEM'))


par(mfrow=c(2,3))
boxplot(STD.comp[,1],STD.mice[,1],STD.saem[,1],main='std for beta0',names=c("no NA","mice",'SAEM'))
abline(h = sqrt(var(EST.saem[,1])), col = "red", lty = 2)
abline(h = sqrt(var(EST.mice[,1])), col = "green", lty = 3)
abline(h = sqrt(var(EST.comp[,1])), col = "blue", lty = 4)
boxplot(STD.comp[,2],STD.mice[,2],STD.saem[,2],main='std for beta1',names=c("no NA","mice",'SAEM'))
abline(h = sqrt(var(EST.saem[,2])), col = "red", lty = 2)
abline(h = sqrt(var(EST.mice[,2])), col = "green", lty = 3)
abline(h = sqrt(var(EST.comp[,2])), col = "blue", lty = 4)
boxplot(STD.comp[,3],STD.mice[,3],STD.saem[,3],main='std for beta2',names=c("no NA","mice",'SAEM'))
abline(h = sqrt(var(EST.saem[,3])), col = "red", lty = 2)
abline(h = sqrt(var(EST.mice[,3])), col = "green", lty = 3)
abline(h = sqrt(var(EST.comp[,3])), col = "blue", lty = 4)
boxplot(STD.comp[,4],STD.mice[,4],STD.saem[,4],main='std for beta3',names=c("no NA","mice",'SAEM'))
abline(h = sqrt(var(EST.saem[,4])), col = "red", lty = 2)
abline(h = sqrt(var(EST.mice[,4])), col = "green", lty = 2)
abline(h = sqrt(var(EST.comp[,4])), col = "blue", lty = 2)
boxplot(STD.comp[,5],STD.mice[,5],STD.saem[,5],main='std for beta4',names=c("no NA","mice",'SAEM'))
abline(h = sqrt(var(EST.saem[,5])), col = "red", lty = 2)
abline(h = sqrt(var(EST.mice[,5])), col = "green", lty = 3)
abline(h = sqrt(var(EST.comp[,5])), col = "blue", lty = 4)
boxplot(STD.comp[,6],STD.mice[,6],STD.saem[,6],main='std for beta5',names=c("no NA","mice",'SAEM'))
abline(h = sqrt(var(EST.saem[,6])), col = "red", lty = 2)
abline(h = sqrt(var(EST.mice[,6])), col = "green", lty = 3)
abline(h = sqrt(var(EST.comp[,6])), col = "blue", lty = 4)



# #-----------------------bootstap------------------
# #---------------Bootstrap--------------------
# #library(boot)
# CEIL.boot.saem =GROUND.boot.saem =matrix(0, nbsim,length(beta.star)+1)
# for (nb in 21:nbsim){
#   print(nb)
#   
#   set.seed(NB)
#   # ----- complete data simulation
#   X.complete <- matrix(rnorm(N*p), nrow=N)%*%chol(Sigma.star) + matrix(rep(mu.star,N), nrow=N, byrow = TRUE)
#   p1 <- 1/(1+exp(-X.complete%*%beta.star-beta0.star))
#   y <- as.numeric(runif(N)<p1)
#   
#   # ------- generating missing data
#   X.obs <- X.complete
#   X.obs[patterns] <- NA
#   
#   #ptm <- proc.time()
#   list.saem=miss.saem(X.obs,y,maxruns=500,tol_em=1e-7,nmcmc=2,ag=1,k1=5, print_iter=FALSE,var_obs_cal=FALSE, ll_obs_cal=FALSE)
#   #time_mcem_ar=proc.time() - ptm
#   #time_mcem_ar=proc.time() - ptm
#   beta.saem = list.saem$beta
#   
#   
#   NB=100
#   betas_boot=matrix(0,NB,length(beta.star)+1)
#   for (ind_boot in 1 :NB) {
#     bootIndice=sample(1 :N, N, replace=T)
#     list.boot=miss.saem(X.obs[bootIndice,],y,maxruns=500,tol_em=1e-7,nmcmc=2,ag=1,k1=5, print_iter=FALSE,var_obs_cal=FALSE, ll_obs_cal=FALSE)
#     betas_boot[ind_boot,] <- list.boot$beta
#   }
# # delta_beta1 <- betas_boot-beta_star[2]
# # d <- quantile(delta_beta1,c(0.025,0.975))
# # ci_boot <- beta_star[2] - c(d[2],d[1])
#   for(ind_p in 1:(p+1)){
#     CEIL.boot.saem[nb,ind_p] = quantile(betas_boot[,ind_p],0.975)
#     GROUND.boot.saem[nb,ind_p] = quantile(betas_boot[,ind_p],0.025)
#   }
# }
