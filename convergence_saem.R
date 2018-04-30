library(MASS)
library(mvtnorm)
library(ggplot2)
library(reshape2)
library(dplyr)
library(RColorBrewer)
theme_set(theme_bw())
source("saem.R")

N <- 1000  # number of subjects
p <- 5     # number of explanatory variables

mu.star <- 1:p  #rep(0,p)  # mean of the explanatory variables
sd <- 1:p # rep(1,p) # standard deviations

C <- matrix(c(   # correlation matrix
  1,   0.8, 0,   0,   0,
  0.8, 1,   0,   0,   0,
  0,   0,   1,   0.3, 0.6,
  0,   0,   0.3, 1,   0.7,
  0,   0,   0.6, 0.7, 1
), nrow=p)

Sigma.star <- diag(sd)%*%C%*%diag(sd) # variance-covariance matrix of the explanatory variables

beta.star <- c(0.5, -0.3, 1, 0, -0.6) # coefficients
beta0.star <- -0.2  # intercept

beta.true = c(beta0.star,beta.star)
p.miss <- 0.10 
patterns = runif(N*p)<p.miss

NB = 5
ag <- c(0.5, 0.7, 1)
k1 <- 50
maxruns=300

BIASBETA1_0.5 = BETA1_0.5 = matrix(0, NB, maxruns+1)
BIASBETA1_0.7 = BETA1_0.7 = matrix(0, NB, maxruns+1)
BIASBETA1_1.0 = BETA1_1.0 = matrix(0, NB, maxruns+1)
BIASBETAAVG1_0.5 = BETAAVG1_0.5 = matrix(0, NB, maxruns+1)
BIASBETAAVG1_0.7 = BETAAVG1_0.7 = matrix(0, NB, maxruns+1)
BIASBETAAVG1_1.0 = BETAAVG1_1.0 = matrix(0, NB, maxruns+1)

seed <- c(1,2,4,6,8,10,12,14,16,18)

for(nb in 1:NB){
  set.seed(seed[nb])  
  print(nb)
  # ----- complete data simulation
  X.complete <- matrix(rnorm(N*p), nrow=N)%*%chol(Sigma.star) + matrix(rep(mu.star,N), nrow=N, byrow = TRUE)
  p1 <- 1/(1+exp(-X.complete%*%beta.star-beta0.star))
  y <- as.numeric(runif(N)<p1)
  
  data.complete <- data.frame(y=y,X.complete)
  model.complete <- glm(y ~.,family=binomial(link='logit'),data=data.complete)
  beta.complete <- model.complete$coefficients
  # ------- generating missing data
  X.obs <- X.complete
  # for(i in c(2,4,5)){
  #   z <- cbind(y,X.complete[,c(1,3)])%*%matrix(sample(1:10, 3, replace=T),ncol=1)        # linear combination 
  #   pr <- 1/(1+exp(-z))         # pass through an inv-logit function
  #   r <- rbinom(N,1,pr)      # bernoulli response variable
  #   X.obs[r==0,i]<-NA
  # }
  X.obs[patterns] <- NA
  #cat('percentage of NA: ', mean(is.na(X.obs[,2])),mean(is.na(X.obs[,4])),mean(is.na(X.obs[,5])),'\n')

  #   list.saem0.5=miss.saem(X.obs,y,maxruns=maxruns,tol_em=1e-50,nmcmc=1,ag=ag[1],k1=k1,print_iter=FALSE)
  # BIASBETA1_0.5[nb,] = list.saem0.5$seqbeta[2,] 
  # BIASBETAAVG1_0.5[nb,] = list.saem0.5$seqbeta_avg[2,] 
  
  list.saem0.7=miss.saem(X.obs=X.obs,y=y,maxruns=maxruns,tol_em=1e-50,nmcmc=5,ag=ag[2],k1=k1,print_iter=FALSE)
  BETA1_0.7[nb,] = list.saem0.7$seqbeta[2,]
  BETAAVG1_0.7[nb,] = list.saem0.7$seqbeta_avg[2,]
  BIASBETA1_0.7[nb,] = list.saem0.7$seqbeta[2,] - list.saem0.7$beta[2]
  BIASBETAAVG1_0.7[nb,] = list.saem0.7$seqbeta_avg[2,] -list.saem0.7$seqbeta_avg[2,maxruns]
  
  
  list.saem1.0=miss.saem(X.obs=X.obs,y=y,maxruns=maxruns,tol_em=1e-50,nmcmc=5,ag=ag[3],k1=k1,print_iter=FALSE)
  BETA1_1.0[nb,] = list.saem1.0$seqbeta[2,]
  BETAAVG1_1.0[nb,] = list.saem1.0$seqbeta_avg[2,] 
  BIASBETA1_1.0[nb,] = list.saem1.0$seqbeta[2,] -list.saem1.0$beta[2]
  BIASBETAAVG1_1.0[nb,] = list.saem1.0$seqbeta_avg[2,] -list.saem1.0$seqbeta_avg[2,maxruns]
}

#----------------------------------------------------

fnames <- c("0.7", "1.0", "0.7Av")
df1 <- as.data.frame(t(BIASBETA1_0.7))
names(df1) <- 1:NB
df1['iteration'] <- 0:(nrow(df1)-1)
df1 <- melt(df1, variable.name="replicate", id.vars = list("iteration")) 
df1['tau'] = fnames[1] 
df2 <- as.data.frame(t(BIASBETA1_1.0))
names(df2) <- 1:NB
df2['iteration'] <- 0:(nrow(df2)-1)
df2 <- melt(df2, variable.name="replicate", id.vars = list("iteration")) 
df2['tau'] = fnames[2] 
df3 <- as.data.frame(t(BIASBETAAVG1_0.7))
names(df3) <- 1:NB
df3['iteration'] <- 0:(nrow(df3)-1)
df3 <- melt(df3, variable.name="replicate", id.vars = list("iteration")) 
df3['tau'] = fnames[3] 

df <- rbind(df1, df2, df3)
df[['tau']] <- factor(df[['tau']], levels=fnames)
levels(df[['tau']]) <- c("tau*' = 0.7'", "tau*' = 1'", "tau*'= 0.7 + Averaging'")

beta2 <- subset(df, iteration==maxruns)
beta1 <- beta2
beta1$iteration <- 0
beta <- rbind(beta1, beta2)

pl <- ggplot(df) + geom_line(aes(iteration,value,color=replicate)) + 
  geom_line(data=beta, aes(iteration, value, color=replicate), linetype=3) +
  facet_grid(~tau, labeller = label_parsed) +  ylab(expression(bias~of~beta[1])) +
  theme(strip.text = element_text(size=12), axis.title=element_text(size=14), 
        legend.position="none")
print(pl)

#---------------------------------------------------------------------------

df1 <- as.data.frame(t(list.saem0.7$seqbeta))
names(df1) <- paste0("beta[",0:5,"]")
df1['iteration'] <- 0:(nrow(df1)-1)
df1 <- melt(df1, variable.name="parameter", id.vars = list("iteration")) 
df1['tau'] = fnames[1] 
df2 <- as.data.frame(t(list.saem1.0$seqbeta))
names(df2) <- paste0("beta[",0:5,"]")
df2['iteration'] <- 0:(nrow(df2)-1)
df2 <- melt(df2, variable.name="parameter", id.vars = list("iteration")) 
df2['tau'] = fnames[2] 
df3 <- as.data.frame(t(list.saem0.7$seqbeta_avg))
names(df3) <- paste0("beta[",0:5,"]")
df3['iteration'] <- 0:(nrow(df3)-1)
df3 <- melt(df3, variable.name="parameter", id.vars = list("iteration")) 
df3['tau'] = fnames[3] 

df <- rbind(df1, df2, df3)
df[['tau']] <- factor(df[['tau']], levels=fnames)
levels(df[['tau']]) <- c("tau*' = 0.7'", "tau*' = 1'", "tau*'= 0.7 + Averaging'")

beta2 <- subset(df, iteration==maxruns)
beta1 <- beta2
beta1$iteration <- 0
beta <- rbind(beta1, beta2)

ldf <- levels(df$parameter)
labl <- list(expression(beta[0]), expression(beta[1]), expression(beta[2]),
             expression(beta[3]), expression(beta[4]), expression(beta[5]) ) 

palette(brewer.pal(6, "Dark2"))
pl <- ggplot(df) + geom_line(aes(iteration,value,color=parameter)) + 
#  geom_line(data=beta, aes(iteration, value, color=replicate)) +
  facet_grid(~tau, labeller = label_parsed) +  ylab(expression(beta)) +
  scale_color_manual(labels = labl, values=1:6) +
  theme(strip.text = element_text(size=12), axis.title=element_text(size=14))
print(pl)


#---------------------------------------------------------------

fnames <- c("0.7", "1.0", "0.7Av")
df1 <- as.data.frame(t(BETA1_0.7))
names(df1) <- 1:NB
df1['iteration'] <- 0:(nrow(df1)-1)
df1 <- melt(df1, variable.name="replicate", id.vars = list("iteration")) 
df1['tau'] = fnames[1] 
df2 <- as.data.frame(t(BETA1_1.0))
names(df2) <- 1:NB
df2['iteration'] <- 0:(nrow(df2)-1)
df2 <- melt(df2, variable.name="replicate", id.vars = list("iteration")) 
df2['tau'] = fnames[2] 
df3 <- as.data.frame(t(BETAAVG1_0.7))
names(df3) <- 1:NB
df3['iteration'] <- 0:(nrow(df3)-1)
df3 <- melt(df3, variable.name="replicate", id.vars = list("iteration")) 
df3['tau'] = fnames[3] 

df <- rbind(df1, df2, df3)
df[['tau']] <- factor(df[['tau']], levels=fnames)
levels(df[['tau']]) <- c("tau*' = 0.7'", "tau*' = 1'", "tau*'= 0.7 + Averaging'")

beta2 <- subset(df, iteration==maxruns)
beta1 <- beta2
beta1$iteration <- 0
beta <- rbind(beta1, beta2)

pl <- ggplot(df) + geom_line(aes(iteration,value,color=replicate)) + 
  geom_line(data=beta, aes(iteration, value, color=replicate), linetype=3) +
  facet_grid(~tau, labeller = label_parsed) +  ylab(expression(beta[1])) +
  theme(strip.text = element_text(size=12), axis.title=element_text(size=14), 
        legend.position="none")
print(pl)
