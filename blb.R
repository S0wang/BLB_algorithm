#######################################
## Example code for BLB Algorithm    ##
## Shuang Wang  March 16, 2021       ##
#######################################
library(glmnet)

#5 Simulations; N = 20,000; J = 100; 
#Yi =Xi+e; Xij ∼t(3),e∼N(0,10) 

sims <- 5
N <- 20000
j <- 100
Xs <- rt(n = N*j, df=3)
es <- rnorm(n=N,mean=0,sd=sqrt(10))

X_mat <- matrix(data=Xs,nrow=N,ncol=j)
X_mat <- cbind(es,X_mat)
X_mat <- cbind(apply(X_mat,1,sum),X_mat)
y = apply(X_mat,1,sum)
X = matrix(data=Xs,nrow=N,ncol=j)
df = data.frame(cbind(y,X))

# ground truth
#gt <- glmnet(y ~ X, alpha=0, lambda = 10^(-5))
gt <- lm(y ~ X)
ci <- confint(gt)
c0 <- mean(ci[,2] - ci[,1])



# ground truth
theta_gt =matrix(0,nrow=j, ncol=2000)
for(i in 1:2000){
  X <-  matrix(data=rt(n = N*j, df=3),nrow=N,ncol=j) 
  es <- rnorm(n=N,mean=0,sd=sqrt(10))
  y = apply(cbind(es,X),1,sum)
  
  gt <- lm(y ~ X-1)
  theta_gt[,i] <- as.vector(gt$coefficients)
}
c0_gt <- apply(theta_gt,1, quantile,prob=0.975) - apply(theta_gt,1, quantile,prob=0.025)
mean(c0_gt)

B=100
rr =  pt = rep(0,B) # Relative error
theta_bt =matrix(0,nrow=j, ncol=B)
Sys.time()
start.time = Sys.time()

for (b in 1:B){
  index <- sample(N, size=N, replace = T)
  newdf <- df[index, ]
  y_bt <- newdf$y
  X_bt <- as.matrix(newdf[,-1])
  f_bt<- lm(y_bt ~ X_bt -1)
  theta_bt[,b] <- f_bt$coefficients
  #CI_bt <- apply(theta_bt[,1:b],1, quantile,prob=0.975) - apply(theta_bt[,1:b],1, quantile,prob=0.025)
  #c1 <- mean( abs(CI_bt-c0)/c0)
  pt[b] <- (Sys.time() - start.time)*60
}

c1 = matrix(0,nrow=j, ncol=B)
for (b in 1:B){
  if(b==1){
    c1[,b] = rep(0,j)
  }else{
    c1[,b] = apply(theta_bt[,1:b],1, quantile,prob=0.975) - apply(theta_bt[,1:b],1, quantile,prob=0.025)
  }
  rr[b] = mean(abs(c1[,b]/mean(c0_gt) -1)) #c0_gt ground truth based on 2000 independent datasets
}
plot(pt[-1], rr[-1] )
plot(c0_gt,type="l")

#theta_bs = rep(0,B)
theta_bs = matrix(0,j,B)
for (b in 1:B){
  if(b==1){
    theta_bs[,b] = theta_bt[,b]
  }else{
    theta_bs[,b] = apply(theta_bt[,1:b],1,mean)
  }
  rr[b] = max(abs(theta_bs[,b-1]/theta_bs[,b] -1))
}

plot(pt, rr,type="l")


###***ground truth 2000 independent datasets***### ground truth
L <- 5     #number of independent datasets within each iteration. 
N <- 20000 #number of subjects
j <- 100   #number of parameters
B <- 100   #number of bootstrap iterations

theta_gt =matrix(0,nrow=j, ncol=2000)
set.seed(123)
for(i in 1:2000){
  X <-  matrix(data=rt(n = N*j, df=3),nrow=N,ncol=j) 
  es <- rnorm(n=N,mean=0,sd=sqrt(10))
  y = apply(cbind(es,X),1,sum)
  
  gt <- lm(y ~ X-1)
  theta_gt[,i] <- as.vector(gt$coefficients)
}
c0_gt <- apply(theta_gt,1, quantile,prob=0.975) - apply(theta_gt,1, quantile,prob=0.025)
mean(c0_gt)

###***simulate 5 independent datasets***###
set.seed(123)
seeds = sample(c(1:500),5)
seeds

df.array = list()
dim(df.array)
for (l in 1:L){
  set.seed(seeds[l])
  X <-  matrix(data=rt(n = N*j, df=3),nrow=N,ncol=j) 
  es <- rnorm(n=N,mean=0,sd=sqrt(10))
  y = apply(cbind(es,X),1,sum)
  df.array[[l]] = data.frame(cbind(y,X))
}

theta_bt_each = matrix(0,j,L)
pt_each =rep(0,L)
for (b in 1:B){
  for (l in 1:L){
    start.time = Sys.time()
    index <- sample(N, size=N, replace = T)
    newdf <- df.array[[l]][index, ]
    y_bt <- newdf$y
    X_bt <- as.matrix(newdf[,-1])
    f_bt<- lm(y_bt ~ X_bt -1)
    theta_bt_each[,l] <- f_bt$coefficients
    pt_each[l] <- (Sys.time() - start.time) 
  }
  theta_bt[,b] <- apply(theta_bt_each,1,mean)
  pt[b] <- mean(pt_each)
}

# CI width for calculating rr
pt.cum = rr.ciw = rep(0,B)
c1 = matrix(0,nrow=j, ncol=B)
for (b in 1:B){
  if(b==1){
    c1[,b] = rep(0,j)
    pt.cum[b] = pt[b]
  }else{
    c1[,b] = apply(theta_bt[,1:b],1, quantile,prob=0.975) - apply(theta_bt[,1:b],1, quantile,prob=0.025)
    pt.cum[b] = pt[b] + pt_cum[b-1]
  }
  rr.ciw[b] = mean(abs(c1[,b]/c0_gt -1))  ##c0_gt ground truth based on 2000 independent datasets
}
plot(pt.cum, rr.ciw,type="l")
plot(c0_gt,type="l")

#Dr. Wu's method for calculating rr
pt_cum = rr.theta = rep(0,B)
theta_bs = matrix(0,j,B)
for (b in 1:B){
  if(b==1){
    theta_bs[,b] = theta_bt[,b]
    pt_cum[b] = pt[b]
  }else{
    theta_bs[,b] = apply(theta_bt[,1:b],1,mean)
    pt_cum[b] = pt[b] + pt_cum[b-1]
  }
  rr.theta[b] = max(abs(theta_bs[,b-1]/theta_bs[,b] -1))
}
plot(pt_cum, rr.theta,type="l")


# correct
TT=5
PT = RR = THETA_bt= rep(0,TT)
THETA_all = array(0,dim=c(dim(theta_bt)[1],dim(theta_bt)[2],20))
#THETA_bt = matrix(0, nrow=j,ncol=TT)
for(tt in 1:TT){
  X_mat <- matrix(data=Xs,nrow=N,ncol=j)
  X_mat <- cbind(es,X_mat)
  X_mat <- cbind(apply(X_mat,1,sum),X_mat)
  y = apply(X_mat,1,sum)
  X = matrix(data=Xs,nrow=N,ncol=j)
  df = data.frame(cbind(y,X))
  
  start.time = Sys.time()
  for (b in 1:B){
    index <- sample(N, size=N, replace = T)
    newdf <- df[index, ]
    y_bt <- newdf$y
    X_bt <- as.matrix(newdf[,-1])
    f_bt<- lm(y_bt ~ X_bt -1)
    theta_bt[,b] <- f_bt$coefficients
  }
  theta_tmp=mean(apply(theta_bt,1,mean))
  THETA_bt[tt] = theta_tmp # bootstrap estimate
  #ci_bt <- apply(theta_bt,1, quantile,prob=0.975) - apply(theta_bt,1, quantile,prob=0.025)
  RR[tt] <- max(abs(THETA_bt[tt-1]/theta_tmp-1)) #mean(abs(ci_bt-ci_bt[tt-1])/ci_bt[tt-1])
  PT[tt] <- (Sys.time() - start.time) 
}
pt.x = rep(0,TT)
PT.tmp=0
for(tt in 1:TT){
  PT.tmp=PT.tmp+PT[tt]
  pt.x[tt] = PT.tmp
}

plot(pt.x,RR)



