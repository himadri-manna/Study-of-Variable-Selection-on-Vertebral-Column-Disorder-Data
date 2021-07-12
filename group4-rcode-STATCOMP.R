library(truncnorm)
library(bayesplot)
library(coda)
library(pracma)
library(mvtnorm)
library(ggplot2)
library(ggpubr)
library(GGally)
library(MASS)


data.set <- read.csv("data.csv")
train<-data.set

for (i in 1:nrow(train)) {
  if(train$class[i]==1) {
    train$class[i]='Abnormal'
  }else{train$class[i]='Normal'}
}
s1<-ggplot(train,aes(x=pelvic_incidence,fill=class))+geom_histogram(alpha=0.5, aes(y=..density..,fill=factor(class))) + labs(x='Pelvic Incidence')  + theme_grey()
s2<-ggplot(train,aes(x=pelvic_tilt,fill=class))+geom_histogram(alpha=0.5, aes(y=..density..,fill=factor(class))) + labs(x="Pelvic Tilt")  + theme_grey()
s3<-ggplot(train,aes(x=lumbar_lordosis_angle,fill=class))+geom_histogram(alpha=0.5, aes(y=..density..,fill=factor(class))) + labs(x="Lumbar Lordosis angle") + theme_grey()
s4<-ggplot(train,aes(x=sacral_slope,fill=class))+geom_histogram(alpha=0.5, aes(y=..density..,fill=factor(class))) + labs(x="Sacral Slope") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) + theme_grey()
s5<-ggplot(train,aes(x=pelvic_readius,fill=class))+geom_histogram(alpha=0.5, aes(y=..density..,fill=factor(class))) + labs(x="Pelvic Radius") + theme_grey()
s6<-ggplot(train,aes(x= grade_of_spondylolisthesis,fill=class))+geom_histogram(alpha=0.5, aes(y=..density..,fill=factor(class))) + labs(x="Grade of Spondylolisthesis") + theme_grey()

ggarrange(s1,s2,s3,s4,s5,s6)

ggpairs(train[,2:7],aes(color = train$class))

X.full <- as.matrix(data.set[-c(116),-1])
X.full[,-7] <- scale(X.full[,-7])
X.full <- cbind(rep(1,dim(X.full)[1]), X.full)

set.seed(743)
sam <- sort(sample(1:dim(X.full)[1], 240, replace = FALSE))

X <- X.full[sam, -8]
X.test<-X.full[-sam,-c(1,8)]
y <- X.full[sam, 8]
X1<-X[,-1]
Xy<-cbind(as.data.frame(X1),y=y)

full<-glm(y~pelvic_incidence+pelvic_tilt+lumbar_lordosis_angle+sacral_slope+pelvic_readius+grade_of_spondylolisthesis,data = Xy,family = binomial(link = "probit"))
backward<-stepAIC(full,direction = 'backward',k=log(nrow(X1)))
backward_model<-glm(y~sacral_slope+pelvic_readius+grade_of_spondylolisthesis,data = Xy,family = binomial(link = "probit"))
summary(backward_model)
# Intercept   beta_4   beta_5   beta_6  
# 1.6009      -0.8726  -0.9249  2.6679
backward_train <- predict.glm(backward_model, newdata = as.data.frame(X1), type = "response")
for (i in 1:length(backward_train)) {
  if(backward_train[i]>=0.5) {
    backward_train[i]=1
  }else{backward_train[i]=0}
  
}
actual_class_train<-X.full[sam, 8]
table(backward_train,actual_class_train)
#training error rate=13.33%
backward_test <- predict.glm(backward_model, newdata = as.data.frame(X.test), type = "response")
for (i in 1:length(backward_test)) {
  if(backward_test[i]>=0.5) {
    backward_test[i]=1
  }else{backward_test[i]=0}
  
}
actual_class<-X.full[-sam, 8]
table(backward_test,actual_class)
#test error rate=18.84%

backward2<-stepAIC(full,direction = 'backward')
backward_model2<-glm(y~pelvic_incidence + sacral_slope + pelvic_readius + grade_of_spondylolisthesis,data = Xy,family = binomial(link = "probit"))
summary(backward_model2)
# Intercept    beta_3    beta_4    beta_5   beta_6  
# 1.6528       0.5596   -1.2673   -0.8208   2.6359
backward_train2 <- predict.glm(backward_model2, newdata = as.data.frame(X1), type = "response")
for (i in 1:length(backward_train2)) {
  if(backward_train2[i]>=0.5) {
    backward_train2[i]=1
  }else{backward_train2[i]=0}
  
}
actual_class_train<-X.full[sam, 8]
table(backward_train2,actual_class_train)
#training error rate=13.75%
backward_test2 <- predict.glm(backward_model2, newdata = as.data.frame(X.test), type = "response")
for (i in 1:length(backward_test2)) {
  if(backward_test2[i]>=0.5) {
    backward_test2[i]=1
  }else{backward_test2[i]=0}
  
}
actual_class<-X.full[-sam, 8]
table(backward_test2,actual_class)
#test error rate=17.39%




gibbs.mod <- function(y, X, rep) {
  n <- dim(X)[1]
  p <- dim(X)[2]
  
  beta_0 <- rep(0, p)
  sigma_0 <- diag(0, nrow = p)
  diag(sigma_0) <- rep(1000, p)
  
  df_beta <- matrix(NA, ncol = p, nrow = rep)
  df_gamma <- matrix(NA, ncol = p, nrow = rep)
  
  z <- rep(1:n)
  gamma <-  c(1, rbinom(p-1,1,prob = 0.5))
  beta <-  c(1:p)
  
  
  for (t in 1:rep) {
    v <- matrix(gamma*beta)
    for (i in 1:n) {
      if(y[i]==1) {
        z[i] <- rtruncnorm(1, mean = t(as.matrix(X[i,])) %*% v, 
                           sd=1, a=0, b=Inf)
      }
      else  {
        z[i] <- rtruncnorm(1, mean = t(as.matrix(X[i,])) %*% v, 
                           sd=1, b=0, a=-Inf)
      }
    }
    
    X_star <- sweep( as.matrix(X), MARGIN=2, as.matrix(gamma),'*')
    b_star_var <- inv( inv(sigma_0) + t(X_star)%*%(X_star) )
    b_star <- b_star_var %*% ( t(X_star)%*% z) # inv(sigma_0)%*% beta_0 +
    beta <- rmvnorm( n=1, mean = b_star, sigma = b_star_var)
    
    g1 <-  function(m){
      gamma[m] = 1
      v_s <- gamma*beta
      val <- exp(-0.5* t(z- X %*% t(v_s)) %*%
                   (z- X %*% t(v_s)))
      return(val)
    }
    g0 <-  function(m){
      gamma[m] = 0
      v_s <- gamma*beta
      val <- exp(-0.5* t(z- X %*% t(v_s)) %*%
                   (z- X %*% t(v_s)))
      return(val)
    }
    
    for (j in 2:p){
      if(g0(j)==0) pi.c.j = 1
      else pi.c.j = (g1(j))/(g1(j)+g0(j))
      gamma[j] = rbinom(1,1,prob = pi.c.j)
    }
    df_beta[t,] <-  beta
    df_gamma[t,] <-  gamma
    if(t %% 10^4 == 0) message(t)
  }
  return(list("Coeff" = df_beta, "inclusion" = df_gamma))
}

thin <-  200
burn <-  10^4
rep <- 510000

set.seed(1234)
chain <- gibbs.mod(y, X, rep )

gamma <- as.matrix(chain[[2]])
colnames(gamma) <-  c("Intercept", "gamma_1", "gamma_2", "gamma_3", "gamma_4", 
                      "gamma_5", "gamma_6")
gamma <- gamma[seq(from = burn+1, to = rep, by = thin),]
colMeans(gamma)

# Intercept   gamma_1   gamma_2   gamma_3   gamma_4   gamma_5   gamma_6 
# 1.0000    0.1156    0.0896    0.0076    0.9564    1.0000    1.0000 

gamma.u <- unique(gamma)

model.post <- matrix(NA, nrow = nrow(gamma), ncol = nrow(gamma.u))
for (j in 1:nrow(gamma.u)) {
  for( i in 1:nrow(gamma)){
    model.post[i,j] <- ifelse( identical(gamma[i,], gamma.u[j,]),
                               1, 0)
  }
}
colMeans(model.post)
#WOO 0.8552 0.0220 0.0480 0.0004 0.0004 0.0244 0.0424 0.0060 0.0004 0.0004 0.0004

gibbs.mod.beta <- function(y, X, rep) {
  n <- dim(X)[1]
  p <- dim(X)[2]
  
  beta_0 <- rep(0, p)
  sigma_0 <- diag(0, nrow = p)
  diag(sigma_0) <- rep(1000, p)
  
  df_beta <- matrix(NA, ncol = p, nrow = rep)
  z <- rep(1:n)
  beta <-  c(1:p)
  
  for (t in 1:rep) {
    for (i in 1:n) {
      if(y[i]==1) {
        z[i] <- rtruncnorm(1, mean = t(as.matrix(X[i,])) %*% matrix(beta, nrow = length(beta)),
                           sd=1, a=0, b=Inf)
      }
      else  {
        z[i] <- rtruncnorm(1, mean = t(as.matrix(X[i,])) %*% matrix(beta, nrow = length(beta)),
                           sd=1, b=0, a=-Inf)
      }
    }
    
    b_star_var <- inv( inv(sigma_0) + t(X) %*% (X) )
    b_star <- b_star_var %*% ( t(X)%*% z) # inv(sigma_0)%*% beta_0 +
    beta <- rmvnorm( n=1, mean = b_star, sigma = b_star_var)
    
    df_beta[t,] <-  beta
    if(t %% 10^4 == 0) message(t)
  }
  return(list("Coeff" = df_beta))
}

X.M <- X[,c(1, 5, 6, 7)]
set.seed(1234)
chain.b <- gibbs.mod.beta(y, X.M, rep )

beta.df <- data.frame(chain.b[[1]])
colnames(beta.df) = c("beta_0", "beta_4", "beta_5", "beta_6")

beta.df <-beta.df[seq(from = burn+1, to = rep, by = 200),]
beta.mcmc <- as.mcmc(beta.df)
mcmc_acf(beta.mcmc)
mcmc_trace(beta.mcmc)
colMeans(beta.df)
# beta_0     beta_4     beta_5     beta_6 
# 1.6681182 -0.8885347 -0.9437431  2.7654438 

apply(beta.df, 2, median)
# beta_0     beta_4     beta_5     beta_6 
# 1.6562169 -0.8906949 -0.9438566  2.7527849 

beta.ci <- HPDinterval(beta.mcmc, prob = 0.95)
beta.ci
#           lower      upper
# beta_0  1.222977  2.1860137
# beta_4 -1.255135 -0.5369650
# beta_5 -1.296979 -0.6201231
# beta_6  1.969314  3.5006822


X.M.train <- as.matrix(X.full[sam, c(1, 5, 6, 7)])
z.pred.train <- X.M.train %*% colMeans(beta.df)
y.pred.train <- ifelse( z.pred.train > 0, 1, 0)

table(Predicted = y.pred.train, Actual = X.full[sam, 8])
mean( X.full[sam, 8] != y.pred.train) 
#             Actual
# Predicted   0   1
#         0  58  18
#         1  14 150
# 0.1333333

X.M.test <- as.matrix(X.full[-sam, c(1, 5, 6, 7)])
z.pred <- X.M.test %*% colMeans(beta.df)
y.pred <- ifelse( z.pred > 0, 1, 0)

table(Predicted = y.pred, Actual = X.full[-sam, 8])
mean( X.full[-sam, 8] != y.pred)  
#           Actual
# Predicted  0  1
#         0 18  4
#         1  9 38
# 0.1884058

pdf("plot-ab.pdf")
mcmc_trace(beta.mcmc)
par(mfrow=c(2,2))
estint <- cumsum(beta.df[,1])/(1:nrow(beta.df)) 
esterr <- sqrt(cumsum((beta.df[,1]-estint)^2))/(1:nrow(beta.df)) 
plot(estint, xlab = "", type = "l", lwd= + 2,
     ylim = mean(beta.df[,1]) + 20*c(-esterr[nrow(beta.df)],esterr[nrow(beta.df)]),
     ylab="Mean and error range", main = expression(beta[0])) 
lines(estint+2*esterr,col="gold",lwd=2) 
lines(estint-2*esterr,col="gold",lwd=2)

estint <- cumsum(beta.df[,2])/(1:nrow(beta.df)) 
esterr <- sqrt(cumsum((beta.df[,2]-estint)^2))/(1:nrow(beta.df)) 
plot(estint, xlab = "", type = "l", lwd= + 2,
     ylim = mean(beta.df[,2]) + 20*c(-esterr[nrow(beta.df)],esterr[nrow(beta.df)]),
     ylab="Mean and error range", main = expression(beta[4])) 
lines(estint+2*esterr,col="gold",lwd=2) 
lines(estint-2*esterr,col="gold",lwd=2)

estint <- cumsum(beta.df[,3])/(1:nrow(beta.df)) 
esterr <- sqrt(cumsum((beta.df[,3]-estint)^2))/(1:nrow(beta.df)) 
plot(estint, xlab = "", type = "l", lwd= + 2,
     ylim = mean(beta.df[,3]) + 20*c(-esterr[nrow(beta.df)],esterr[nrow(beta.df)]),
     ylab="Mean and error range", main = expression(beta[5])) 
lines(estint+2*esterr,col="gold",lwd=2) 
lines(estint-2*esterr,col="gold",lwd=2)

estint <- cumsum(beta.df[,4])/(1:nrow(beta.df)) 
esterr <- sqrt(cumsum((beta.df[,4]-estint)^2))/(1:nrow(beta.df)) 
plot(estint, xlab = "", type = "l", lwd= + 2,
     ylim = mean(beta.df[,4]) + 20*c(-esterr[nrow(beta.df)],esterr[nrow(beta.df)]),
     ylab="Mean and error range", main = expression(beta[6])) 
lines(estint+2*esterr,col="gold",lwd=2) 
lines(estint-2*esterr,col="gold",lwd=2)
dev.off()

################# Model Selection using Bayesian Lasso
gibbs.lasso=function(y, X, rep, a, b, burn){
  n = dim(X)[1]
  p = dim(X)[2]
  
  
  betadraw = matrix(0,nrow=rep,ncol=p)
  
  
  beta = rep(1,p)
  lambda2 = 1
  tau2 = rep(1,p)
  z = rep(0,n)
  for (i in 1:n) {
    if(y[i]==1) z[i] = 1
    else z[i] = -1
  }
  
  
  for (iter in 1:rep) {
    lambda2 = rgamma(1, p+a, b+sum(tau2)/2)
    
    shape = lambda2
    mean =  sqrt(shape/(beta^2))
    tau2 = c(1/rinvgauss(p, mean = mean, shape = shape))
    
    
    XtX = t(X) %*% X
    beta.cov = inv(XtX + diag(1/tau2, nrow = p))
    Xtz = t(X) %*% z
    beta.mean = beta.cov %*% Xtz
    beta = rmvnorm(1, mean = beta.mean, sigma = beta.cov)
    
    z.mean = X %*% matrix(beta, nrow = length(beta))
    for (i in 1:n) {
      z[i] = ifelse(y[i]==1, rtruncnorm(1, a=0, b=Inf, mean = z.mean[i], sd=1), 
                    rtruncnorm(1, a=-Inf, b=0, mean = z.mean[i], sd=1))
    }
    
    betadraw[iter,] = beta
    if(iter %% 10000 == 0) message(iter)
  }
  sample = list("Coefficients" = betadraw[-seq(burn),])
  return(sample)
}





thin <-  120
burn <-  10^4
rep <- 2.5*10^5


set.seed(1234)
chain <- gibbs.lasso(y, X, rep, 0.5, 0.5, burn)



beta.df <- data.frame(chain[[1]])
colnames(beta.df) = c("beta_0", "beta_1", "beta_2", "beta_3", "beta_4", "beta_5", "beta_6")

cand = seq(1, nrow(beta.df), by=thin)
beta.df <-beta.df[cand,]
beta.mcmc <- as.mcmc(beta.df)
mcmc_acf(beta.mcmc)
mcmc_trace(beta.mcmc)
colMeans(beta.df)
# beta_0      beta_1      beta_2      beta_3      beta_4      beta_5      beta_6 
#1.56276808 -0.09035723  0.40756949 -0.14285184 -0.61122904 -0.76579466  2.49977260  
apply(beta.df, 2, median)
#   beta_0     beta_1     beta_2     beta_3     beta_4     beta_5     beta_6 
#1.55055256 -0.08378184  0.39595142 -0.13514894 -0.60033051 -0.76154425  2.49646424     
beta.ci <- HPDinterval(beta.mcmc, prob = 0.95)
beta.ci
#              lower      upper
#beta_0     1.1224095  2.0635469
#beta_1    -1.6965093  1.7962278
#beta_2    -0.6344210  1.4592239
#beta_3    -0.6051095  0.3375358
#beta_4    -1.9786189  0.6408568
#beta_5    -1.1066248 -0.4304665
#beta_6     1.7320914  3.2665078









length = (rep-burn)/thin
par(mfrow=c(2,2))
int=beta.df[,1]
estint=cumsum(int)/(1:length)
esterr=sqrt(cumsum((int-estint)^2))/(1:length) 
plot(estint, ylab="Mean and error range",type="l",lwd= + 2,ylim=mean(int)+20*c(-esterr[length],esterr[length]),
     xlab="",main=expression(beta[0])) 
lines(estint+2*esterr,col="gold",lwd=2) 
lines(estint-2*esterr,col="gold",lwd=2)

beta1=beta.df[,2]
estint=cumsum(beta1)/(1:length)
esterr=sqrt(cumsum((beta1-estint)^2))/(1:length) 
plot(estint, ylab="Mean and error range",type="l",lwd= + 2,ylim=mean(beta1)+20*c(-esterr[length],esterr[length]),
     xlab="",main=expression(beta[1])) 
lines(estint+2*esterr,col="gold",lwd=2) 
lines(estint-2*esterr,col="gold",lwd=2)

beta2=beta.df[,3]
estint=cumsum(beta2)/(1:length)
esterr=sqrt(cumsum((beta2-estint)^2))/(1:length) 
plot(estint, ylab="Mean and error range",type="l",lwd= + 2,ylim=mean(beta2)+20*c(-esterr[length],esterr[length]),
     xlab="",main=expression(beta[2])) 
lines(estint+2*esterr,col="gold",lwd=2) 
lines(estint-2*esterr,col="gold",lwd=2)

beta3=beta.df[,4]
estint=cumsum(beta3)/(1:length)
esterr=sqrt(cumsum((beta3-estint)^2))/(1:length) 
plot(estint, ylab="Mean and error range",type="l",lwd= + 2,ylim=mean(beta3)+20*c(-esterr[length],esterr[length]),
     xlab="",main=expression(beta[3])) 
lines(estint+2*esterr,col="gold",lwd=2) 
lines(estint-2*esterr,col="gold",lwd=2)

beta4=beta.df[,5]
estint=cumsum(beta4)/(1:length)
esterr=sqrt(cumsum((beta4-estint)^2))/(1:length) 
plot(estint, ylab="Mean and error range",type="l",lwd= + 2,ylim=mean(beta4)+20*c(-esterr[length],esterr[length]),
     xlab="",main=expression(beta[4])) 
lines(estint+2*esterr,col="gold",lwd=2) 
lines(estint-2*esterr,col="gold",lwd=2)

beta5=beta.df[,6]
estint=cumsum(beta5)/(1:length)
esterr=sqrt(cumsum((beta5-estint)^2))/(1:length) 
plot(estint, ylab="Mean and error range",type="l",lwd= + 2,ylim=mean(beta5)+20*c(-esterr[length],esterr[length]),
     xlab="",main=expression(beta[5])) 
lines(estint+2*esterr,col="gold",lwd=2) 
lines(estint-2*esterr,col="gold",lwd=2)

beta6=beta.df[,7]
estint=cumsum(beta6)/(1:length)
esterr=sqrt(cumsum((beta6-estint)^2))/(1:length) 
plot(estint, ylab="Mean and error range",type="l",lwd= + 2,ylim=mean(beta6)+20*c(-esterr[length],esterr[length]),
     xlab="",main=expression(beta[6])) 
lines(estint+2*esterr,col="gold",lwd=2) 
lines(estint-2*esterr,col="gold",lwd=2)








beta.hat = rep(0, dim(X)[2])
for (j in 1:dim(X)[2]) {
  if(beta.ci[j,1]>0 || beta.ci[j,2]<0) beta.hat[j] = colMeans(beta.df)[j]
}






X.test <- as.matrix(X.full[-sam, -8])

z.pred.test <- X.test %*% as.matrix(beta.hat, nrow = length(beta.hat))
p.pred.test <- pnorm(z.pred.test, lower.tail = T)

y.pred.test <- ifelse( p.pred.test >= 0.5, 1, 0)

Z.pred.train <- X %*% as.matrix(beta.hat, nrow = length(beta.hat))
p.pred.train <- pnorm(Z.pred.train, lower.tail = T)

y.pred.train <- ifelse( p.pred.train >= 0.5, 1, 0)


table(Actual = X.full[-sam, 8], Predicted = y.pred.test)
#      Predicted
#Actual   0   1
#     0  21  6
#     1   7 35
mean( X.full[-sam, 8] != y.pred.test) # 0.1884058 -> test error


table(Actual = X.full[sam, 8], Predicted = y.pred.train)
#      Predicted
#Actual   0   1
#     0  60  12
#     1  38 130
mean( X.full[sam, 8] != y.pred.train) # 0.2083333 -> training error
