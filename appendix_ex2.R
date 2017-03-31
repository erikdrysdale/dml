
# Call in necessary libraries
library(magrittr); library(hdm); library(stringr); library(glmnet)

##################################################
# ------------- SIMULATION SCENARIOS ----------- #

# Set up the dgp for the logistic setting
logit.dgp <- function(pos=T,dbin=F) {
  n <- 200
  p <- 200
  px <- 10
  X <- matrix(rnorm(n*p), ncol=p)
  A <- 2
  rho <- 2
  if (pos) {
    beta <- A*c((1/(1:px)^rho),rep(0,p-px))
    gamma <- beta
  } else {
    beta <- -A*c((1/(1:px)^rho),rep(0,p-px))
    gamma <- -beta
  }
  # Keep the same intercept
  beta0 <- -1/2
  gamma0 <- -1/2
  lambda <- 1
  if (dbin) {
    P.d <- 1/(1+exp(-(gamma0 + X %*% gamma)))
    d <- rbinom(length(P.d), size=1, prob=P.d)
  } else {
  d <- gamma0 + (X %*% gamma) + rnorm(n)
  }
  P.y <- 1/(1+exp(-(beta0 + lambda*d + X %*% beta)))
  y <- rbinom(length(P.y), size=1, prob=P.y)
  ydX <- cbind(y,d,X)
  return(ydX)
}

# Number of simulations
nsim <- 2500

##################################################

# ------------- SIMULATION 1A: POSITIVE BETA WITH CONTINUOUS TREATMENT VARIABLE ----------- #

sim1a.df <- matrix(NA,nrow = nsim,ncol=4)
colnames(sim1a.df) <- c('dml1','dml2','dml3','lasso')
set.seed(1)
for (k in 1:nsim) {
  # Simulate data
  ydX <- logit.dgp(pos=T,dbin=F) %>% data.frame %>% set_colnames(c('y','d',str_c('x',1:(ncol(.)-2))))
  x <- ydX[,-c(1:2)] %>% as.matrix
  d <- as.matrix(ydX[,2], ncol = 1)
  y <- as.matrix(ydX[,1], ncol = 1)
  n <- nrow(x)
  p <- ncol(x)

  # ------ DML MODELS 1 & 2: DATA-DRIVEN LAMBDA'S ------- #
  # Find the x coefficients that are associated with y
  la1 <- 1.1/2 * sqrt(n) * qnorm(1 - 0.05/(max(n, (p + 1) * log(n))))
  l1 <- rlassologit(y ~ cbind(d,x), post = TRUE, intercept = TRUE, penalty = list(lambda.start = la1))
  xvar <-  which(l1$index[-1])
  # then x with d
  la2 <- rep(2.2 * sqrt(n) * qnorm(1 - 0.05/(max(n, p * log(n)))), p)
  l2 <- rlasso(x, d, post = TRUE, intercept = TRUE,
               penalty = list(homoscedastic = "none", lambda.start = la2, c = 1.1, gamma = 0.1))
  dvar <-  which(l2$index)
  all.x <- sort(union(xvar,dvar))
  # Now get the estimate of d
  # ---- USE THE lambda based on the robust lasso approach ----- #
  alpha1 <- coef(rlassologit(x=cbind(d,x[,all.x]),y=y,post=F,intercept=T))[2]
  alpha2 <- coef(glm(y~cbind(d,x[,all.x]),family=binomial(link='logit')))[2]

  # ------ DML MODEL 3: USE CV TO SELECT LAMBDA ------- #
  alpha.cv <- cv.glmnet(x=cbind(d,x[,all.x]),y=y,family='binomial',nfolds=10,
                         penalty.factor=c(0,rep(1,ncol(x))))
  alpha3 <- coef(glmnet(x=cbind(d,x[,all.x]),y=y,family='binomial',lambda=alpha.cv$lambda.min,
                          penalty.factor=c(0,rep(1,length(all.x)))))[2]

  # ------ LASSO MODEL 1: CV-DRIVEN LAMBDA'S ------- #
  # LASSO with no penalty on D
  lasso1.cv <- cv.glmnet(x=cbind(d,x),y=y,family='binomial',type.measure = 'class',nfolds=10,
                         penalty.factor=c(0,rep(1,ncol(x))))
  d.lasso1 <- coef(glmnet(x=cbind(d,x),y=y,family='binomial',lambda=lasso1.cv$lambda.min,
                          penalty.factor=c(0,rep(1,ncol(x)))))[2]

  # Store
  sim1a.df[k,] <- c(alpha1,alpha2,alpha3,d.lasso1)

  # Print status
  if (mod(k,50)==0) { print(k) }
}

# ------------- SIMULATION 1B: NEGATIVE BETA WITH CONTINUOUS TREATMENT VARIABLE ----------- #

sim1b.df <- matrix(NA,nrow = nsim,ncol=4)
colnames(sim1b.df) <- c('dml1','dml2','dml3','lasso')
set.seed(1)
for (k in 1:nsim) {
  # Simulate data
  ydX <- logit.dgp(pos=F,dbin=F) %>% data.frame %>% set_colnames(c('y','d',str_c('x',1:(ncol(.)-2))))
  x <- ydX[,-c(1:2)] %>% as.matrix
  d <- as.matrix(ydX[,2], ncol = 1)
  y <- as.matrix(ydX[,1], ncol = 1)
  n <- nrow(x)
  p <- ncol(x)

  # ------ DML MODELS 1 & 2: DATA-DRIVEN LAMBDA'S ------- #
  # Find the x coefficients that are associated with y
  la1 <- 1.1/2 * sqrt(n) * qnorm(1 - 0.05/(max(n, (p + 1) * log(n))))
  l1 <- rlassologit(y ~ cbind(d,x), post = TRUE, intercept = TRUE, penalty = list(lambda.start = la1))
  xvar <-  which(l1$index[-1])
  # then x with d
  la2 <- rep(2.2 * sqrt(n) * qnorm(1 - 0.05/(max(n, p * log(n)))), p)
  l2 <- rlasso(x, d, post = TRUE, intercept = TRUE,
               penalty = list(homoscedastic = "none", lambda.start = la2, c = 1.1, gamma = 0.1))
  dvar <-  which(l2$index)
  all.x <- sort(union(xvar,dvar))
  # Now get the estimate of d
  # ---- USE THE lambda based on the robust lasso approach ----- #
  alpha1 <- coef(rlassologit(x=cbind(d,x[,all.x]),y=y,post=F,intercept=T))[2]
  alpha2 <- coef(glm(y~cbind(d,x[,all.x]),family=binomial(link='logit')))[2]

  # ------ DML MODEL 3: USE CV TO SELECT LAMBDA ------- #
  alpha.cv <- cv.glmnet(x=cbind(d,x[,all.x]),y=y,family='binomial',nfolds=10,
                         penalty.factor=c(0,rep(1,ncol(x))))
  alpha3 <- coef(glmnet(x=cbind(d,x[,all.x]),y=y,family='binomial',lambda=alpha.cv$lambda.min,
                          penalty.factor=c(0,rep(1,length(all.x)))))[2]

  # ------ LASSO MODEL 1: CV-DRIVEN LAMBDA'S ------- #
  # LASSO with no penalty on D
  lasso1.cv <- cv.glmnet(x=cbind(d,x),y=y,family='binomial',type.measure = 'class',nfolds=10,
                         penalty.factor=c(0,rep(1,ncol(x))))
  d.lasso1 <- coef(glmnet(x=cbind(d,x),y=y,family='binomial',lambda=lasso1.cv$lambda.min,
                          penalty.factor=c(0,rep(1,ncol(x)))))[2]

  # Store
  sim1b.df[k,] <- c(alpha1,alpha2,alpha3,d.lasso1)

  # Print status
  if (mod(k,50)==0) { print(k) }
}

# ------------- SIMULATION 2A: POSITIVE BETA WITH BINARY TREATMENT VARIABLE ----------- #

sim2a.df <- matrix(NA,nrow = nsim,ncol=4)
colnames(sim2a.df) <- c('dml1','dml2','dml3','lasso')
set.seed(1)
for (k in 1:nsim) {
  # Simulate data
  ydX <- logit.dgp(pos=T,dbin=T) %>% data.frame %>% set_colnames(c('y','d',str_c('x',1:(ncol(.)-2))))
  x <- ydX[,-c(1:2)] %>% as.matrix
  d <- as.matrix(ydX[,2], ncol = 1)
  y <- as.matrix(ydX[,1], ncol = 1)
  n <- nrow(x)
  p <- ncol(x)

  # ------ DML MODELS 1 & 2: DATA-DRIVEN LAMBDA'S ------- #
  # Find the x coefficients that are associated with y
  la1 <- 1.1/2 * sqrt(n) * qnorm(1 - 0.05/(max(n, (p + 1) * log(n))))
  l1 <- rlassologit(y ~ cbind(d,x), post = TRUE, intercept = TRUE, penalty = list(lambda.start = la1))
  xvar <-  which(l1$index[-1])
  # then x with d
  la2 <- 1.1/2 * sqrt(n) * qnorm(1 - 0.05/(max(n, (p) * log(n))))
  l2 <- rlassologit(y ~ x, post = TRUE, intercept = TRUE, penalty = list(lambda.start = la2))
  dvar <-  which(l2$index)
  all.x <- sort(union(xvar,dvar))
  # Now get the estimate of d
  # ---- USE THE lambda based on the robust lasso approach ----- #
  alpha1 <- coef(rlassologit(x=cbind(d,x[,all.x]),y=y,post=F,intercept=T))[2]
  alpha2 <- coef(glm(y~cbind(d,x[,all.x]),family=binomial(link='logit')))[2]

  # ------ DML MODEL 3: USE CV TO SELECT LAMBDA ------- #
  alpha.cv <- cv.glmnet(x=cbind(d,x[,all.x]),y=y,family='binomial',nfolds=10,
                         penalty.factor=c(0,rep(1,ncol(x))))
  alpha3 <- coef(glmnet(x=cbind(d,x[,all.x]),y=y,family='binomial',lambda=alpha.cv$lambda.min,
                          penalty.factor=c(0,rep(1,length(all.x)))))[2]

  # ------ LASSO MODEL 1: CV-DRIVEN LAMBDA'S ------- #
  # LASSO with no penalty on D
  lasso1.cv <- cv.glmnet(x=cbind(d,x),y=y,family='binomial',type.measure = 'class',nfolds=10,
                         penalty.factor=c(0,rep(1,ncol(x))))
  d.lasso1 <- coef(glmnet(x=cbind(d,x),y=y,family='binomial',lambda=lasso1.cv$lambda.min,
                          penalty.factor=c(0,rep(1,ncol(x)))))[2]

  # Store
  sim2a.df[k,] <- c(alpha1,alpha2,alpha3,d.lasso1)

  # Print status
  if (mod(k,50)==0) { print(k) }
}

# ------------- SIMULATION 2B: NEGATIVE BETA WITH BINARY TREATMENT VARIABLE ----------- #

sim2b.df <- matrix(NA,nrow = nsim,ncol=4)
colnames(sim2b.df) <- c('dml1','dml2','dml3','lasso')
set.seed(1)
for (k in 1:nsim) {
  # Simulate data
  ydX <- logit.dgp(pos=F,dbin=T) %>% data.frame %>% set_colnames(c('y','d',str_c('x',1:(ncol(.)-2))))
  x <- ydX[,-c(1:2)] %>% as.matrix
  d <- as.matrix(ydX[,2], ncol = 1)
  y <- as.matrix(ydX[,1], ncol = 1)
  n <- nrow(x)
  p <- ncol(x)

  # ------ DML MODELS 1 & 2: DATA-DRIVEN LAMBDA'S ------- #
  # Find the x coefficients that are associated with y
  la1 <- 1.1/2 * sqrt(n) * qnorm(1 - 0.05/(max(n, (p + 1) * log(n))))
  l1 <- rlassologit(y ~ cbind(d,x), post = TRUE, intercept = TRUE, penalty = list(lambda.start = la1))
  xvar <-  which(l1$index[-1])
  # then x with d
  la2 <- 1.1/2 * sqrt(n) * qnorm(1 - 0.05/(max(n, (p) * log(n))))
  l2 <- rlassologit(y ~ x, post = TRUE, intercept = TRUE, penalty = list(lambda.start = la2))
  dvar <-  which(l2$index)
  all.x <- sort(union(xvar,dvar))
  # Now get the estimate of d
  # ---- USE THE lambda based on the robust lasso approach ----- #
  alpha1 <- coef(rlassologit(x=cbind(d,x[,all.x]),y=y,post=F,intercept=T))[2]
  alpha2 <- coef(glm(y~cbind(d,x[,all.x]),family=binomial(link='logit')))[2]

  # ------ DML MODEL 3: USE CV TO SELECT LAMBDA ------- #
  alpha.cv <- cv.glmnet(x=cbind(d,x[,all.x]),y=y,family='binomial',nfolds=10,
                         penalty.factor=c(0,rep(1,ncol(x))))
  alpha3 <- coef(glmnet(x=cbind(d,x[,all.x]),y=y,family='binomial',lambda=alpha.cv$lambda.min,
                          penalty.factor=c(0,rep(1,length(all.x)))))[2]

  # ------ LASSO MODEL 1: CV-DRIVEN LAMBDA'S ------- #
  # LASSO with no penalty on D
  lasso1.cv <- cv.glmnet(x=cbind(d,x),y=y,family='binomial',type.measure = 'class',nfolds=10,
                         penalty.factor=c(0,rep(1,ncol(x))))
  d.lasso1 <- coef(glmnet(x=cbind(d,x),y=y,family='binomial',lambda=lasso1.cv$lambda.min,
                          penalty.factor=c(0,rep(1,ncol(x)))))[2]

  # Store
  sim2b.df[k,] <- c(alpha1,alpha2,alpha3,d.lasso1)

  # Print status
  if (mod(k,50)==0) { print(k) }
}

# Save all the data in a list and write
save.list <- list(sim1a=sim1a.df,sim1b=sim1b.df,sim2a=sim2a.df,sim2b=sim2b.df)
