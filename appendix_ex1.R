#########################################################################
############# ----------- ML CAUSALITY SCRIPT 2 ----------- #############
#########################################################################

rm(list=ls())

# Call in the CRAN packages
ll <- c('tidyverse','magrittr','cowplot','stringr','scales','hdm','glmnet',
        'survival')
sapply(ll,function(l) require(l,character.only = T))

##################################################
##### ----- EXERCISE SIMULATIONS: EX1 ----- ######
##################################################

# -------- CASE 0: NO COVARIATES EXCEPT AN INTERCEPT ----------- #

# base R defines the exponential as:
# S(ti|bi) = exp(- ti * bi); bi = b0 + lam*di

# Define a simple simulation function
ex1a.sim <- function(b0,lam,t.co) {
  n1 <- 100
  n2 <- 100
  # Create the treatment indicator
  d <- rbinom(n=n1+n2,size=1,prob=0.5)
  # Create the individual level means
  gami <- b0 + lam*d
  # Generate the survival distributions
  ti <- rexp(n=n1+n2,rate=gami)
  # Determine if observation is dead or alive over the time period of analysis
  # Encore 1 as survival
  y <- ifelse(ti>=t.co,1,0)
  # Put into a data.frame
  df <- tibble(y=y,d=d)
  # Run the logit
  logit <- glm(y~d,family = binomial(link='logit'))
  # Return the coefficients
  return(coef(logit))
}

# Define number of simulations and run it
nsim <- 1000
sim1a.df <- matrix(NA,nrow=nsim,ncol=2)
# So we should expect a mean of 0.5 and 1 for each group, respectively
b0 <- 2
lam <- -1
t.co <- 1/2
set.seed(1)
for (k in 1:nsim) { sim1a.df[k,] <- ex1a.sim(b0,lam,t.co) }
# Get the coefficient values
lam1a.mu <- apply(sim1a.df,2,mean)
b0a <- lam1a.mu[1]
lama <- lam1a.mu[2]

# We see the fitted probabilities are the same
1/(1+exp(-(b0a+lama)))
exp(-t.co*(b0+lam))

# Define the odds of survival for both
odds.d1 <- exp(-t.co*(b0+lam*1))/(1-exp(-t.co*(b0+lam*1)))
odds.d0 <- exp(-t.co*(b0+lam*0))/(1-exp(-t.co*(b0+lam*0)))

odds.d1/odds.d0
exp(lam1a.mu[2])

# -------- CASE 1: ONE COVARIATE WITH CONFOUNDING ----------- #

# base R defines the exponential as:
# S(ti|bi) = exp(- ti * bi); bi = b0 + psi*di

# Define the dgp with the confounders
dgp.1b <- function(g0,g1,psi,mu.d0,mu.d1,t.co,n0,n1) {
  # Create the treatment indicator
  d <- rep(c(0,1),times=c(n0,n1))
  # Create the covariate for each
  x.d1 <- rnorm(n=n1,mean=mu.d1,sd = 1)
  x.d0 <- rnorm(n=n0,mean=mu.d0,sd = 1)
  x <- c(x.d0,x.d1)
  # Create the individual level means
  gami <- g0 + g1*x + psi*d
  # Generate the survival distributions
  ti <- rexp(n=(n0+n1),rate=gami)
  # Determine if observation is dead or alive over the time period of analysis: encode 1 as survival
  y <- ifelse(ti>=t.co,1,0)
  # Put into a data.frame
  df <- tibble(y,d,x)
  # df %>% group_by(d) %>% summarise_each(funs(mean(.)))
  return(df)
}

# Define a simple simulation function
ex1b.sim <- function(g0,g1,psi,mu.d0,mu.d1,t.co,n0,n1) {
  # Run the dgp
  df <- dgp.1b(g0,g1,psi,mu.d0,mu.d1,t.co,n0,n1)
  # Run the logit
  logit <- glm(y~d+x,family = binomial(link='logit'),data=df)
  # Return the coefficients
  return(coef(logit))
}

# Define number of simulations and run it
nsim <- 5000
sim1b.df <- matrix(NA,nrow=nsim,ncol=3)
# Define the gamma
g0 <- 1; n0 <- 100
g1 <- 0.5; n1 <- 100
mu.d0 <- 4
mu.d1 <- 5
psi <- -1
t.co <- 1/3
set.seed(1)
for (k in 1:nsim) { sim1b.df[k,] <- ex1b.sim(g0,g1,psi,mu.d0,mu.d1,t.co,n0,n1) }
# Get the coefficient values
(lam1b.mu <- apply(sim1b.df,2,mean))
coefb <- t(data.frame(lam1b.mu)) %>% data.frame %>% set_colnames(c('b0','lam','b1'))

# Now run the simulation studies for the 50/50 case where everyone gets the treatment....
# Change in probability for the people who got the treatment
npeep <- 50000
# Generate the survival distributions
set.seed(1)
g1.treated <- rexp(n=npeep,rate=g0 + g1*rnorm(npeep,mu.d1) + psi*1) %>% na.omit
g1.notreat <- rexp(n=npeep,rate=g0 + g1*rnorm(npeep,mu.d1) + psi*0) %>% na.omit
g0.treated <- rexp(n=npeep,rate=g0 + g1*rnorm(npeep,mu.d0) + psi*1) %>% na.omit
g0.notreat <- rexp(n=npeep,rate=g0 + g1*rnorm(npeep,mu.d0) + psi*0) %>% na.omit
ga.treated <- rexp(n=npeep,rate=g0 + g1*rnorm(npeep,(mu.d0*n0+mu.d1*n1)/(n0+n1)) + psi*1) %>% na.omit
ga.notreat <- rexp(n=npeep,rate=g0 + g1*rnorm(npeep,(mu.d0*n0+mu.d1*n1)/(n0+n1)) + psi*0) %>% na.omit
# Store in list and then get the average effect
att.list <- list(g1.treated,g1.notreat,g0.treated,g0.notreat) %>%
              lapply(.,function(qq) ifelse(qq>t.co,1,0) %>% mean) %>% unlist
att.list2 <- list(ga.treated,ga.notreat) %>% lapply(.,function(qq) ifelse(qq>t.co,1,0) %>% mean) %>% unlist
# Odds increase of survival for treatment group
att.t <- (att.list[1]/(1-att.list[1]))/(att.list[2]/(1-att.list[2]))
att.nt <- (att.list[3]/(1-att.list[3]))/(att.list[4]/(1-att.list[4]))
# Next
lam.hata <- (att.t*n1 + att.nt*n0)/(n0+n1)
lam.hatb <- (att.list2[1]/(1-att.list2[1]))/(att.list2[2]/(1-att.list2[2]))
exp(lam1b.mu[2])

# Get the average value of x
mu.x <- (mu.d0+mu.d1)/2

# Now plot it
lam.b.df <- tibble(lam.hat=sim1b.df[,2])
# gg it
gg.simb <- ggplot(lam.b.df,aes(x=lam.hat)) + 
  geom_histogram(color='blue',fill='lightgrey',bins=30) + 
  geom_vline(xintercept = log(lam.hatb),color='red',linetype=2,size=1.5) + 
  labs(x=expression(hat(lambda)),y='Frequency',
       subtitle='Red line shows estimate of true effect',
       caption='5000 simulations')

# Do our own simulation
set.seed(1)
ni <- 100
# Create the treatment indicator
d <- rep(c(0,1),each=ni)
# Create the covariate for each
x.d1 <- rnorm(n=ni,mean=mu.d1,sd = 1)
x.d0 <- rnorm(n=ni,mean=mu.d0,sd = 1)
x <- c(x.d0,x.d1)
# Create the individual level means
gami <- g0 + g1*x + psi*d
# Generate the survival distributions
ti <- rexp(n=ni*2,rate=gami)
# Determine if observation is dead or alive over the time period of analysis: encode 1 as survival
y <- ifelse(ti>=t.co,1,0)
# Put into a data.frame
df <- tibble(y,d,x,ti)  
sim.Surv <- with(df,Surv(time=ti,event = rep(1,length(ti))))
sim.Survfit <- survfit(sim.Surv~d,data = df %>% mutate(d=ifelse(d==1,'Treatment','Non-treatment')))
library(survminer)
gg.survb <- ggsurvplot(sim.Survfit,legend.title = NULL,
           legend.labs = c('Non-treatment','Treatment'),
           legend='bottom') #c(0.5,0.75) 
