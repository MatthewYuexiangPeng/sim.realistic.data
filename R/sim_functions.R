# change path
setwd("D:/OneDrive - UW/Documents/GitHub/sim.realistic.data")

library(mvtnorm)

.ordgendata <- function(n, sigma, quants.norm){
  retval = mvtnorm::rmvnorm(n = n, sigma = sigma)
  for (i in 1:ncol(sigma)) {
    retval[, i] = cut(x = retval[, i], breaks = c(-1/0, quants.norm[[i]]), 
                      right = FALSE)
  }
  retval - 1  
}

.generate.cov <- function(n, P, Common.P, coef.AonB, coef.XonZ)
{
  B <- bindata::rmvbin(n, margprob=P, commonprob=Common.P)
  
  ### Categorical variables A
  ## Obtain correct predicted probabilities from multinomial coefficients by person
  P.A <- exp(cbind(rep(1,n),B)%*%t(coef.AonB[[1]]) )
  P.A <- cbind(rep(1,n),P.A)
  P.A <- P.A/apply(P.A,1,sum) # gives same result as fitted(multinom(age ~ z + x))
  
  # Generate Confounders A
  A.indicator <- t(apply(P.A,1,FUN=function(y){rmultinom(n=1,size=1,prob=y)})) 
  #A <- factor( as.character(c(1,2,3,4)%*%t(A.indicator)) )
  A.indicator <- A.indicator[,-1]
  
  # For multiple categorical variables
  n.A <- length(coef.AonB)
  if(n.A>1)
  {
    for(i in 2:length(coef.AonB))
    {
      P.A2 <- exp( cbind(rep(1,n),A.indicator,B)%*%t(coef.AonB[[i]]) )
      P.A2 <- cbind(rep(1,n),P.A2)
      P.A2 <- P.A2/apply(P.A2,1,sum) # gives same result as fitted(multinom(age ~ z + x))
      
      A.indicator.tmp <- t(apply(P.A2,1,FUN=function(y){rmultinom(n=1,size=1,prob=y)}))
      #A <- cbind(A, factor( as.character(c(1,2,3,4)%*%t(A.indicator.tmp)) ))
      A.indicator <- cbind(A.indicator, A.indicator.tmp[,-1])
    }
  }
  
  ### Confounders Z
  Z.model.data <- as.matrix(cbind(B, A.indicator))
  
  ### Exposure variable X
  X <- rbinom(n,size=1,
              p=1 / (1 + exp(-(coef.XonZ[1] + Z.model.data %*% coef.XonZ[2:length(coef.XonZ)]))))
  
  return(list(X=X, B=B, A.indicator=A.indicator))
}

.generate.cov.ord <- function(n, P.ord, Quant.norm, Corr.norm, coef.XonZ)
{
  
  Z <- .ordgendata(n, sigma=Corr.norm, quants.norm=Quant.norm)
  n.B <- length(which(unlist(lapply(P.ord,FUN=function(x){length(x)})==2)))
  n.A <- ncol(Z) - n.B
  
  if (n.B > 0 ) {
    B <- Z[, 1:n.B, drop=FALSE]
    colnames(B) <- paste0('bin', 1:n.B)
  }
  
  if(n.A>0)
  {
    A <- Z[,(n.B+1):ncol(Z)]
    A.indicator <- NULL
    catlabs <- NULL
    for(i in 1:n.A)
    {
      dummy <- NULL
      levels <- sort(unique(A[,i]))
      for(level in levels[-1]) {
        catlabs <- c(catlabs, paste0('cat',i,'_','level',level))
        dummy <- cbind(dummy, A[,i]==level)
      }
      A.indicator <- cbind(A.indicator, dummy)
    }
    colnames(A.indicator) <- catlabs
  }
  
  ### Confounders Z
  Z.model.data <- as.matrix(cbind(B, A.indicator))
  
  ### Exposure variable X
  X <- rbinom(n,size=1,
              p=1 / (1 + exp(-(coef.XonZ[1] + 
                                 Z.model.data %*% coef.XonZ[2:length(coef.XonZ)]))))
  
  return(list(X=X, B=B, A.indicator=A.indicator))
}

.generate.cov.chain <- function(n, coef.chain, coef.AonB, coef.XonZ,names.cat=NULL)
{
  n.B = length(coef.chain)
  B = rbinom(n=n,size=1,prob=coef.chain[[1]])
  for(i in 2:n.B){
    P.B = (1+exp(- (cbind(rep(1,n),B)%*%(coef.chain[[i]])) ))^(-1)
    B = cbind(B, rbinom(n=n,size=1,prob=P.B))
  }
  colnames(B) <- names(coef.chain)
  ### Categorical variables A
  ## Obtain correct predicted probabilities from multinomial coefficients by person
  P.A <- exp(cbind(rep(1,n),B)%*%t(coef.AonB[[1]]) )
  P.A <- cbind(rep(1,n),P.A)
  P.A <- P.A/apply(P.A,1,sum) # gives same result as fitted(multinom(age ~ z + x))
  
  # Generate Confounders A
  A.indicator <- t(apply(P.A,1,FUN=function(y){rmultinom(n=1,size=1,prob=y)})) 
  #A <- factor( as.character(c(1,2,3,4)%*%t(A.indicator)) )
  A.indicator <- A.indicator[,-1]
  
  # For multiple categorical variables
  n.A <- length(coef.AonB)
  if(n.A>1)
  {
    for(i in 2:length(coef.AonB))
    {
      P.A2 <- exp( cbind(rep(1,n),B)%*%t(coef.AonB[[i]]) ) # Matthew: exp( cbind(rep(1,n),A.indicator,B)%*%t(coef.AonB[[i]]) ) Question: want to regress A[2] on A[1]?
      P.A2 <- cbind(rep(1,n),P.A2)
      P.A2 <- P.A2/apply(P.A2,1,sum) # gives same result as fitted(multinom(age ~ z + x))
      
      A.indicator.tmp <- t(apply(P.A2,1,FUN=function(y){rmultinom(n=1,size=1,prob=y)}))
      #A <- cbind(A, factor( as.character(c(1,2,3,4)%*%t(A.indicator.tmp)) ))
      A.indicator <- cbind(A.indicator, A.indicator.tmp[,-1])
    }
  }
  
  colnames(A.indicator) <- names.cat
  
  ### Confounders Z
  Z.model.data <- as.matrix(cbind(B, A.indicator))
  
  ### Exposure variable X
  X <- rbinom(n,size=1,
              p=1 / (1 + exp(-(coef.XonZ[1] + Z.model.data %*% coef.XonZ[2:length(coef.XonZ)]))))
  
  
  return(list(X=X, B=B, A.indicator=A.indicator))
}

#' Title
#'
#' @param binary_outcome 
#' @param E 
#' @param Y 
#' @param X 
#' @param B 
#' @param A 
#' @param prescription.mode 
#' @param my.presc.K 
#' @param tie.method 
#' @param interact 
#'
#' @return
#' @export
#'
#' @examples
summary.stat <- function(binary_outcome=FALSE,E=NULL,Y,X,B,A,prescription.mode=seq(30,trunc,by=30),
                         my.presc.K=1,tie.method="efron",interact=FALSE)
{
  if(binary_outcome == FALSE){  ## Matthew added
    if (is.null(E)){
      stop("E is required for survival data")
    }
  }
  
  glm_control <- glm.control(epsilon = 1e-8, maxit = 25, trace = FALSE)
  cox.ctrl <- survival::coxph.control(eps=1e-09, toler.chol=.Machine$double.eps^0.75, 
                                      iter.max=20, toler.inf=sqrt(1e-09), outer.max=10)
  
  ### (Correlated) binary covariates B
  B <- as.matrix(B)
  n.B <- ncol(B)
  
  #### Marginal mean, i.e. prevalence of each binary covariate
  P <- apply(as.matrix(B),2,function(XX){as.numeric(table(XX)[[2]]/sum(table(XX)))})
  #### Common probability, i.e. number of one's shared by each pair of binary variables divided by sample size
  Common.P <- t(as.matrix(B))%*%as.matrix(B)/nrow(B)
  
  ### Categorical variables A
  #### Intercept and coefficients from multinomial regression for categorical variables, fitted on all binary covariates B. 
  #### Number of categorical
  A <- as.matrix(A)
  n.A <- ncol(A)
  
  ## Create list containing distributions of levels for 
  ## each categorical variable
  A.dist <- apply(A,2,function(XX){as.numeric(table(XX)/sum(table(XX)))})
  
  ### P.orddata and Corr.orddata
  P.ord = append(lapply(P,FUN=function(x){c(1-x,x)}),
                 A.dist)
  
  if (is.numeric(A)) {
    Corr.ord = cor(cbind(B,A))
  } else {
    Corr.ord = cor(cbind(B,as.factor(A)))
  }
  
  ### Categorical variables A
  #### Intercept and coefficients from multinomial regression for categorical variables, fitted on all binary covariates B. 
  #### Number of categorical
  #### A -> model.matrix type "A.indicator"
  if(n.A>0)
  {
    A.indicator <- NULL
    catlabs <- NULL
    for(i in 1:n.A)
    {
      dummy <- NULL
      levels <- sort(unique(A[,i]))
      for(level in levels[-1]) {
        catlabs <- c(catlabs, paste0('cat',i,'_','level',level))
        dummy <- cbind(dummy, A[,i]==level)
      }
      A.indicator <- cbind(A.indicator, dummy)
    }
    colnames(A.indicator) <- catlabs
  }
  
  
  ####################################  
  #### Chain of coef for B
  coef.chain = vector("list", n.B)
  names(coef.chain) <- colnames(B)
  coef.chain[[1]] = mean(B[,1])
  for(i in 2:n.B){
    coef.chain[[i]] = coef(glm(B[,i]~B[,1:(i-1)],family="binomial"))
  }
  
  
  #### Coef for relationship between binary and categorical confounders
  coef.AonB <- vector("list", n.A)
  names(coef.AonB) <- colnames(A)
  # Multinomial model for the first categorical variable and binary confounders
  coef.AonB[[1]] <- coef(nnet::multinom(A[,1] ~ B,trace=F))
  # Fits Multinomial Model if more than one categorical variable
  # Matthew: add the following
  if(n.A>1)
  {
    for(i in 2:n.A)
    {
      coef.AonB[[i]] <- coef(nnet::multinom(A[,i] ~ B,trace=F))
    }
  }
  
  ####################################
  
  ### Exposure variable X|B,A (Propensity Score Model Coefficients)
  #### Intercept and coefficients from logistic regression
  #dat <- data.frame(y=X,cbind(B,A))
  ps.fit <- glm.fit(cbind(1,B,A.indicator), X, control=glm_control, family=binomial())
  class(ps.fit) <- 'glm'
  coef.XonZ <- coef(ps.fit)
  ps.by.x <- cbind(X, fitted(ps.fit))
  ps.vcov <- vcov(ps.fit)
  probs <- predict(ps.fit, type = "response")
  ps.c <- mean(sample(probs[X == 1L], 1000000L, TRUE) > sample(probs[X == 0L], 1000000L, TRUE))
  
  if(binary_outcome == FALSE){
    ####  Censoring distribution -- Weibull shape and scale parameters estimated through parametric survival regression:
    #### Note: use only censoring data to estimate proportion among the censored.
    #### Unique censoring distribution shape
    #### Note: If one decide to add extra distribution properties, (1) and (2) should fit model using data eliminating the unique observations, 
    #### e.g. observations with E <- the top two most frequent presc pattern bumps
    
    ## Proportion with a prescription bump amongst those censored
    P.presc <- NULL
    for( i in 1:length(prescription.mode))
    {
      P.presc<-cbind(P.presc,mean(E[Y==0]==prescription.mode[i]))
    }
    P.presc <- data.frame(P.presc)
    names(P.presc) <- as.character(prescription.mode)
    
    #### CENSORING DISTRIBUTION #### 
    
    # Covariate matrix
    if (interact) {
      interact.X.covs <- cbind(B,A.indicator)*X
      colnames(interact.X.covs) <- paste("X",colnames(interact.X.covs),sep="_")
      Z <- cbind(X, cbind(B,A.indicator,interact.X.covs))
    } else {
      Z <- cbind(X, B, A.indicator)
    }
    
    
    # Simple Censoring No Bumps
    simple.cens <- survival::survreg(survival::Surv(E, 1-Y) ~ 1, dist="weib")
    simple.coef.cens <- coef(simple.cens)
    simple.scale.cens <- simple.cens$scale
    simple.vcov.cens <- vcov(simple.cens)
    names(simple.coef.cens) <- c("Intercept")
    
    # Covariate Censoring no Bumps
    rownames(Z) <- paste(1:dim(Z)[[1]])
    cov.cens <- survival::survreg(survival::Surv(E, 1-Y) ~ Z, dist="weib")
    cov.coef.cens <- coef(cov.cens)
    cov.scale.cens <- cov.cens$scale
    cov.vcov.cens <- vcov(cov.cens)
    names(cov.coef.cens) <- c("Intercept",colnames(Z))
    
    #Simple Censoring With Bumps
    # Select only the top my.presc.K bumps
    prescription.mode.topK <- prescription.mode[tail(order(unlist(P.presc)),my.presc.K)]
    P.presc.topK <- P.presc[tail(order(unlist(P.presc)),my.presc.K)]
    
    ## Bump Censoring Types
    ind<- (!E%in%prescription.mode.topK)
    
    simplebump.cens <- survival::survreg(survival::Surv(E[ind], 1-Y[ind]) ~ 1, dist="weib")
    simplebump.coef.cens <- coef(simplebump.cens)
    simplebump.scale.cens <- simplebump.cens$scale
    names(simplebump.coef.cens) <- c("Intercept")
    simplebump.vcov.cens <- vcov(simplebump.cens)
    
    # Covariate Censoring with Bumps
    covbump.cens <- survival::survreg(survival::Surv(E[ind], 1-Y[ind]) ~ Z[ind,], dist="weib")
    covbump.coef.cens <- coef(covbump.cens)
    covbump.scale.cens <- covbump.cens$scale
    names(covbump.coef.cens) <- c("Intercept",colnames(Z))
    covbump.vcov.cens <- vcov(covbump.cens)
    
    ### EVENT DISTRIBUTION
    #### Weibull shape and scale parameters estimated through parametric survival regression, 
    adj.event <- survival::survreg(survival::Surv(E, Y) ~ Z, dist="weib")
    
    if(sum(is.nan(coef(adj.event)))>0) {
      adj.event <- survival::survreg(survival::Surv(E, Y) ~ Z, dist="exp")
      adj.coef.event <- coef(adj.event)
      adj.scale.event <- 1
    } else {
      adj.coef.event <- coef(adj.event)
      adj.scale.event <- adj.event$scale
    }
    names(adj.coef.event) <- c("Intercept",colnames(Z))
    adj.vcov.event <- vcov(adj.event)
    
    #### Hazard ratio estimated through Cox PH model, w/adjustment for treatment X and covariates Z
    require(survival)
    cox.adjusted <- coxph(Surv(E, Y) ~ . , data=as.data.frame(model.matrix(~Z)[,-1]),
                          method=tie.method,x=TRUE)
    # class(cox.adjusted) <- "coxph"
    cox.coef.adjusted <- coef(cox.adjusted)
    names(cox.coef.adjusted) <- colnames(Z)
    cox.vcov <- vcov(cox.adjusted)
    
    ## Cox censoring model for plasmode simulation
    cox.adjusted.cens <-  coxph(Surv(E, 1-Y) ~ . , data=as.data.frame(model.matrix(~Z)[,-1]),
                                method=tie.method,x=TRUE)
    # class(cox.adjusted.cens) <- "coxph"
    
    #### Additional information
    control.rate <- sum(Y[X==0])/sum(E[X==0])
    compare.rate <- sum(Y[X==1])/sum(E[X==1])
    control.events <- sum(Y[X==0])
    compare.events <- sum(Y[X==1])
    N.X <- table(X)
    P.time <- c(sum(E[X==0]), sum(E[X==1]))/table(X)
  } else {
    #### Matthew added: Binary outcome
    #### Logistic regression for the event distribution
    
    ps.fit <- glm.fit(cbind(1,X,B,A.indicator), Y, control=glm_control, family=binomial())
    class(ps.fit) <- 'glm'
    coef.Yon1 <- coef(ps.fit)[1]
    coef.YonX <- coef(ps.fit)[2]
    coef.YonZ <- coef(ps.fit)[-c(1,2)]
    
    # matthew edited here 
    # dat2 <- data.frame(y=Y,cbind(X,B,A.indicator))
    # coef.Yon1 <- coef(glm(y~.,data=dat2,family="binomial"))[1]
    # coef.YonX <- coef(glm(y~.,data=dat2,family="binomial"))[2]
    # coef.YonZ <- coef(glm(y~.,data=dat2,family="binomial"))[-c(1,2)]
    
    control.events <- sum(Y[X==0])
    compare.events <- sum(Y[X==1])
    N.X <- table(X)
  }  
  
  # retrun different objects
  if(binary_outcome == FALSE){  ## Matthew added
    return(list(n=length(Y),N.X=N.X, P.time=P.time,
                control.events=control.events,compare.events=compare.events,
                control.rate=control.rate,compare.rate=compare.rate,
                P=P,Common.P=Common.P,P.ord=P.ord,
                Corr.ord=Corr.ord,coef.XonZ=coef.XonZ,Coef.bin=coef.chain, Coef.cat=coef.AonB,
                propensity=ps.by.x, 
                propensity.vcov=ps.vcov,P.presc=P.presc, cStat=ps.c,
                prescription.mode=prescription.mode,P.presc.topK=P.presc.topK,
                prescription.mode.topK=prescription.mode.topK,
                simple.coef.cens=simple.coef.cens,
                simple.scale.cens=simple.scale.cens, simple.vcov.cens=simple.vcov.cens,
                simplebump.coef.cens=simplebump.coef.cens,
                simplebump.scale.cens=simplebump.scale.cens,
                simplebump.vcov.cens=simplebump.vcov.cens,
                cov.coef.cens=cov.coef.cens,cov.scale.cens=cov.scale.cens,
                cov.vcov.cens=cov.vcov.cens,
                covbump.coef.cens=covbump.coef.cens,
                covbump.scale.cens=covbump.scale.cens,
                covbump.vcov.cens=covbump.vcov.cens,
                adj.coef.event=adj.coef.event,adj.scale.event=adj.scale.event,
                adj.vcov.event=adj.vcov.event,cox.coef.adjusted=cox.coef.adjusted,cox.vcov=cox.vcov,
                cox.fit.event=cox.adjusted,cox.fit.cens=cox.adjusted.cens
    ))
  } else {
    return(list(n=length(Y),N.X=N.X,
                control.events=control.events,compare.events=compare.events,
                P=P,Common.P=Common.P,P.ord=P.ord,
                Corr.ord=Corr.ord,coef.XonZ=coef.XonZ,Coef.bin=coef.chain, Coef.cat=coef.AonB,
                propensity=ps.by.x, 
                propensity.vcov=ps.vcov,cStat=ps.c,
                coef.Yon1=coef.Yon1,coef.YonX=coef.YonX,coef.YonZ=coef.YonZ  ## Matthew added
    ))
  }
}
