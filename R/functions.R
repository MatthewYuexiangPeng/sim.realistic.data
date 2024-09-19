
library(mvtnorm)

#### Generate data

generate.cov <- function(n, P, Common.P, coef.AonB, coef.XonZ)
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

generate.cov.ord <- function(n, P.ord, Quant.norm, Corr.norm, coef.XonZ)
{

  Z <- ordgendata(n, sigma=Corr.norm, quants.norm=Quant.norm)
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

generate.cov.chain <- function(n, coef.chain, coef.AonB, coef.XonZ,names.cat=NULL)
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

## OrdToNorm (appropriated from orddata)
## Description: Computes multivariate normal correlation
## and quantiles from categorical correlation and marginal
## probabilities.
##
## Inputs:
##    probs = list of marginal probabilites for categorical
##            and binary variables.
##    cor   = Pearson correlation matrix compute on observed
##            categorical and binary covariates.
##
## Outputs:
##    Type  = List
##    Values
##      Corr.norm = Multivariate normal correlation matrix.
##      quants.norm = Quantiles for thresholding marginal normals.
##

OrdToNorm <- function (probs, Cor) 
{
  q = length(probs)
  categ_probs = 0
  cumul_probs = list(0)
  quant_probs = list(0)
  means = 0
  vars = 0
  var.wt = function(x, w) {
    m = weighted.mean(x = x, w = w)
    sum((x[1:length(x)] - m)^2 * w[1:length(x)])
  }
  for (i in 1:q) {
    categ_probs[i] = length(probs[[i]])
    cumul_probs[[i]] = cumsum(1:categ_probs[i]/10^12 + probs[[i]])
    cumul_probs[[i]][categ_probs[i]] = 1
    quant_probs[[i]] = qnorm(p = cumul_probs[[i]], mean = 0, 
                             sd = 1)
    means[i] = weighted.mean(x = 1:categ_probs[i], w = probs[[i]])
    vars[i] = var.wt(x = 1:categ_probs[i], w = probs[[i]])
  }
  Cor_norm = Cor
  for (i in 1:(q - 1)) {
    for (j in (i + 1):q) {
      gridd = rep(0, times = 201)
      for (ci in 1:(categ_probs[i] - 1)) {
        for (cj in 1:(categ_probs[j] - 1)) {
          for (steps in -100:100) {
            gridd[101 + steps] = gridd[101 + steps] + 
              mvtnorm::pmvnorm(upper = c(quant_probs[[i]][ci], 
                                quant_probs[[j]][cj]), corr = matrix(2, 
                                                                     2, data = c(1, steps/100, steps/100, 
                                                                                 1)))[1]
          }
        }
      }
      f = suppressWarnings(approxfun(y = -100:100/100, 
                                     x = gridd))
      Cor_norm[i, j] = Cor_norm[j, i] = f(Cor[i, j] * sqrt(vars[i] * 
                                                             vars[j]) + means[i] * means[j] - categ_probs[i] * 
                                            categ_probs[j] + categ_probs[j] * sum(cumul_probs[[i]][1:(categ_probs[i] - 
                                                                                                        1)]) + categ_probs[i] * sum(cumul_probs[[j]][1:(categ_probs[j] - 
                                                                                                                                                          1)]))
    }
  }
  return(list('corr.norm'=Cor_norm,'quants.norm'=quant_probs))
  
}



## ordgendata (appropriated from orddata)
## Description: Simulates matrix of correlated categorical and binary covariates.
##
## Inputs:
##    n       = Sample size
##            
##    sigma   = Correlation matrix for multivariate normal.
##  
##    quants.norm = Quatile thresholds for binning marginal normals.
##
## Outputs:
##    Type  = List
##    Values
##      
##      
##
ordgendata <- function(n, sigma, quants.norm){
  retval = mvtnorm::rmvnorm(n = n, sigma = sigma)
  for (i in 1:ncol(sigma)) {
    retval[, i] = cut(x = retval[, i], breaks = c(-1/0, quants.norm[[i]]), 
                      right = FALSE)
  }
  retval - 1  
}


###################################################################################
#### summary.stat
#### DESCRIPTION: Obtains the site summary information needed for data generation 
#### INPUTS:
#### E = observed time (in days)
#### Y = event indicator (1 or 0)
#### X = exposure (1 or 0)
#### B = Binary confounders
#### A = Categorical confounders (categories of variables not indicators variables)
#### prescription.mode = possible different 30 day prescription patterns that we may observe
####		(e.g. prescription.mode = seq(30,180, by=30))
#### my.presc.K = number of prescription bumps to include in censoring dist (e.g. my.presc.K = 2)
#### OUTPUTS:
#### 
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


###################################################################################
#### generate.data
#### DESCRIPTION:  
#### INPUTS:
#### censtype = "simple", "simplebump", "cov", "covbump"
####
#### OUTPUTS:
#### 
generate.data <- function(binary_outcome=FALSE, logHR.X.site=NULL, n, P, Common.P, coef.AonB=NULL, coef.XonZ,
                          coef.cens=NULL, scale.cens=NULL,
                          coef.event=NULL, scale.event=NULL,
                          censtype="simple", trunc=366,
                          P.presc.topK=NULL, prescription.mode.topK=NULL,
                          method=1, Corr.norm=NULL, Quant.norm=NULL, P.ord=NULL,
                          coef.chain=NULL, user.data=NULL,noX=FALSE,
                          strat.var.cens=NULL,strat.var.event=NULL,
                          coef.Yon1, coef.YonX, coef.YonZ)  ## Matthew added
{
  # # test 9.10
  # binary_outcome=FALSE
  # logHR.X.site=NULL
  # coef.AonB=NULL
  # coef.cens=NULL
  # scale.cens=NULL
  # coef.event=NULL
  # scale.event=NULL
  # censtype="simple"
  # trunc=366
  # P.presc.topK=NULL
  # prescription.mode.topK=NULL
  # method=3
  # Corr.norm=NULL
  # Quant.norm=NULL
  # P.ord=NULL
  # coef.chain=NULL
  # user.data=NULL
  # noX=FALSE
  # strat.var.cens=NULL
  # strat.var.event=NULL
  # 
  # 
  # binary_outcome = TRUE
  # logHR.X.site=NULL
  # n=SS$n
  # P=SS$P
  # Common.P=SS$Common.P
  # coef.XonZ=SS$coef.XonZ
  # coef.chain=SS$Coef.bin
  # coef.AonB=SS$Coef.cat
  # method=method
  # Corr.norm=SS$Corr.norm
  # Quant.norm=SS$Quants.norm
  # P.ord=SS$P.ord
  # coef.Yon1=SS$coef.Yon1
  # coef.YonX=SS$coef.YonX
  # coef.YonZ=SS$coef.YonZ
  # 
  # ####
  
  ### set up specified HR of X
  if(binary_outcome == FALSE){  ## Matthew added
    if (!noX) {
      coef.event[grep("X",names(coef.event))] <- -logHR.X.site*scale.event
    }
    if (is.null(logHR.X.site|coef.cens|scale.cens|coef.event|scale.event)){
      stop("For binary outcome, logHR.X.site, coef.cens, scale.cens, coef.event, and scale.event must be specified")
    }
    ### set up specified HR of X
    coef.event[grep("X",names(coef.event))] <- -logHR.X.site*scale.event
  } else {
    if (is.null(coef.Yon1|coef.YonX|coef.YonZ)){
      stop("For binary outcome, coef.Yon1, coef.YonX, and coef.YonZ must be specified")
    }
  }
  
  ### Generate covariates via specified method

  
  if (method == 1) {
    
    exppluscovs <- generate.cov.ord(n, P.ord, Quant.norm, Corr.norm, coef.XonZ)
    
  } else if (method==2) {
    
    exppluscovs <- generate.cov(n, P, Common.P, coef.AonB, coef.XonZ)
    
  } else if (method==3) {
    
    exppluscovs <- generate.cov.chain(n, coef.chain, coef.AonB, coef.XonZ)
    
  } else if (method==4) {
    exppluscovs <- list(B=NULL,A.indicator=NULL)
    exppluscovs$A.indicator <- user.data[,-c(1)]
    if(noX) {
      exppluscovs$X <- NULL
    } else {
      exppluscovs$X <- user.data[,1]
    }
    
  }
  
  ### Confounders Z
  Z.model.data <- cbind(exppluscovs$B, exppluscovs$A.indicator)
  
  ### Exposure variable X
  X <- exppluscovs$X
  
  if(binary_outcome == FALSE){
    ### Censoring
    censorT<-NULL
    u <- runif(n)
  
    if(censtype%in%c("simple","simplebump"))
    {
      censorT <- ceiling((-log(u) * exp(coef.cens / scale.cens)) ^ (scale.cens))
      censorT.1 <- ceiling((-log(runif(n)) * exp(coef.cens / scale.cens)) ^ (scale.cens))
      censorT.0 <- ceiling((-log(runif(n)) * exp(coef.cens / scale.cens)) ^ (scale.cens))
    }
    if(censtype%in%c("cov","covbump"))
    {
      if (length(scale.cens)==1) {
        
        censorT <- ceiling((-log(u)*exp( cbind(1,X,Z.model.data) %*% coef.cens/scale.cens ))^(scale.cens))
        if (!is.null(X)) {
          censorT.1 <- ceiling((-log(runif(n))*exp( cbind(rep(1,n),1,Z.model.data) %*% coef.cens/scale.cens ))^(scale.cens))
          censorT.0 <- ceiling((-log(runif(n))*exp( cbind(rep(1,n),0,Z.model.data) %*% coef.cens/scale.cens ))^(scale.cens))
        } else {
          censorT.1 <- ceiling((-log(runif(n))*exp( cbind(rep(1,n),Z.model.data) %*% coef.cens/scale.cens ))^(scale.cens))
          censorT.0 <- ceiling((-log(runif(n))*exp( cbind(rep(1,n),Z.model.data) %*% coef.cens/scale.cens ))^(scale.cens))
        }
        
      } else {
        
        censorT <- rep(0,times=n)
        censorT[strat.var.cens==0] <- ceiling((-log(u[strat.var.cens==0])*exp( 
          cbind(1,X,Z.model.data)[strat.var.cens==0,] %*% 
            coef.cens/scale.cens[[1]] ))^(scale.cens[[1]]))
        censorT[strat.var.cens==1] <- ceiling((-log(u[strat.var.cens==1])*exp( 
          cbind(1,X,Z.model.data)[strat.var.cens==1,] %*% 
            coef.cens/scale.cens[[2]] ))^(scale.cens[[2]]))
      }
      
    }
    if(is.null(censorT)) stop("censoring model is misspecified")
    
    ### Event time
    ue <- runif(n)
    
    if (length(scale.event)==1) {
      eventT <- ceiling((-log(ue) * exp( cbind(1,X,Z.model.data) %*% coef.event / scale.event )) ^ (scale.event))
      if (!is.null(X)) {
        eventT.1 <- ceiling((-log(runif(n)) * exp( cbind(1,1,Z.model.data) %*% coef.event / scale.event )) ^ (scale.event))
        eventT.0 <- ceiling((-log(runif(n)) * exp( cbind(1,0,Z.model.data) %*% coef.event / scale.event )) ^ (scale.event))
      } else {
        eventT.1 <- ceiling((-log(runif(n)) * exp( cbind(1,Z.model.data) %*% coef.event / scale.event )) ^ (scale.event))
        eventT.0 <- ceiling((-log(runif(n)) * exp( cbind(1,Z.model.data) %*% coef.event / scale.event )) ^ (scale.event))
      }
      
    } else {
      eventT <- rep(0,times=n)
      eventT[strat.var.event==0] <- ceiling((-log(ue[strat.var.event==0])*exp( 
        cbind(1,X,Z.model.data)[strat.var.event==0,] %*% 
          coef.event/scale.event[[1]] ))^(scale.event[[1]]))
      eventT[strat.var.event==1] <- ceiling((-log(ue[strat.var.event==1])*exp( 
        cbind(1,X,Z.model.data)[strat.var.event==1,] %*% 
          coef.event/scale.event[[2]] ))^(scale.event[[2]]))
    }
    
    
    ### Observed time and event indicator
    if(censtype%in%c("simplebump","covbump"))
    {
  
      trunc.presc <- apply(rmultinom(n,1,c(P.presc.topK,'NA' = 1-sum(P.presc.topK))),
                           2,function(XX) {c(prescription.mode.topK,NA)[which(XX==1)]})
      CensorTbump <- censorT
      CensorTbump[!is.na(trunc.presc)] <- trunc.presc[!is.na(trunc.presc)]
      
      E <- apply(cbind(eventT,trunc,CensorTbump),1,min,na.rm=T)
      Y   <- as.numeric(ifelse(E == eventT,1,0))
      
      CensorTbump.1 <- censorT.1
      CensorTbump.1[!is.na(trunc.presc)] <- trunc.presc[!is.na(trunc.presc)]
      CensorTbump.0 <- censorT.0
      CensorTbump.0[!is.na(trunc.presc)] <- trunc.presc[!is.na(trunc.presc)]
      
      E.1 <- apply(cbind(eventT.1,trunc,CensorTbump.1),1,min,na.rm=T)
      E.0 <- apply(cbind(eventT.0,trunc,CensorTbump.0),1,min,na.rm=T)
      
      Y.1 <- as.numeric(ifelse(E.1 == eventT.1,1,0))
      Y.0 <- as.numeric(ifelse(E.0 == eventT.0,1,0))
      
    }
    else
    {
      E <- apply(cbind(eventT,trunc,censorT),1,min,na.rm=T)
      Y   <- as.numeric(ifelse(E == eventT,1,0))
      
      E.1 <- apply(cbind(eventT.1,trunc,censorT.1),1,min,na.rm=T)
      E.0 <- apply(cbind(eventT.0,trunc,censorT.0),1,min,na.rm=T)
      
      Y.1 <- as.numeric(ifelse(E.1 == eventT.1,1,0))
      Y.0 <- as.numeric(ifelse(E.0 == eventT.0,1,0))
    }
  } else {
    #### Matthew added: Binary outcome
    Y <- rbinom(n,size=1,
                p=1 / (1 + exp(-(coef.Yon1 + X * coef.YonX + Z.model.data %*% coef.YonZ))))
    Y.1 <- rbinom(n,size=1,
                  p=1 / (1 + exp(-(coef.Yon1 + coef.YonX + Z.model.data %*% coef.YonZ))))
    Y.0 <- rbinom(n,size=1,
                  p=1 / (1 + exp(-(coef.Yon1 + Z.model.data %*% coef.YonZ))))
  }
  
  ## Marginal sample needs to be same size as original
  marg.x1 <- sample(1:n,size=sum(X))
  
  if(binary_outcome == FALSE){  ## Matthew added
    return(list(marginal.dat=rbind(cbind(x=1, y=Y.1, obst=E.1)[marg.x1,],cbind(x=0, y=Y.0, obst=E.0)[-marg.x1,]),
                B=exppluscovs$B,A.indicator=exppluscovs$A.indicator, 
           X=exppluscovs$X,Y=Y,E=E)
    )
  } else {
    return(list(marginal.dat=rbind(cbind(x=1, y=Y.1)[marg.x1,],cbind(x=0, y=Y.0)[-marg.x1,]),
                B=exppluscovs$B,A.indicator=exppluscovs$A.indicator, 
           X=exppluscovs$X,Y=Y)
    )
  }
}

###################################################################################
#### generate.data.full
#### DESCRIPTION:  
#### INPUTS:
#### 
####
#### OUTPUTS:
#### 
generate.data.full <- function(binary_outcome = FALSE,Summ.Stat,hetero=0,censtype="simple", trunc=365,method=1)
{
  # # test 9.10
  # method = 3
  # binary_outcome = TRUE
  # censtype=censtype
  # trunc=366
  # ####
  
  glm_control <- glm.control(epsilon = 1e-8, maxit = 25, trace = FALSE)
  n.sites<-length(Summ.Stat)
  
  ## Data.Site: simulated data across sites
  Data.Simulated <- NULL
  Data.Marg <- NULL
  
  ## Data.Pool.forPS: pooled data for est pooled PS
  Data.Pool.forPS <- NULL
  
  ## Loop through sites to generate site specific data
  for(i in 1:n.sites)
  {
    SS <- NULL
    ## read in summary statistics SS
    SS <- Summ.Stat[[i]]
    ## generate site specific data
    DS<-NULL
    
    # method 2, 3, 4 not working
    if(binary_outcome == FALSE){
      if(censtype=="simple")
      {
        DS <- generate.data(logHR.X.site=SS$logHR.X,
                            n=SS$n, P=SS$P, Common.P=SS$Common.P, coef.XonZ=SS$coef.XonZ,
                            coef.cens=SS$simple.coef.cens, scale.cens=SS$simple.scale.cens,
                            coef.event=c(SS$intercept,SS$adj.coef.event[-1]), scale.event=SS$adj.scale.event,
                            censtype=censtype,trunc=trunc,coef.chain=SS$Coef.bin, coef.AonB=SS$Coef.cat,
                            P.presc.topK=NULL, prescription.mode.topK=NULL,
                            method=method, Corr.norm=SS$Corr.norm, Quant.norm=SS$Quants.norm, P.ord=SS$P.ord)	
      }
      
      if(censtype=="simplebump")
      {
        DS <- generate.data(logHR.X.site=SS$logHR.X,
                            n=SS$n, P=SS$P, Common.P=SS$Common.P, coef.XonZ=SS$coef.XonZ,
                            coef.cens=SS$simplebump.coef.cens, scale.cens=SS$simplebump.scale.cens,
                            coef.event=c(SS$intercept,SS$adj.coef.event[-1]), scale.event=SS$adj.scale.event,
                            censtype=censtype,trunc=trunc,coef.chain=SS$Coef.bin, coef.AonB=SS$Coef.cat,
                            P.presc.topK=SS$P.presc.topK, 
                            prescription.mode.topK=SS$prescription.mode.topK,
                            method=method, Corr.norm=SS$Corr.norm, Quant.norm=SS$Quants.norm, P.ord=SS$P.ord)
        
      }
      
      if(censtype=="cov")
      {
        DS <- generate.data(logHR.X.site=SS$logHR.X,
                            n=SS$n, P=SS$P, Common.P=SS$Common.P, coef.XonZ=SS$coef.XonZ,
                            coef.cens=SS$cov.coef.cens, scale.cens=SS$cov.scale.cens,
                            coef.event=c(SS$intercept,SS$adj.coef.event[-1]), scale.event=SS$adj.scale.event,
                            censtype=censtype,trunc=trunc,coef.chain=SS$Coef.bin, coef.AonB=SS$Coef.cat,
                            P.presc.topK=NULL, prescription.mode.topK=NULL,
                            method=method, Corr.norm=SS$Corr.norm, Quant.norm=SS$Quants.norm, P.ord=SS$P.ord)	
      }
      
      if(censtype=="covbump")
      {
        DS <- generate.data(logHR.X.site=SS$logHR.X,
                            n=SS$n, P=SS$P, Common.P=SS$Common.P, coef.XonZ=SS$coef.XonZ,
                            coef.cens=SS$covbump.coef.cens, scale.cens=SS$covbump.scale.cens,
                            coef.event=c(SS$intercept,SS$adj.coef.event[-1]), scale.event=SS$adj.scale.event,
                            censtype=censtype,trunc=trunc,coef.chain=SS$Coef.bin, coef.AonB=SS$Coef.cat,
                            P.presc.topK=SS$P.presc.topK, prescription.mode.topK=SS$prescription.mode.topK,
                            method=method, Corr.norm=SS$Corr.norm, Quant.norm=SS$Quants.norm, P.ord=SS$P.ord)	
      }
      if(is.null(DS)) stop("Censoring type is misspecified")
    } else { 
      DS <- generate.data(binary_outcome = TRUE, logHR.X.site=NULL, n=SS$n, P=SS$P, Common.P=SS$Common.P, coef.XonZ=SS$coef.XonZ,
                          coef.chain=SS$Coef.bin, coef.AonB=SS$Coef.cat,
                          method=method, Corr.norm=SS$Corr.norm, Quant.norm=SS$Quants.norm, P.ord=SS$P.ord,
                          coef.Yon1=SS$coef.Yon1, coef.YonX=SS$coef.YonX, coef.YonZ=SS$coef.YonZ) # Matthew added
    }
    
    ## save site specific data
    if (binary_outcome == FALSE) {
      if (n.sites > 1) {
  
        Data.Simulated <- rbind(Data.Simulated, data.frame(B=DS$B,
                                                      A=DS$A.indicator, 
                                                      X=DS$X,Y=DS$Y,E=DS$E,
                                                      site=i))
        Data.Marg <- rbind(Data.Marg,data.frame(site=i,DS$marginal.dat))
        
      } else {
        Data.Simulated <- data.frame(B=DS$B, A=DS$A.indicator, X=DS$X,Y=DS$Y,E=DS$E)
        Data.Marg <- data.frame(DS$marginal.dat)
      }
    } else {
      if (n.sites > 1) {
        Data.Simulated <- rbind(Data.Simulated, data.frame(B=DS$B, 
                                                           A=DS$A.indicator, 
                                                           X=DS$X,Y=DS$Y,
                                                           site=i))
        Data.Marg <- rbind(Data.Marg,data.frame(site=i,DS$marginal.dat))
      } else {
        Data.Simulated <- data.frame(B=DS$B, A=DS$A.indicator, X=DS$X,Y=DS$Y)
        Data.Marg <- data.frame(DS$marginal.dat)
      }
    }
    
  }

  return(list(Data.Simulated=Data.Simulated,Data.Marginal=Data.Marg))
}


###################################################################################
#### method
#### DESCRIPTION:  
#### INPUTS:
#### 
####
#### OUTPUTS:
#### 
method <- function(data.full, cov.cols, cox.ctrl=cox.ctrl,
                   tie.method=tie.method,truth.logHR.X,
                   n.PSstr=5,bs.df=5,bs.degree=3,test.side=1,alpha=0.05)
{
  data.marg <- data.full$Data.Marginal
  data.full <- data.full$Data.Simulated
  X <- data.full$X
  Y <- data.full$Y
  E <- data.full$E
  S <- data.full$site
  Z <- cbind(as.matrix(data.full[,cov.cols]), model.matrix(~factor(S))[,-1]) 
  Z.nosite <- as.matrix(data.full[,cov.cols])

  ## Two interaction models for marginal hazard

  ## Site interacts will all Z only main effect of X
  marg.nonstrat <- survival::coxph.fit(model.matrix(~ X + Z.nosite*factor(S))[,-c(1)],
                                   survival::Surv(E, Y),
                                   strata=NULL, control=cox.ctrl,
                                   method=tie.method, rownames=NULL)
  BetaX.marg.nonstrat <- coef(marg.nonstrat)[[1]]
  cover.marg.nonstrat <- ci.coverage(marg.nonstrat, truth=truth.logHR.X)

  ## Site interacts will all Z but no main effect of site. Instead
  ## site is used as a stratification variable. Only main effect of X
  marg.strat <- survival::coxph.fit(model.matrix(~X + Z.nosite + Z.nosite:factor(S))[,-c(1)],
                                       survival::Surv(E, Y),
                                       strata=S, control=cox.ctrl,
                                       method=tie.method, rownames=NULL)
  BetaX.marg.strat <- coef(marg.strat)[[1]]
  cover.marg.strat <- ci.coverage(marg.strat, truth=truth.logHR.X)
  
  ## Simulated marginal cox
  marg <- survival::coxph.fit(model.matrix(~ data.marg$x + factor(data.marg$site))[,-c(1)],
                                    survival::Surv(data.marg$obst, data.marg$y),
                                    strata=NULL, control=cox.ctrl,
                                    method=tie.method, rownames=NULL)
  BetaX.marg <- coef(marg)[[1]]
  cover.marg <- ci.coverage(marg, truth=truth.logHR.X)

  #### Mantel-Haenzel Estimator
  mh.fit          <- coxmh(E=E,Y=Y,X=X,Z=Z[,cov.cols],S=S,test.side=2)
  BetaX.mh.samp   <- mh.fit$SampWeight[[1]][[1]]
  BetaX.mh.var    <- mh.fit$VarWeight[[1]][[1]]
  cover.mh.samp   <- (truth.logHR.X >= mh.fit$SampWeight[[3]][[1]]) & 
    (truth.logHR.X <= mh.fit$SampWeight[[3]][[2]])
  cover.mh.var    <- (truth.logHR.X >= mh.fit$VarWeight[[3]][[1]]) & 
    (truth.logHR.X <= mh.fit$VarWeight[[3]][[2]])
  signorm.mh.samp <- mh.fit$SampWeight[[5]][[1]]
  signorm.mh.var  <- mh.fit$VarWeight[[5]][[1]]
  Score.mh.samp   <- mh.fit$SampWeight[[4]][[1]]
  Score.mh.var    <- mh.fit$VarWeight[[4]][[1]]
 
  #### Compute pooled and site-specific propensity scores
  PS.list <- comp.pscore(X=data.full$X,Z=data.full[,cov.cols],S=data.full$site)
  
  #### Set pooled p-scores
  PS.Pool <- PS.list$PS.pooled
  #### Set site-specific p-scores
  PS.Site.list <- PS.list$PS.site
  #### Count the number of sites
  n.sites <- length(PS.Site.list)

  #### Site-specific p-scores, p-score categories & b-splines
  PS.Site <- PSstr.Site <- PSspline.Site <- NULL ;

  for(i in 1:n.sites)
  {
    PS.Site <- c(PS.Site,PS.Site.list[[i]])
    PSstr.Site <- c(PSstr.Site, cut(PS.Site.list[[i]],include.lowest=T,
                                    breaks=quantile(PS.Site.list[[i]],
                                    seq(0,1,1/n.PSstr),names=F),
                                    labels=FALSE))
    PSspline.Site <- rbind(PSspline.Site, splines::bs(PS.Site.list[[i]],
                                                      df=bs.df, degree=bs.degree))
  }
  
  ## Pooled b-splines
  PSspline <- splines::bs(PS.Pool, df=bs.df, degree=bs.degree)
  ## Pooled p-score categories
  PSstr <- cut(PS.Pool, include.lowest=T, 
               breaks=quantile(PS.Pool, seq(0,1,1/n.PSstr), names=F))
  ## Site-specific p-score strata
  PSstr.by.site <- as.character(PSstr.Site+S*100000)
  
  Z.pstr <- model.matrix(~ PSstr)[,c(-1)]
  Z.bspline <- PSspline
  
  ## site-specific PS.spline/strata
  Z.pstr.site <- model.matrix(~ factor(PSstr.Site)*factor(S))[,-1]
  Z.bspline.site <- model.matrix(~ PSspline.Site*factor(S))[,-1]
  Z.bspline.site.str <- Z.bspline.site[,-c((bs.df+1):(bs.df + n.sites - 1))]
  
  #### Mantel-Haenzel Estimator with b-splines
  mh.bs.fit <- coxmh(E=E,Y=Y,X=X,Z=PSspline.Site,S=S,test.side,alpha)
  BetaX.mh.bs.samp <- mh.bs.fit$SampWeight[[1]][[1]]
  BetaX.mh.bs.var <-mh.bs.fit$VarWeight[[1]][[1]]
  cover.mh.bs.samp <- truth.logHR.X >= mh.bs.fit$SampWeight[[3]][[1]] & 
    truth.logHR.X <= mh.bs.fit$SampWeight[[3]][[2]]
  cover.mh.bs.var <- truth.logHR.X >= mh.bs.fit$VarWeight[[3]][[1]] & 
    truth.logHR.X <= mh.bs.fit$VarWeight[[3]][[2]]
  signorm.mh.bs.samp <- mh.bs.fit$SampWeight[[5]][[1]]
  signorm.mh.bs.var <- mh.bs.fit$VarWeight[[5]][[1]]
  Score.mh.bs.samp <- mh.bs.fit$SampWeight[[4]][[1]]
  Score.mh.bs.var <- mh.bs.fit$VarWeight[[4]][[1]]
  
  ###################################
  #### Cox PH Regression Methods #### 
  
  ## [Unadj] Cox PH regression without adjusting for confounders ##
  unadj.fit <- survival::coxph.fit(as.matrix(as.numeric(X)), survival::Surv(E, Y), 
                                   strata=NULL, control=cox.ctrl, method=tie.method, rownames=NULL)
  BetaX.unadj <- coef(unadj.fit)
  cover.unadj <- ci.coverage(unadj.fit, truth=truth.logHR.X)

  
  #######################################################
  #### POOLED METHODS (E.G. NON-DISTRIBUTED SETTING) ####
  
  #### [Adj Confounders+Site] Regression on categorical confounders and site

  ## NULL Model
  BetaZ.adjcon <- as.numeric(survival::coxph.fit(as.matrix(Z), survival::Surv(E, Y),
                                          strata=NULL,control=cox.ctrl,
                                          method=tie.method,
                                          rownames=NULL)[[c('coefficients')]])
  
  ## Score statistic ##
  if (sum(!is.finite(BetaZ.adjcon)) == 0) {
    ScT.adjcon<-ScTt.cox(delta=Y, time=E, X=X, Z=Z, beta=BetaZ.adjcon)
    if (is.finite(ScT.adjcon)) {
      if (test.side==1) {
        signorm.adjcon <- ScT.adjcon >= qnorm(1-alpha)
      } else if (test.side==2) {
        signorm.adjcon <- ScT.adjcon <= qnorm(alpha)
      }
    } else {
      ScT.adjcon <- NA
      signorm.adjcon <- NA
    }
    
  } else {
    ScT.adjcon <- NA
    signorm.adjcon <- NA
  }

  ## Estimation model
  adjcon.fit <- survival::coxph.fit(cbind(X,Z), survival::Surv(E, Y),
                                    strata=NULL, control=cox.ctrl, method=tie.method,rownames=NULL)
  BetaX.adjcon <- coef(adjcon.fit)[[1]]
  cover.adjcon <- ci.coverage(adjcon.fit, truth=truth.logHR.X)
  
  
  #### [Adj PS Indicators] Regression on Propensity Score Indicators ####

  ## NULL MOdel
  BetaZ.adjpsind <- as.numeric(survival::coxph.fit(Z.pstr, survival::Surv(E, Y),
                                            strata=NULL, 
                                            control=cox.ctrl,
                                            method=tie.method,
                                            rownames=NULL)[[c('coefficients')]])
  
  ## Score statistic ##
  if (sum(!is.finite(BetaZ.adjpsind)) == 0) 
  {
    ScT.adjpsind<-ScTt.cox(delta=Y, time=E, X=X, Z=Z.pstr, beta=BetaZ.adjpsind)
    if (is.finite(ScT.adjpsind)) {
      if (test.side==1) {
        signorm.adjpsind <- ScT.adjpsind >= qnorm(1-alpha)
      } else if (test.side==2) {
        signorm.adjpsind <- ScT.adjpsind <= qnorm(alpha)
      }
    } else {
      ScT.adjpsind <- NA
      signorm.adjpsind <- NA
    }
   
  } else
  {
    ScT.adjpsind <- NA
    signorm.adjpsind <- NA
  }
  
  ## Estimation model
  padj.fit <- survival::coxph.fit(cbind(X, Z.pstr), survival::Surv(E, Y),
                                   strata=NULL, 
                                   control=cox.ctrl, 
                                  method=tie.method, rownames=NULL)
  
  BetaX.adjpsind <- coef(padj.fit)[[1]]
  cover.adjpsind <- ci.coverage(padj.fit, truth=truth.logHR.X)
  
  #### [Adj PS B-Splines] Regression on Propensity Score (includes confounders)
  #### B-splines 

  ## NULL Model
  BetaZ.adjpsbs <- as.numeric(survival::coxph.fit(Z.bspline, 
                                            survival::Surv(E, Y),
                                            strata=NULL,
                                            control=cox.ctrl,
                                            method=tie.method,
                                            rownames=NULL)[[c('coefficients')]])
  
  ## Score statistic ##
  if (sum(!is.finite(BetaZ.adjpsbs)) == 0) {
    ScT.adjpsbs <- ScTt.cox(delta=Y, time=E, X=X, Z=Z.bspline, beta=BetaZ.adjpsbs)
    if (is.finite(ScT.adjpsbs)) {
      if (test.side==1) {
        signorm.adjpsbs <- ScT.adjpsbs >= qnorm(1-alpha)
      } else if (test.side==2) {
        signorm.adjpsbs <- ScT.adjpsbs <= qnorm(alpha)
      }
    } else {
      ScT.adjpsbs <- NA
      signorm.adjpsbs <- NA
    }
   
  } else {
    ScT.adjpsbs <- NA
    signorm.adjpsbs <- NA
  }
  
  ## Estimation model
  bspline.fit <- survival::coxph.fit(cbind(X, Z.bspline),survival::Surv(E, Y), 
                                      strata=NULL, 
                                      control=cox.ctrl, 
                                      method=tie.method,rownames=NULL)
  BetaX.adjpsbs <- coef(bspline.fit)[[1]]
  cover.adjpsbs <- ci.coverage(bspline.fit, truth=truth.logHR.X)
  
  #### [Stratify Site Adj Confounders] Statify on Site and Regress on 
  #### categorical confounders

  ## NULL Model
  BetaZ.strS.adjcon <- as.numeric(survival::coxph.fit(Z.nosite, 
                                          survival::Surv(E, Y),
                                          strata=S,
                                          control=cox.ctrl,
                                          method=tie.method,
                                          rownames=NULL)[[c('coefficients')]])
  
  ## Score statistic ##
  if (sum(!is.finite(BetaZ.strS.adjcon)) == 0) {
    
    ScT.strS.adjcon<-ScTt.cox.strat(delta=Y,time=E,X=X,Z=Z.nosite,
                                    beta=BetaZ.strS.adjcon,strata=S)
    if (is.finite(ScT.strS.adjcon)) {
      if (test.side==1) {
        signorm.strS.adjcon <- ScT.strS.adjcon >= qnorm(1-alpha)
      } else if (test.side==2) {
        signorm.strS.adjcon <- ScT.strS.adjcon <= qnorm(alpha)
      }
    } else {
      ScT.strS.adjcon <- NA
      signorm.strS.adjcon <- NA
    }
   
  } else {
    ScT.strS.adjcon <- NA
    signorm.strS.adjcon <- NA
  }
  
  ## Estimation model
  strS.adjcon.fit <- survival::coxph.fit(cbind(X,Z.nosite), survival::Surv(E, Y),
                                         strata=S, 
                                         control=cox.ctrl, 
                                         method=tie.method,
                                         rownames=NULL)
  
  BetaX.strS.adjcon <- coef(strS.adjcon.fit)[[1]]
  cover.strS.adjcon <- ci.coverage(strS.adjcon.fit, truth=truth.logHR.X)
  
  #### [Stratify PS] Stratify on Propensity Score (includes confounders 
  #### and site) categories 

  ## Score statistic ##
  ## No confounders in model so beta under the null is 0
  BetaZ.null <- 0
  ScT.strPS <- ScTt.cox.strat(delta=Y,time=E,X,Z=rep(1,length(X)),
                              beta=BetaZ.null,strata=PSstr)
  if (is.finite(ScT.strPS)) {
    if (test.side==1) {
      signorm.strPS <- ScT.strPS >= qnorm(1-alpha)
    } else if (test.side==2) {
      signorm.strPS <- ScT.strPS <= qnorm(alpha)
    }
  } else {
    signorm.strPS <- NA
    ScT.strPS <- NA
  }
  
  ## Estimation model
  psstrat.fit <- survival::coxph.fit(as.matrix(as.numeric(X)), 
                                     survival::Surv(E, Y),
                                     strata=PSstr, 
                                     control=cox.ctrl, 
                                     method=tie.method,rownames=NULL)
  
  BetaX.strPS <- coef(psstrat.fit)[[1]]
  cover.strPS <- ci.coverage(psstrat.fit, truth=truth.logHR.X)
  
  #####################################################
  #### SITE-SPECIFIC METHODS (DISTRIBUTED SETTING) ####
  
  #### [Adj Site-PS Indicators] Regression on Site-Specific Propensity Score
  #### Indicators and adjust for site and interactions with site ##

  BetaZ.adjpsind.site <- as.numeric(survival::coxph.fit(Z.pstr.site, 
                                            survival::Surv(E, Y),
                                            strata=NULL, 
                                            control=cox.ctrl, 
                                            method=tie.method,
                                            rownames=NULL)[[c('coefficients')]])
  
  ## Score statistic ##
  if (sum(!is.finite(BetaZ.adjpsind.site)) == 0) {
    
    ScT.adjpsind.site <- ScTt.cox(delta=Y, time=E, X=X, Z=Z.pstr.site, 
                                  beta=BetaZ.adjpsind.site)
    if (is.finite(ScT.adjpsind.site)) {
      if (test.side==1) {
        signorm.adjpsind.site <- ScT.adjpsind.site >= qnorm(1-alpha)
      } else if (test.side==2) {
        signorm.adjpsind.site <- ScT.adjpsind.site <= qnorm(alpha)
      }
    } else {
      ScT.adjpsind.site <- NA
      signorm.adjpsind.site <- NA
    }
    
  } else {
    ScT.adjpsind.site <- NA
    signorm.adjpsind.site <- NA
  }
  
  ## Estimation model site ##
  padj.site.fit <- survival::coxph.fit(cbind(X, Z.pstr.site), survival::Surv(E, Y),
                                       strata=NULL, 
                                       control=cox.ctrl, method=tie.method, 
                                       rownames=NULL)
  
  BetaX.adjpsind.site <- coef(padj.site.fit)[[1]]
  cover.adjpsind.site <- ci.coverage(padj.site.fit, truth=truth.logHR.X)
  
  #### [Adj Site-PS B-splines] Regression on Site-Specific Propensity Score 
  #### B-Splines and adjust for site and interactions with site 
 
  BetaZ.adjpsbs.site <- as.numeric(survival::coxph.fit(Z.bspline.site, 
                                            survival::Surv(E, Y),
                                            strata=NULL, 
                                            control=cox.ctrl,
                                            method=tie.method, 
                                            rownames=NULL)[[c('coefficients')]])
  
  # Score statistic #
  if (sum(!is.finite(BetaZ.adjpsbs.site)) == 0) {
    
    ScT.adjpsbs.site <- ScTt.cox(delta=Y, time=E, X=X, Z=Z.bspline.site, 
                                 beta=BetaZ.adjpsbs.site)
    if (is.finite(ScT.adjpsbs.site)) {
      if (test.side==1) {
        signorm.adjpsbs.site <- ScT.adjpsbs.site >= qnorm(1-alpha)
      } else if (test.side==2) {
        signorm.adjpsbs.site <- ScT.adjpsbs.site <= qnorm(alpha)
      }
    } else {
      ScT.adjpsbs.site <- NA
      signorm.adjpsbs.site <- NA
    }
    
  } else {
    ScT.adjpsbs.site <- NA
    signorm.adjpsbs.site <- NA
  }
  
  
  ## Estimation model ##
  bspline.site.fit <- survival::coxph.fit(cbind(X, Z.bspline.site),survival::Surv(E, Y), 
                                           strata=NULL, 
                                           control=cox.ctrl, 
                                           method=tie.method,rownames=NULL)
  BetaX.adjpsbs.site <- coef(bspline.site.fit)[[1]]
  cover.adjpsbs.site <- ci.coverage(bspline.site.fit, truth=truth.logHR.X)
  
  #### [Stratify Site Adj Site-PS B-splines] Stratify on Site and regress on 
  #### Site-Specific Propensity scores B-Splines and include interactions with 
  #### site
  ############################
  BetaZ.strS.adjpsbs.site <- try(as.numeric(survival::coxph.fit(Z.bspline.site.str, 
                                            survival::Surv(E, Y),
                                            strata=S, 
                                            control=cox.ctrl, 
                                            method=tie.method, 
                                            rownames=NULL)[[c('coefficients')]]),
                                 silent = TRUE)
  
  if (inherits(BetaZ.strS.adjpsbs.site, "try-error")) {
    BetaZ.strS.adjpsbs.site <- Inf
  }
  
  # Score statistic #
  if (sum(!is.finite(BetaZ.strS.adjpsbs.site)) == 0) {
    
    ScT.strS.adjpsbs.site <- ScTt.cox.strat(delta=Y, time=E, X=X, 
                                            Z=Z.bspline.site.str, 
                                            beta=BetaZ.strS.adjpsbs.site,strata=S)
    if (is.finite(ScT.strS.adjpsbs.site)) {
      if (test.side==1) {
        signorm.strS.adjpsbs.site <- ScT.strS.adjpsbs.site >= qnorm(1-alpha)
      } else if (test.side==2) {
        signorm.strS.adjpsbs.site <- ScT.strS.adjpsbs.site <= qnorm(alpha)
      }
    } else {
      ScT.strS.adjpsbs.site <- NA
      signorm.strS.adjpsbs.site <- NA
    }

  } else {
    ScT.strS.adjpsbs.site <- NA
    signorm.strS.adjpsbs.site <- NA
  }
  
  ## Estimation model ##
  ##### Rita: no S main effect in here! Only S*bspline
  bspline.site.str.fit <- try(survival::coxph.fit(cbind(X, Z.bspline.site.str),survival::Surv(E, Y), 
                                              strata=S, 
                                              control=cox.ctrl,
                                              method=tie.method,rownames=NULL),
                              silent = TRUE )
  if (inherits(bspline.site.str.fit, "try-error")) {
    bspline.site.str.fit <- NULL
  }
  
  if (!is.null(bspline.site.str.fit)) {
    BetaX.strS.adjpsbs.site <- coef(bspline.site.str.fit)[[1]]
    cover.strS.adjpsbs.site <- ci.coverage(bspline.site.str.fit, 
                                           truth=truth.logHR.X)
  } else {
    BetaX.strS.adjpsbs.site <- NA
    cover.strS.adjpsbs.site <- c(NA,NA,NA)
  }
  
  
  #### [Stratify Site+Site-PS] Stratify on Site and regress on Site-Specific 
  #### Propensity scores B-Splines and include interactions with site 
  
  ## Score statistic ##
  ## No confounders in model so no beta under the null is 0
  BetaZ.null <- 0
  ScT.strSps.site <- ScTt.cox.strat(delta=Y,time=E,X,Z=rep(1,length(X)),
                                    beta=BetaZ.null,strata=PSstr.by.site)
  if (is.finite(ScT.strSps.site)) {
    if (test.side==1) {
      signorm.strSps.site <- ScT.strSps.site >= qnorm(1-alpha)
    } else if (test.side==2) {
      signorm.strSps.site <- ScT.strSps.site <= qnorm(alpha)
    }
  } else {
    ScT.strSps.site <- NA
    signorm.strSps.site <- NA
  }
  

  ## Estimation model ##
  psstrat.site.fit <- survival::coxph.fit(as.matrix(as.numeric(X)), 
                                          survival::Surv(E, Y),
                                          strata=PSstr.by.site, 
                                          control=cox.ctrl, 
                                          method=tie.method,rownames=NULL)
  
  BetaX.strSps.site <- coef(psstrat.site.fit)[[1]]
  cover.strSps.site <- ci.coverage(psstrat.site.fit, truth=truth.logHR.X)
  
  return(list(
    N.exp=sum(Y[X==1]), 
    N.unexp=sum(Y[X==0]),

    Signorm.adjcon=signorm.adjcon,
    Signorm.adjpsind=signorm.adjpsind,
    Signorm.adjpsbs=signorm.adjpsbs,
    Signorm.strS.adjcon=signorm.strS.adjcon,
    Signorm.strPS=signorm.strPS,
    Signorm.adjpsind.site=signorm.adjpsind.site,
    Signorm.adjpsbs.site=signorm.adjpsbs.site,
    Signorm.strS.adjpsbs.site=signorm.strS.adjpsbs.site,
    Signorm.strSps.site=signorm.strSps.site,
    Signorm.mh.samp = signorm.mh.samp, 
    Signorm.mh.var = signorm.mh.var,
    Signorm.mh.bs.samp = signorm.mh.bs.samp, 
    Signorm.mh.bs.var = signorm.mh.bs.var,
    
    Score.adjcon=ScT.adjcon,
    Score.adjpsind=ScT.adjpsind,
    Score.adjpsbs=ScT.adjpsbs,
    Score.strS.adjcon=ScT.strS.adjcon,
    Score.strPS=ScT.strPS,
    Score.adjpsind.site=ScT.adjpsind.site,
    Score.adjpsbs.site=ScT.adjpsbs.site,
    Score.strS.adjpsbs.site=ScT.strS.adjpsbs.site,
    Score.strSps.site=ScT.strSps.site,
    Score.mh.samp = Score.mh.samp,
    Score.mh.var = Score.mh.var,
    Score.mh.bs.samp = Score.mh.bs.samp,
    Score.mh.bs.var = Score.mh.bs.var,
    
    BetaX.unadj=BetaX.unadj,
    BetaX.adjcon=BetaX.adjcon,
    BetaX.adjpsind=BetaX.adjpsind,
    BetaX.adjpsbs=BetaX.adjpsbs,
    BetaX.strS.adjcon=BetaX.strS.adjcon,
    BetaX.strPS=BetaX.strPS,
    BetaX.adjpsind.site=BetaX.adjpsind.site,
    BetaX.adjpsbs.site=BetaX.adjpsbs.site,
    BetaX.strS.adjpsbs.site=BetaX.strS.adjpsbs.site,
    BetaX.strSps.site=BetaX.strSps.site,
    BetaX.marg.nonstrat= BetaX.marg.nonstrat,
    BetaX.marg.strat=BetaX.marg.strat,
    BetaX.mh.samp = BetaX.mh.samp,
    BetaX.mh.var = BetaX.mh.var,
    BetaX.mh.bs.samp = BetaX.mh.bs.samp,
    BetaX.mh.bs.var = BetaX.mh.bs.var,
    BetaX.marg = BetaX.marg,

    Cover.unadj=cover.unadj,
    Cover.adjcon=cover.adjcon,
    Cover.adjpsind=cover.adjpsind,
    Cover.adjpsbs=cover.adjpsbs,
    Cover.strS.adjcon=cover.strS.adjcon,
    Cover.strPS=cover.strPS,
    Cover.adjpsind.site=cover.adjpsind.site,
    Cover.adjpsbs.site=cover.adjpsbs.site,
    Cover.strS.adjpsbs.site=cover.strS.adjpsbs.site,
    Cover.strSps.site=cover.strSps.site,
    Cover.mh.samp = cover.mh.samp,
    Cover.mh.var = cover.mh.var,
    Cover.mh.bs.samp = cover.mh.bs.samp,
    Cover.mh.bs.var = cover.mh.bs.var,
    Cover.marg = cover.marg,
    
    BetaZ.adjcon=BetaZ.adjcon, 
    BetaZ.adjpsind=BetaZ.adjpsind, 
    BetaZ.adjpsbs=BetaZ.adjpsbs, 
    BetaZ.strS.adjcon=BetaZ.strS.adjcon, 
    BetaZ.adjpsind.site=BetaZ.adjpsind.site,
    BetaZ.adjpsbs.site=BetaZ.adjpsbs.site,  
    BetaZ.strS.adjpsbs.site=BetaZ.strS.adjpsbs.site,
    PS.model.pooled = PS.list$Beta.pooled,
    PS.model.site = PS.list$Beta.site
  )) 
}


###################################################################################
#### methodDist
#### DESCRIPTION:  
#### INPUTS:
#### 
####
#### OUTPUTS:
#### 
methodDist <- function(data.full,
                       cov.cols,
                       cox.ctrl = cox.ctrl,
                       tie.method = tie.method,
                       truth.logHR.X,
                       n.PSstr = 5,
                       bs.df = 5,
                       bs.degree = 3,
                       time.int = 7,test.side=1,alpha=0.05)
{
  data.full <- data.full$Data.Simulated
  X <- data.full$X
  Y <- data.full$Y
  E <- data.full$E
  E.dis <- ceiling(E / time.int) * time.int
  S <- data.full$site
  Z <-
    cbind(as.matrix(data.full[, cov.cols]), model.matrix( ~ factor(S))[, -1])
  Z.nosite <- as.matrix(data.full[, cov.cols])
  
  
  #### Compute pooled and site-specific propensity scores
  PS.list <-
    comp.pscore(X = data.full$X, Z = data.full[, cov.cols], S = data.full$site)
  
  #### Set pooled p-scores
  PS.Pool <- PS.list$PS.pooled
  #### Set site-specific p-scores
  PS.Site.list <- PS.list$PS.site
  #### Count the number of sites
  n.sites <- length(PS.Site.list)
  
  #### Site-specific p-scores, p-score categories & b-splines
  PS.Site <- PSstr.Site <- NULL
  
  
  for (i in 1:n.sites)
  {
    PS.Site <- c(PS.Site, PS.Site.list[[i]])
    PSstr.Site <-
      c(
        PSstr.Site,
        cut(
          PS.Site.list[[i]],
          include.lowest = T,
          breaks = quantile(PS.Site.list[[i]],
                            seq(0, 1, 1 / n.PSstr), names =
                              F),
          labels = FALSE
        )
      )
  }
  

  ## Pooled p-score categories
  PSstr <- cut(PS.Pool,
               include.lowest = T,
               breaks = quantile(PS.Pool, seq(0, 1, 1 / n.PSstr), names =
                                   F))
  ## Site-specific p-score strata
  PSstr.by.site <- as.character(PSstr.Site + S * 100000)
  
  Z.pstr <- model.matrix( ~ PSstr)[, c(-1)]
  
  ## site-specific PS.spline/strata
  
  Z.pstr.site <- model.matrix( ~ factor(PSstr.Site) * factor(S))[, -1]

  ## [Unadj] Cox PH regression without adjusting for confounders ##
  unadj.fit <- survival::coxph.fit(as.matrix(as.numeric(X)), survival::Surv(E.dis, Y), 
                                   strata=NULL, control=cox.ctrl, method=tie.method, rownames=NULL)
  BetaX.unadj <- coef(unadj.fit)
  cover.unadj <- ci.coverage(unadj.fit, truth=truth.logHR.X)
  
  #### [Adj Confounders+Site] Regression on categorical confounders and site
 
  ## NULL Model
  BetaZ.adjcon <-
    as.numeric(
      survival::coxph.fit(
        as.matrix(Z),
        survival::Surv(E.dis, Y),
        strata = NULL,
        control = cox.ctrl,
        method = tie.method,
        rownames = NULL
      )[[c('coefficients')]]
    )
  
  ## Score statistic ##
  if (sum(is.na(BetaZ.adjcon) | is.infinite(BetaZ.adjcon)) == 0) {
    
    ScT.adjcon <- ScTt.cox(delta = Y,time = E.dis,X = X, Z = Z,beta = BetaZ.adjcon)
    
    if (test.side==1) {
      signorm.adjcon <- ScT.adjcon >= qnorm(1-alpha)
    } else if (test.side==2) {
      signorm.adjcon <- ScT.adjcon <= qnorm(alpha)
    }
  } else {
    ScT.adjcon <- -99
    signorm.adjcon <- NA
  }
  
  ## Estimation model
  adjcon.fit <- survival::coxph.fit(
    cbind(X, Z),
    survival::Surv(E.dis, Y),
    strata = NULL,
    control = cox.ctrl,
    method = tie.method,
    rownames = NULL
  )
  BetaX.adjcon <- coef(adjcon.fit)[[1]]
  cover.adjcon <- ci.coverage(adjcon.fit, truth = truth.logHR.X)
  
  
  #####################################################
  #### SITE-SPECIFIC METHODS (DISTRIBUTED SETTING) ####
  
  #### [Adj Site-PS Indicators] Regression on Site-Specific Propensity Score
  #### Indicators and adjust for site and interactions with site ##
  
  BetaZ.adjpsind.site <- as.numeric(
    survival::coxph.fit(
      Z.pstr.site,
      survival::Surv(E.dis, Y),
      strata = NULL,
      control = cox.ctrl,
      method = tie.method,
      rownames = NULL
    )[[c('coefficients')]]
  )
  
  ## Score statistic ##
  if (sum(is.na(BetaZ.adjpsind.site) | is.infinite(BetaZ.adjpsind.site)) == 0) {
    
    ScT.adjpsind.site <- ScTt.cox(delta = Y,time = E.dis,X = X,Z = Z.pstr.site,beta = BetaZ.adjpsind.site)
    
    if (test.side==1) {
      signorm.adjpsind.site <- ScT.adjpsind.site >= qnorm(1-alpha)
    } else if (test.side==2) {
      signorm.adjpsind.site <- ScT.adjpsind.site <= qnorm(alpha)
    }
  } else {
    ScT.adjpsind.site <- -99
    signorm.adjpsind.site <- NA
  }
  
  ## Estimation model site ##
  padj.site.fit <-
    survival::coxph.fit(
      cbind(X, Z.pstr.site),
      survival::Surv(E.dis, Y),
      strata = NULL,
      control = cox.ctrl,
      method = tie.method,
      rownames = NULL
    )
  
  BetaX.adjpsind.site <- coef(padj.site.fit)[[1]]
  cover.adjpsind.site <-
    ci.coverage(padj.site.fit, truth = truth.logHR.X)

  
  #### [Stratify Site+Site-PS] Stratify on Site and on Site-Specific
  #### Propensity scores 
  
  ## Score statistic ##
  ## No confounders in model so no beta under the null is 0
  BetaZ.null <- 0
  ScT.strSps.site <- ScTt.cox.strat(delta = Y,time = E.dis,X,
                                    Z = rep(1, length(X)),beta = BetaZ.null,strata = PSstr.by.site)
  
  if (test.side==1) {
    signorm.strSps.site <- ScT.strSps.site >= qnorm(1-alpha)
  } else if (test.side==2) {
    signorm.strSps.site <- ScT.strSps.site <= qnorm(alpha)
  }
  
  ## Estimation model ##
  psstrat.site.fit <- survival::coxph.fit(
    as.matrix(as.numeric(X)),
    survival::Surv(E.dis, Y),
    strata = PSstr.by.site,
    control = cox.ctrl,
    method = tie.method,
    rownames = NULL
  )
  
  BetaX.strSps.site <- coef(psstrat.site.fit)[[1]]
  cover.strSps.site <-
    ci.coverage(psstrat.site.fit, truth = truth.logHR.X)
  
  return(list(
    N.exp=sum(Y[X==1]), 
    N.unexp=sum(Y[X==0]),
    
    Signorm.adjcon=signorm.adjcon,
    Signorm.adjpsind.site=signorm.adjpsind.site,
    Signorm.strSps.site=signorm.strSps.site,
    
    Score.adjcon=ScT.adjcon,
    Score.adjpsind.site=ScT.adjpsind.site,
    Score.strSps.site=ScT.strSps.site,
    
    BetaX.unadj=BetaX.unadj,
    BetaX.adjcon=BetaX.adjcon,
    BetaX.adjpsind.site=BetaX.adjpsind.site,
    BetaX.strSps.site=BetaX.strSps.site,
    
    Cover.unadj=cover.unadj,
    Cover.adjcon=cover.adjcon,
    Cover.adjpsind.site=cover.adjpsind.site,
    Cover.strSps.site=cover.strSps.site,
    
    BetaZ.adjcon=BetaZ.adjcon, 
    BetaZ.adjpsind.site=BetaZ.adjpsind.site,
    PS.model.pooled = PS.list$Beta.pooled,
    PS.model.site = PS.list$Beta.site
  ))
}


###################################################################################
#### methodOneOff
#### DESCRIPTION:  
#### INPUTS:
#### 
####
#### OUTPUTS:
#### 
methodOneOff <- function(data.full, cov.cols, cox.ctrl=cox.ctrl,
                         tie.method=tie.method,truth.logHR.X,
                         n.PSstr=5,bs.df=5,bs.degree=3,test.side=1,alpha=0.05)
{
  data.marg <- data.full$Data.Marginal
  data.full <- data.full$Data.Simulated
  X <- data.full$X
  Y <- data.full$Y
  E <- data.full$E
  S <- data.full$site
  Z <- cbind(as.matrix(data.full[,cov.cols]), model.matrix(~factor(S))[,-1]) 
  Z.nosite <- as.matrix(data.full[,cov.cols])
  
  ## Two interaction models for marginal hazard
  
  ## Site interacts will all Z only main effect of X
  marg.nonstrat <- survival::coxph.fit(model.matrix(~ X + Z.nosite*factor(S))[,-c(1)],
                                       survival::Surv(E, Y),
                                       strata=NULL, control=cox.ctrl,
                                       method=tie.method, rownames=NULL)
  BetaX.marg.nonstrat <- coef(marg.nonstrat)[[1]]
  cover.marg.nonstrat <- ci.coverage(marg.nonstrat, truth=truth.logHR.X)
  
  ## Site interacts will all Z but no main effect of site. Instead
  ## site is used as a stratification variable. Only main effect of X
  marg.strat <- survival::coxph.fit(model.matrix(~X + Z.nosite + Z.nosite:factor(S))[,-c(1)],
                                    survival::Surv(E, Y),
                                    strata=S, control=cox.ctrl,
                                    method=tie.method, rownames=NULL)
  BetaX.marg.strat <- coef(marg.strat)[[1]]
  cover.marg.strat <- ci.coverage(marg.strat, truth=truth.logHR.X)
  
  ## Simulated marginal cox
  marg <- survival::coxph.fit(model.matrix(~ data.marg$x + factor(data.marg$site))[,-c(1)],
                              survival::Surv(data.marg$obst, data.marg$y),
                              strata=NULL, control=cox.ctrl,
                              method=tie.method, rownames=NULL)
  BetaX.marg <- coef(marg)[[1]]
  cover.marg <- ci.coverage(marg, truth=truth.logHR.X)
  
  #### Compute pooled and site-specific propensity scores
  PS.list <- comp.pscore(X=data.full$X,Z=data.full[,cov.cols],S=data.full$site)
  
  #### Set pooled p-scores
  PS.Pool <- PS.list$PS.pooled
  #### Set site-specific p-scores
  PS.Site.list <- PS.list$PS.site
  #### Count the number of sites
  n.sites <- length(PS.Site.list)
  
  #### Site-specific p-scores, p-score categories & b-splines
  PS.Site <- PSstr.Site <- PSspline.Site <- NULL ;
  
  for(i in 1:n.sites)
  {
    PS.Site <- c(PS.Site,PS.Site.list[[i]])
    PSstr.Site <- c(PSstr.Site, cut(PS.Site.list[[i]],include.lowest=T,
                                    breaks=quantile(PS.Site.list[[i]],
                                                    seq(0,1,1/n.PSstr),names=F),
                                    labels=FALSE))
    PSspline.Site <- rbind(PSspline.Site, splines::bs(PS.Site.list[[i]],
                                                      df=bs.df, degree=bs.degree))
  }
  
  ## Pooled b-splines
  PSspline <- splines::bs(PS.Pool, df=bs.df, degree=bs.degree)
  ## Pooled p-score categories
  PSstr <- cut(PS.Pool, include.lowest=T, 
               breaks=quantile(PS.Pool, seq(0,1,1/n.PSstr), names=F))
  ## Site-specific p-score strata
  PSstr.by.site <- as.character(PSstr.Site+S*100000)
  
  Z.pstr <- model.matrix(~ PSstr)[,c(-1)]
  Z.bspline <- PSspline
  
  ## site-specific PS.spline/strata
  
  Z.pstr.site <- model.matrix(~ factor(PSstr.Site)*factor(S))[,-1]
  Z.bspline.site <- model.matrix(~ PSspline.Site*factor(S))[,-1]
  Z.bspline.site.str <- Z.bspline.site[,-c((bs.df+1):(bs.df + n.sites - 1))]
  
  ###################################
  #### Cox PH Regression Methods #### 
  
  ## [Unadj] Cox PH regression without adjusting for confounders ##
  unadj.fit <- survival::coxph.fit(as.matrix(as.numeric(X)), survival::Surv(E, Y), 
                                   strata=NULL, control=cox.ctrl, method=tie.method, rownames=NULL)
  BetaX.unadj <- coef(unadj.fit)
  cover.unadj <- ci.coverage(unadj.fit, truth=truth.logHR.X)
  
  
  #######################################################
  #### POOLED METHODS (E.G. NON-DISTRIBUTED SETTING) ####
  
  #### [Adj Confounders+Site] Regression on categorical confounders and site
  
  ## NULL Model
  BetaZ.adjcon <- as.numeric(survival::coxph.fit(as.matrix(Z), survival::Surv(E, Y),
                                                 strata=NULL,control=cox.ctrl,
                                                 method=tie.method,
                                                 rownames=NULL)[[c('coefficients')]])
  
  ## Score statistic ##
  if (sum(!is.finite(BetaZ.adjcon)) == 0) {
    ScT.adjcon<-ScTt.cox(delta=Y, time=E, X=X, Z=Z, beta=BetaZ.adjcon)
    if (is.finite(ScT.adjcon)) {
      if (test.side==1) {
        signorm.adjcon <- ScT.adjcon >= qnorm(1-alpha)
      } else if (test.side==2) {
        signorm.adjcon <- ScT.adjcon <= qnorm(alpha)
      }
    } else {
      ScT.adjcon <- NA
      signorm.adjcon <- NA
    }
    
  } else {
    ScT.adjcon <- NA
    signorm.adjcon <- NA
  }
  
  ## Estimation model
  adjcon.fit <- survival::coxph.fit(cbind(X,Z), survival::Surv(E, Y),
                                    strata=NULL, control=cox.ctrl, method=tie.method,rownames=NULL)
  BetaX.adjcon <- coef(adjcon.fit)[[1]]
  cover.adjcon <- ci.coverage(adjcon.fit, truth=truth.logHR.X)
  
  
  #### [Stratify PS] Stratify on Propensity Score (includes confounders 
  #### and site) categories 
  
  ## Score statistic ##
  ## No confounders in model so beta under the null is 0
  BetaZ.null <- 0
  ScT.strPS <- ScTt.cox.strat(delta=Y,time=E,X,Z=rep(1,length(X)),
                              beta=BetaZ.null,strata=PSstr)
  if (is.finite(ScT.strPS)) {
    if (test.side==1) {
      signorm.strPS <- ScT.strPS >= qnorm(1-alpha)
    } else if (test.side==2) {
      signorm.strPS <- ScT.strPS <= qnorm(alpha)
    }
  } else {
    signorm.strPS <- NA
    ScT.strPS <- NA
  }
  
  ## Estimation model
  psstrat.fit <- survival::coxph.fit(as.matrix(as.numeric(X)), 
                                     survival::Surv(E, Y),
                                     strata=PSstr, 
                                     control=cox.ctrl, 
                                     method=tie.method,rownames=NULL)
  
  BetaX.strPS <- coef(psstrat.fit)[[1]]
  cover.strPS <- ci.coverage(psstrat.fit, truth=truth.logHR.X)
  
  
  
  return(list(
    N.exp=sum(Y[X==1]), 
    N.unexp=sum(Y[X==0]),
    
    Signorm.adjcon=signorm.adjcon,
    Signorm.strPS=signorm.strPS,
    
    Score.adjcon=ScT.adjcon,
    Score.strPS=ScT.strPS,
    
    BetaX.unadj=BetaX.unadj,
    BetaX.adjcon=BetaX.adjcon,
    BetaX.strPS=BetaX.strPS,
    BetaX.marg.nonstrat=BetaX.marg.nonstrat,
    BetaX.marg.strat=BetaX.marg.strat,
    
    Cover.unadj=cover.unadj,
    Cover.adjcon=cover.adjcon,
    Cover.strPS=cover.strPS,
    cover.marg.nonstrat,
    cover.marg.strat,
    
    BetaZ.adjcon=BetaZ.adjcon, 
    
    PS.model.pooled = PS.list$Beta.pooled,
    PS.model.site = PS.list$Beta.site
  )) 
}

## Mantel-Hanzel Estimator
##
##
##

coxmh <- function (E,Y,X,Z,S,test.side=1,alpha=0.05){

  n.sites <- length(table(S))
  
  dat.site <- sapply(1:n.sites,function(XX) NULL)
  
  if (mode(S)=='character') {
    names.site <- names(sort(table(S),decreasing=TRUE))
    for (i in 1:n.sites) {
      dat.site[[i]] <- list('E'=E[S==names.site[[i]]],'Y'=Y[S==names.site[[i]]],
                            'X'=X[S==names.site[[i]]],'Z'=Z[S==names.site[[i]],])
    }
  } else {
    for (i in 1:n.sites) {
      dat.site[[i]] <- list('E'=E[S==i],'Y'=Y[S==i],'X'=X[S==i],'Z'=Z[S==i,])
    }
  }

  
  ### Null model for MH
  BetaZ.mh <- lapply(dat.site, function(XX) {
      tryCatch(
    {
      coef(survival::coxph.fit(as.matrix(XX$Z), survival::Surv(XX$E, XX$Y), 
                               strata=NULL, control=cox.ctrl, method=tie.method, 
                               rownames=NULL))
    },
    error = function(e)  NULL)
  })

  
  ### Model for site-specific effect estimates.
  Model.mh <- lapply(dat.site, function(XX) {
    tryCatch(
      {
        survival::coxph.fit(cbind(XX$X, as.matrix(XX$Z)), survival::Surv(XX$E, XX$Y), 
                          strata=NULL, control=cox.ctrl, 
                          method=tie.method, rownames=NULL)
        }, 
      error = function(e) NULL )
    })
    
  
  ### If both Null Model and Estimation Modle fit without error,
  ### add results to site-specific MH vectors: Score Stat, BetaX, Var(BetaX) 
  ### and site sample size.
  ScT.mh <- BetaX.mh <- Var.mh <- N.mh <- NULL
  for (i in 1:n.sites) {
    if (!is.null(Model.mh[[i]]) & !is.null(BetaZ.mh[[i]])) {
      ScT.mh<-c(ScT.mh, ScTt.cox(delta=dat.site[[i]]$Y, time=dat.site[[i]]$E,
                                 X=dat.site[[i]]$X, Z=dat.site[[i]]$Z,
                                 beta=BetaZ.mh[[i]]))
      BetaX.mh <- c(BetaX.mh, coef(Model.mh[[i]])[[1]])
      Var.mh <- c(Var.mh, Model.mh[[i]]$var[1,1])
      N.mh <- c(N.mh, length(dat.site[[i]]$Y))
      
   } #  else {
    #   ScT.mh<-c(ScT.mh, NA)
    #   BetaX.mh <- c(BetaX.mh, NA)
    #   Var.mh <- c(Var.mh, NA)
    #   N.mh <- c(N.mh, NA)
    # } 
  }

  ## Combined MH Score estimate
  if (!is.null(ScT.mh)) {
    ## Sample size weighted estimate
    N.mh <- as.numeric(N.mh)
    ScT.mh.comb <- ScT.mh%*%N.mh/sqrt(sum(N.mh^2))
    BetaX.mh.comb <- BetaX.mh%*%N.mh/sum(N.mh)
    Var.mh.comb <- Var.mh%*%(N.mh*N.mh/sum(N.mh)^2)
    ci.mh <- BetaX.mh.comb + c(-1.959964, 1.959964)*sqrt(Var.mh.comb)
    
    if (test.side==1) {
      signorm.mh.comb <- ScT.mh.comb >= qnorm(1-alpha)
    } else if (test.side==2) {
      signorm.mh.comb <- ScT.mh.comb <= qnorm(alpha)
    }
    
    
    ## Inverse variance weighted estimate
    ScT.mh.comb.var <- ScT.mh%*%(1/Var.mh)/sqrt(sum(as.numeric((1/Var.mh)*(1/Var.mh))))
    BetaX.mh.comb.var <- BetaX.mh%*%(1/Var.mh)/sum(1/Var.mh)
    Var.mh.comb.var <- 1/sum(1/Var.mh)
    ci.mh.var <- BetaX.mh.comb.var + c(-1.959964, 1.959964)*sqrt(Var.mh.comb.var)
    
    if (test.side==1) {
      signorm.mh.comb.var <- ScT.mh.comb.var >= qnorm(1-alpha)
    } else if (test.side==2) {
      signorm.mh.comb.var <- ScT.mh.comb.var <= qnorm(alpha)
    }
  } else {
    ## Sample size weighted estimate
    ScT.mh.comb <- -99
    BetaX.mh.comb <- -99
    Var.mh.comb <- -99

    ## Inverse variance weighted estimate
    ScT.mh.comb.var <- -99
    BetaX.mh.comb.var <- -99
    Var.mh.comb.var <- -99
  }
  
  ## End site-specific models for MH
  return(list(SampWeight = list(BetaX.mh.comb, Var.mh.comb, ci.mh, ScT.mh.comb, signorm.mh.comb),
              VarWeight=list(BetaX.mh.comb.var,
              Var.mh.comb.var,ci.mh.var,ScT.mh.comb.var, signorm.mh.comb.var)))
}

comp.pscore <- function(X,Z,S) {

  glm_control <- glm.control(epsilon = 1e-8, maxit = 25, trace = FALSE)
  s.dummy <- model.matrix(~factor(S))[,-1]
  fit.pooled <- glm.fit(cbind(Int=1,Z,s.dummy), X, control=glm_control, family=binomial())
  beta.pooled <- coef(fit.pooled)
  ps.pooled <- fitted(fit.pooled)
  
  n.sites <- length(table(S))
  
  if (n.sites > 1) {
    ps.site <- beta.site <- sapply(1:n.sites, function(XX) NULL)
    for (i in 1:n.sites) {
      fit.site <- NULL
      fit.site <- glm.fit(cbind(Int=1,Z[S==i,]), X[S==i], control=glm_control, family=binomial())
      ps.site[[i]] <- fitted(fit.site)
      beta.site[[i]] <- coef(fit.site)
    }
    return(list(PS.pooled=ps.pooled, Beta.pooled=beta.pooled, PS.site=ps.site,
                Beta.site=beta.site))
  }
  return(list(PS.pooled=ps.pooled, Beta.pooled=beta.pooled, PS.site=ps.pooled,
              Beta.site=beta.pooled))
}

ScTt.cox <- function(delta,time,X,Z,beta){
  delta <- as.matrix(delta)
  time <- as.matrix(time)
  X <- as.matrix(X)
  Z <- as.matrix(Z)
  
  # Unique Times
  utimes <- unique(time[delta==1])
  nutimes <- table(time[delta==1])
  
  # for logrank num and den
  num <- 0
  den <- 0
  Zbeta <- exp(Z%*%beta)
  
  for (j in 1:length(utimes)) {
    Vlt <- (time>=utimes[j])
    
    # Logrank test
    Sum1 <- t(Vlt*Zbeta)%*%X
    Sum2 <- c(t(Vlt)%*%Zbeta)
    
    num <- num+t(delta[time==utimes[j]])%*%(X[time==utimes[j],]-Sum1/Sum2)
    den <- den+sum(delta[time==utimes[j]])*((Sum1/Sum2)-(Sum1^2/Sum2^2))
  }
  
  # Logrank  
  LR<-num/sqrt(den)
  LR[is.na(LR)]<-0
  
  return(LR)
}

ci.coverage <- function(surv.fit, truth){
  
  coeff.ci <- surv.fit[['coefficients']][[1]] + c(-1.959964, 1.959964)*sqrt(surv.fit$var[1,1])
  c((truth >= coeff.ci[1]) & (truth <= coeff.ci[2]), coeff.ci)
}

ScTt.cox.strat <- function(delta,time,X,Z,beta=NULL,strata){
  
  delta <- as.matrix(delta)
  time <- as.matrix(time)
  X <- as.matrix(X)
  Z <- as.matrix(Z)
  
  S <- sort(unique(strata))
  #S <- S[table(strata,Y)[,2]>0] ##########Rita changed****************************************
  S <- S[table(strata,delta)[,2]>0]
  ScTt.S <- NULL
  
  for(i in S){
    
    #ScTt.S <- rbind(ScTt.S,ScTt.cox.more(delta[strata==i],time[strata==i],X[strata==i],Z[strata==i,],beta=beta))
    ScTt.S <- rbind(ScTt.S,ScTt.cox2(delta[strata==i],time[strata==i],X[strata==i],Z[strata==i,],beta=beta))
  }
  
  ScTt.strat <- apply(ScTt.S,2,sum)
  
  
  return(ScTt.strat[[1]]/sqrt(ScTt.strat[[2]]))
  #return(ScTt.strat[[1]]/sqrt(length(S)))
}

ScTt.cox2 <- function(delta,time,X,Z,beta) {
  
  delta <- as.matrix(delta)
  time <- as.matrix(time)
  X <- as.matrix(X)
  Z <- as.matrix(Z)
  
  # Unique Times
  utimes <- unique(time[delta==1])
  nutimes <- table(time[delta==1])
  
  # for logrank num and den
  num <- 0
  den <- 0
  Zbeta <- exp(Z%*%beta)
  
  for (j in 1:length(utimes)) {
    Vlt <- (time>=utimes[j])
    
    # Logrank test
    Sum1 <- t(Vlt*Zbeta)%*%X
    Sum2 <- c(t(Vlt)%*%Zbeta)
    
    num <- num+t(delta[time==utimes[j]])%*%(X[time==utimes[j],]-Sum1/Sum2)
    den <- den+sum(delta[time==utimes[j]])*((Sum1/Sum2)-(Sum1^2/Sum2^2))
  }
  
  # Logrank  
  #LR<-num/sqrt(den)
  #LR[is.na(LR)]<-0
  
  return(c(num,den))
}


#### 
## Extract and combine raw data results
##
extract_res <- function(res, ncores, nsim, labels, nsites) {
  temp_main <- sapply(1:length(labels), function(XX) {NULL})
  
  for (k in 1:length(labels)) {
    
    temp_labels <- sapply(1:(ncores*nsim), function(XX) {NULL})
    
    counter <- 1
    
    for (i in 1:ncores) {
      
      for (j in 1:nsim) {
        
        temp_labels[[counter]] <- res[[i]][[j]][[which(labels==labels[k])]]
        
        counter <- counter + 1
      }
      
    }
    
    if (k %in% grep("Cover.", labels)) {
      temp_main[[k]] <- do.call('rbind', temp_labels)[,1]
    } else if (k %in% grep("BetaZ.", labels)) {
      temp_main[[k]] <- do.call('cbind',temp_labels)
    } else if (k %in% grep("BetaX.",labels)) {
      temp_main[[k]] <- do.call('cbind',temp_labels)
    } else if (k %in% grep("PS.model.pooled",labels)) {
      temp_main[[k]] <- do.call('rbind',temp_labels)
    } else if (k %in% grep("PS.model.site",labels)) {
      temp_ps_site <- sapply(1:nsites, function(XX) {NULL})
      for (l in 1:nsites) {
        temp_ps_site[[l]] <- do.call("rbind", lapply(temp_labels, function(XX) XX[[l]]))
      }
      temp_main[[k]] <- temp_ps_site
    } else {
      temp_main[[k]] <- do.call('c',temp_labels)
    }  
    
  }
  names(temp_main) <- labels
  temp_main
}


make.plot <- function (XX) {
  par(mar=c(5.1,4.1,3.1,2.1))
  max_y <- 1.04*max(c(hist(XX$propensity[XX$propensity[,1]==1,2],plot=FALSE, breaks=30)$density,
                      hist(XX$propensity[XX$propensity[,1]==0,2],plot=FALSE, breaks=30)$density))
  hist(XX$propensity[XX$propensity[,1]==1,2], col=rgb(255/255,99/255,71/255,0.4), 
       xlim=c(0,1),ylim=c(0,max_y), xlab='Propensity Score', 
       ylab="Density", freq=FALSE, main='', breaks=30)
  hist(XX$propensity[XX$propensity[,1]==0,2], col=rgb(70/255,130/255,180/255,0.4), 
       xlim=c(0,1), freq=FALSE, add=TRUE, breaks=30)
  text(0.60, .95*max_y, c(toupper(XX[[1]])), cex=1.0, pos=4)
  text(0.65, .86*max_y, "Rivaroxaban", cex=1.0, pos=4)
  text(0.65, .78*max_y, "Warfarin", cex=1.0, pos=4)
  points(.63, 0.86*max_y, pch=15, cex=1.5,col=rgb(255/255,99/255,71/255,0.4))
  points(.63, 0.78*max_y, pch=15, cex=1.5,col=rgb(70/255,130/255,180/255,0.4))
  text(0.65,0.65*max_y, paste0("c-statistic (AUC): ",XX$cStat),pos=4)
}


## Function for cross validation
cross.val.weib <- function(Y,E,Z,k,n.pts,strat.var.event=NULL,strat.var.cens=NULL) {
  
  val.group <- sample(1:k,size=length(Y), prob=rep(1,times=k)/k,
                      replace=TRUE)
  E.sim <- as.numeric(rep(NA,length(val.group)))
  Y.sim <- as.numeric(rep(NA,length(val.group)))
  
  for(i in 1:k) {
    indx <- val.group!=i
    simulated <- NULL
    test.cens <- NULL
    test.event<- NULL
    pt.masses <- NULL
   
    ## Fit model to training data
    pt.props <- tail(sort(table(E[indx & Y==0])/length(E[indx & Y==0])),n.pts)
    pt.masses <- as.numeric(names(pt.props))
    
    # ## Proportion with a prescription bump amongst those censored
    # pt.props <- NULL
    # for( i in 1:length(pt.masses))
    # {
    #   pt.props<-cbind(pt.props,mean(E[Y==0 & indx]==pt.masses[i]))
    # }
    # pt.props <- data.frame(pt.props)
    # names(pt.props) <- as.character(pt.masses)
    print(pt.props)
    
    ind <- (!E%in%pt.masses)
    mode(Z) <- "integer"
    if (is.null(strat.var.cens)) {
      test.cens <- survival::survreg(survival::Surv(E[indx & ind],
                                   1-Y[indx & ind]) ~ Z[indx & ind,], 
                                   dist="weib")

    } else {
      test.cens <- survival::survreg(survival::Surv(E[indx & ind],
                                     1-Y[indx & ind]) ~ Z[indx & ind,] + 
                                       strata(strat.var.cens[indx & ind]), 
                                     dist="weib")

    }
    
    if (is.null(strat.var.event)) {

      test.event <- tryCatch(survival::survreg(survival::Surv(E[indx], 
                                                     Y[indx]) ~ Z[indx,], 
                                      dist="weib"),
                             error=function(e) NULL)
      if (is.null(test.event)>0) {
        test.event <- survival::survreg(survival::Surv(E[indx], 
                                                       Y[indx]) ~ Z[indx,], 
                                        dist="exp")
      }
    } else {

      test.event <- tryCatch(survival::survreg(survival::Surv(E[indx], 
                                                     Y[indx]) ~ Z[indx,] + 
                                        strata(strat.var.event[indx]), 
                                      dist="weib"),
                             error=function(e) NULL)
    if (is.null(test.event)>0) {
      test.event <- survival::survreg(survival::Surv(E[indx], 
                                                     Y[indx]) ~ Z[indx,] + 
                                        strata(strat.var.event[indx]), 
                                      dist="exp")
    }
      
    }
    
    
    ## Simulate values using on test data using fit from training set
    simulated <- generate.data(n=length(Y[!indx]),
                               user.data=cbind(1,Z)[!indx,],
                               coef.cens=test.cens$coefficients,
                               scale.cens=test.cens$scale,
                               coef.event=test.event$coefficients,
                               scale.event=test.event$scale,
                               censtype="covbump", trunc=NULL,P.presc.topK=pt.props,
                               prescription.mode.topK=pt.masses,
                               method=4,noX=TRUE,strat.var.cens=strat.var.cens[!indx],
                               strat.var.event=strat.var.event[!indx])

    ## Create new follow up time for the test set
    E.sim[!indx] <- simulated$E
    Y.sim[!indx] <- simulated$Y
  }
  
  list(E=E,Y=Y,E.sim=E.sim,Y.sim=Y.sim,
       Test=suppressWarnings(ks.test(E,E.sim, alternative="two.sided")))
}

bisection2 <- function(x.l, x.r, tol=1e-09,
                       logHR.X=logHR.X, control.rate=control.rate,
                       n=60000, P=P, Common.P=Common.P, coef.AonB=coef.AonB, coef.XonZ=coef.XonZ,
                       coef.cens, scale.cens,
                       coef.event, scale.event,
                       censtype="simple", trunc=366,
                       P.presc.topK=NULL, prescription.mode.topK=NULL,
                       Corr.norm=NULL, Quant.norm=NULL, P.ord=NULL){
  if (x.l >= x.r) {
    cat("error: x.l >= x.r \n")
    return(NULL)
  }
  
  simdat.l <- generate.data(
    logHR.X.site=logHR.X,
    n=n, P=P, Common.P=Common.P, coef.AonB=coef.AonB, coef.XonZ=coef.XonZ,
    coef.cens=coef.cens, scale.cens=scale.cens,
    coef.event=c(x.l,coef.event[-1]), scale.event=scale.event,
    censtype=censtype, trunc=trunc,
    P.presc.topK=P.presc.topK, prescription.mode.topK=prescription.mode.topK,
    Corr.norm=Corr.norm, Quant.norm=Quant.norm, P.ord=P.ord)
  simdat.r <- generate.data(
    logHR.X.site=logHR.X,
    n=n, P=P, Common.P=Common.P, coef.XonZ=coef.XonZ,
    coef.cens=coef.cens, scale.cens=scale.cens,
    coef.event=c(x.r,coef.event[-1]), scale.event=scale.event,
    censtype=censtype, trunc=trunc,
    P.presc.topK=P.presc.topK, prescription.mode.topK=prescription.mode.topK,
    Corr.norm=Corr.norm, Quant.norm=Quant.norm, P.ord=P.ord)
  f.l <- sum(simdat.l$Y[simdat.l$X==0])/sum(simdat.l$E[simdat.l$X==0]) - control.rate
  f.r <- sum(simdat.r$Y[simdat.r$X==0])/sum(simdat.r$E[simdat.r$X==0]) - control.rate
  
  if (f.l == 0) {
    return(x.l)
  }else if (f.r == 0) {
    return(x.r)
  }else if (f.l * f.r > 0) {
    cat("error: ftn(x.l) * ftn(x.r) > 0 \n")
    return(NULL)
  }
  while ((x.r - x.l) > tol) {
    x.m <- (x.l + x.r)/2
    simdat.m <- generate.data(
      logHR.X.site=logHR.X,
      n=n, P=P, Common.P=Common.P, coef.XonZ=coef.XonZ,
      coef.cens=coef.cens, scale.cens=scale.cens,
      coef.event=c(x.m,coef.event[-1]), scale.event=scale.event,
      censtype=censtype, trunc=trunc,
      P.presc.topK=P.presc.topK, prescription.mode.topK=prescription.mode.topK,
      Corr.norm=Corr.norm, Quant.norm=Quant.norm, P.ord=P.ord)
    f.m <- sum(simdat.m$Y[simdat.m$X==0])/sum(simdat.m$E[simdat.m$X==0]) - control.rate
    
    if (f.m == 0) {
      return(x.m)
    }
    else if (f.l * f.m < 0) {
      x.r <- x.m
      f.r <- f.m
    }
    else {
      x.l <- x.m
      f.l <- f.m
    }
    print(paste0(x.l,", ",x.r))
  }
  return((x.l + x.r)/2)
}

cox.est <- function(X,Y,E,S,Z) {

  site.ind <- model.matrix(~factor(S))[,-1]
  
  cox.pooled <- as.numeric(survival::coxph.fit(as.matrix(cbind(X,Z,site.ind)), 
                                               survival::Surv(E,as.integer(Y)),
                                               strata=NULL,control=cox.ctrl,
                                               method=tie.method,
                                               rownames=NULL)[[c('coefficients')]])
  
  cox.site <- foreach(i=sort(unique(S))) %do% {
    site.curr <- S==i
    cov.mat.site <- as.matrix(cbind(X,Z) + 0.0)[site.curr,] 
    tryCatch(
      {
        
        as.numeric(survival::coxph.fit(cov.mat.site,
                                       survival::Surv(E[site.curr],Y[site.curr]),
                                       strata=NULL,control=cox.ctrl,
                                       method=tie.method,rownames=NULL)[[c('coefficients')]])
      }, error = function(e) NULL)
  }
  list(Pooled=cox.pooled,Site.Specific=cox.site)
}


# Matthew: similar to cox.est, for binary case
logistic.est <- function(X, Y, S, Z) {
  
  site.ind <- model.matrix(~factor(S))[,-1]
  
  logistic.pooled <- as.numeric(glm(Y ~ as.matrix(cbind(X, Z, site.ind)), family = binomial)$coefficients)
  
  logistic.site <- foreach(i = sort(unique(S))) %do% {
    site.curr <- S == i
    cov.mat.site <- as.matrix(cbind(X, Z) + 0.0)[site.curr, ]
    
    tryCatch(
      {
        as.numeric(glm(Y[site.curr] ~ cov.mat.site, family = binomial)$coefficients)
      }, error = function(e) NULL
    )
  }
  
  list(Pooled = logistic.pooled, Site.Specific = logistic.site)
}

# Code and documentation for the simulation setup described in
# "Plasmode simulation for the evaluation of pharmacoepidemiologic
# methods in complex healthcare databases"
# Software: Jessica M Franklin
hdSimSetup <- function(x,treatVar,
                       smod,cmod, effectRR = 1, MM = 1,
                       size = nrow(x), eventRate = NULL,trunc=NULL) {
  require(survival)
  # x = datset on which sims are based
  # idVar = name of id variable
  # outcomeVar = name of outcome variable
  # timeVar = name of the follow-up time variable
  # treatVar = name of treatment variable
  # form = RHS of formula used for outcome simulation - should look like
  # "~ C1 + C2 + .". Can include anything allowed by coxph.
  # effectRR = the desired treatment effect relative risk
  # MM = multiplier of confounder effects on outcome on
  # the log-scale
  # nsim = number of desired outcome vectors
  # size = desired size of simulated cohort studies (i.e., # of individuals)
  # eventRate = desired average event rate -- default is the event
  # rate observed in the base dataset
  
  n <- nrow(x)
  
  sidx <- sapply(c(treatVar),
                 function(v) which(names(x) == v))
  names(x)[sidx] <- c("ZX")
  # y1 <- Surv(x$TIME, x$OUTCOME)
  # y2 <- Surv(x$TIME, !x$OUTCOME)
  # form1 <- as.formula(paste("y1 ~", form))
  # form2 <- as.formula(paste("y2 ~", form))
  # 
  # # estimate survival and censoring models
  # smod <- coxph(form1, x = TRUE, data = x)
  fit <- survfit(smod,colMeans(x))
  s0 <- fit$surv # survival curve for average patient
  ts <- fit$time
  nts <- length(ts)
  # cmod <- coxph(form2, data = x)
  fit <- survfit(cmod,colMeans(x))
  c0 <- fit$surv # censoring curve for average patient
  # find event rate in base cohort (if everyone was followed to end of study)
  Xb <- as.vector(as.matrix(x) %*% coef(smod))
  mx <- colMeans(x)
  xb0 <- mx %*% coef(smod)
  s0end <- min(s0)
  if(is.null(eventRate)) eventRate <- 1-mean(s0end^exp(Xb - xb0))
  # find delta value needed to get approximate desired event rate under new
  # parameters

  bnew <- replace(MM*coef(smod), names(coef(smod)) == "ZX", log(effectRR))
  Xbnew <- as.vector(as.matrix(x) %*% bnew)
  # sXend <- s0end^(exp(Xb - xb0))
  # fn <- function(d) mean(sXend^d) - (1 - eventRate)
  # delta <- uniroot(fn, lower = 0, upper = 20)$root
  # # setup n X nts matrix of individual survival and censoring curves under new
  # # parameters
  # exponent <- delta*exp(Xbnew - xb0)
  # Sx <- do.call(cbind,lapply(s0, function(s) s^exponent))
  # Xbnew <- as.vector(smod$x %*% coef(cmod))
  # xb0 <- mx %*% coef(cmod)
  delta <- 1
  exponent <- delta*exp(Xbnew - xb0)
  Sx <- do.call(cbind,lapply(s0, function(s) s^exponent))
  Xbnew <- as.vector(as.matrix(x) %*% coef(cmod))
  xb0 <- mx %*% coef(cmod)
  exponent <- delta*exp(Xbnew - xb0)
  Cx <- do.call(cbind,lapply(c0, function(s) s^exponent))
  
  #### sample and simulate
  
  # event time
  u <- runif(size, 0, 1)
  # the first time survival drops below u
  w <- apply(Sx < u, 1, function(x) which(x)[1])
  stime <- ts[w]
  # for any individuals with survival that never drops below u,
  # replace with arbitrary time beyond last observed event/censoring time
  w <- Sx[,nts] > u
  stime[w] <- max(ts) + 1
  
  # censoring time
  u <- runif(size, 0, 1)
  # the first time censor-free survival drops below u
  w <- apply(Cx < u, 1, function(x) which(x)[1])
  ctime <- ts[w]
  # for any individuals with censor-free survival that never drops below u,
  # replace with hard censor time at last observed event/censoring time
  w <- Cx[,nts] > u
  ctime[w] <- max(ts)
  
  # put it together
  tnew <- pmin(stime, ctime)
  names(tnew) <- "E"
  tnew[tnew>=trunc] <- trunc
  ynew <- as.integer(stime == tnew)
  
  data.frame('Y'=ynew, 'E'=tnew)
}

