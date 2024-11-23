# Hidden functions to get summary statistics ----
.ordtonorm <- function (probs, Cor) {
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
#####


# Functions to get summary statistics ----
#' Title
#'
#' @param E a numeric vector of event times
#' @param Y a binary vector of outcomes
#' @param X a binary vector of exposures
#' @param B a matrix of binary covariates, with each column representing a covariate
#' @param A a matrix of categorical covariates, with each column representing a covariate
#' @param prescription.mode a numeric vector of prescription modes for exposures
#' @param my.presc.K a numeric value of the number of prescription modes to be used
#' @param tie.method specifies the method for handling ties in the Cox proportional hazards model, with options including 'efron', 'breslow', and 'exact'
#' @param method specify the type of summary statistics abstracted corresponding to the data generation method,
#' with options including 'all', 1, 2, 3. 'all' will return all summary statistics can be used for generate data with three methods,
#' while 1, 2, 3 corresponding to the multivariate normal thresholding for A and B, the multivariate normal thresholding only for B,
#' and the chain of regression methods, respectively. See more details in the manuscript.
#' @param censtype specify the type of censoring, with options including 'simple', 'simplebump', 'cov', and 'covbump',
#' see more details in the manuscript
#'
#' @return a list of summary statistics
#' @export
#'
#' @examples
#' data(example_data)
#' E <- example_data$E
#' Y <- example_data$Y
#' X <- example_data$X
#' B <- as.matrix(example_data[, c("B.1", "B.2", "B.3", "B.4")])
#' A <- as.matrix(example_data[, c("A1", "A2")])
#' summstat.survival <- get.summstat.survival(E, Y, X, B, A, method = "all")
#' summstat.survival
get.summstat.survival <- function(E,Y,X,B,A,prescription.mode=seq(30,365,by=30),
                         my.presc.K=1,tie.method="efron",method="all",censtype="simple")
  {
  # # test
  # A=C

  if (is.null(E)){
    stop("E is required for survival data")
  }

  glm.ctrl <- glm.control(epsilon = 1e-8, maxit = 25, trace = FALSE)
  cox.ctrl <- survival::coxph.control(eps=1e-09, toler.chol=.Machine$double.eps^0.75,
                                      iter.max=20, toler.inf=sqrt(1e-09), outer.max=10)

  ### (Correlated) binary covariates B
  B <- as.matrix(B)
  n.B <- ncol(B)

  #### Marginal mean, i.e. prevalence of each binary covariate
  P <- apply(as.matrix(B),2,function(XX){as.numeric(table(XX)[[2]]/sum(table(XX)))})
  #### Common probability, i.e. number of one's shared by each pair of binary variables divided by sample size
  Common.P <- t(as.matrix(B))%*%as.matrix(B)/nrow(B)
  Corr.B <- cor(B)

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
  # 9.24: catlabs may be redundant
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
  coef.chain <- vector("list", n.B)
  names(coef.chain) <- colnames(B)
  coef.chain[[1]] <- mean(B[,1])
  for (i in 2:n.B) {
    coef.chain[[i]] <- tryCatch({
      coef(glm(B[,i] ~ B[, 1:(i-1)], family = "binomial"))
    }, error = function(e) {
      message("Chain regression failed for B[", i, "]: ", e$message)
      NA
    })
  }

  #### Coef for relationship between binary and categorical confounders

  coef.AonB <- vector("list", n.A)
  names(coef.AonB) <- colnames(A)

  # For the first categorical variable, regress on all B variables
  coef.AonB[[1]] <- tryCatch({
    coef(nnet::multinom(A[, 1] ~ B, trace = FALSE))
  }, error = function(e) {
    message("Model fitting failed for A[1]: ", e$message)
    NA
  })

  # For subsequent categorical variables, regress on the indicators of previous categorical variables and B
  if (n.A > 1) {
    for (i in 2:n.A) {
      # Get the corresponding indicator columns for all previous categorical variables
      start_index <- 1
      end_index <- sum(sapply(1:(i - 1), function(k) length(unique(A[, k])) - 1))

      # Extract the relevant indicator columns
      previous_A_indicators <- A.indicator[, start_index:end_index, drop = FALSE]

      # Create a formula using previous indicators and B as predictors
      formula <- as.formula(paste("A[, i] ~ . + B"))  # Use . to represent all predictors

      # Fit the multinomial logistic regression model, catch errors due to mismatched rows or other issues
      coef.AonB[[i]] <- tryCatch({
        model_data <- data.frame(A = A[, i], previous_A_indicators, B)
        model <- nnet::multinom(A ~ ., data = model_data, trace = FALSE)
        coef(model)
      }, warning = function(w) {
        message("Warning during model fitting for A[", i, "]: ", w$message)
        NA
      }, error = function(e) {
        message("Model fitting failed for A[", i, "]: ", e$message)
        NA
      })
    }
  }


  ####################################

  ### Exposure variable X|B,A (Propensity Score Model Coefficients)
  #### Intercept and coefficients from logistic regression
  ps.fit <- glm.fit(cbind(1,B,A.indicator), X, control=glm.ctrl, family=binomial())
  class(ps.fit) <- 'glm'
  coef.XonZ <- coef(ps.fit)
  ps.by.x <- cbind(X, fitted(ps.fit))
  ps.vcov <- vcov(ps.fit)
  probs <- predict(ps.fit, type = "response")
  ps.c <- mean(sample(probs[X == 1L], 1000000L, TRUE) > sample(probs[X == 0L], 1000000L, TRUE))


  ####  Censoring distribution -- Weibull shape and scale parameters estimated through parametric survival regression:
  #### Note: use only censoring data to estimate proportion among the censored.
  #### Unique censoring distribution shape
  #### Note: If one decide to add extra distribution properties, (1) and (2) should fit model using data eliminating the unique observations,
  #### e.g. observations with E <- the top two most frequent presc pattern bumps

  ## Proportion with a prescription bump amongst those censored
  P.presc <- NULL
  # P.presc <- sapply(prescription.mode, function(pm) mean(E[Y==0] == pm)) # 9.24: changed below

  for( i in 1:length(prescription.mode))
  {
    P.presc<-cbind(P.presc,mean(E[Y==0]==prescription.mode[i]))
  }

  P.presc <- data.frame(P.presc)
  names(P.presc) <- as.character(prescription.mode)

  #### CENSORING DISTRIBUTION ####

  # Covariate matrix
  Z <- cbind(X, B, A.indicator)

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
  # 9.24: add the cox.ctrl
  cox.adjusted <- survival::coxph(survival::Surv(E, Y) ~ . , data=as.data.frame(model.matrix(~Z)[,-1]),
                        method=tie.method,x=TRUE,control=cox.ctrl)
  # class(cox.adjusted) <- "coxph"
  cox.coef.adjusted <- coef(cox.adjusted)
  names(cox.coef.adjusted) <- colnames(Z)
  cox.vcov <- vcov(cox.adjusted)

  ## Cox censoring model for plasmode simulation
  cox.adjusted.cens <-  survival::coxph(survival::Surv(E, 1-Y) ~ . , data=as.data.frame(model.matrix(~Z)[,-1]),
                              method=tie.method,x=TRUE,control=cox.ctrl)
  # class(cox.adjusted.cens) <- "coxph"

  #### Additional information
  # control.rate <- sum(Y[X==0])/sum(E[X==0])
  # compare.rate <- sum(Y[X==1])/sum(E[X==1])
  # control.events <- sum(Y[X==0])
  # compare.events <- sum(Y[X==1])
  # N.X <- table(X)
  # P.time <- c(sum(E[X==0]), sum(E[X==1]))/table(X)

  # summary statistics for the ord method
  norm.spec <- .ordtonorm(probs=P.ord, Cor=Corr.ord)

  # for the
  norm.spec.B <- .ordtonorm(probs=lapply(P,FUN=function(x){c(1-x,x)}), Cor=Corr.B)

  # return specific variables for different censor types/ methods
  output_vars = switch(as.character(method),
                       "all" = switch(censtype,
                                      "simple" = list(coef.XonZ=coef.XonZ, # method 1+2+3
                                                      Coef.cat=coef.AonB, # method 2+3
                                                      P.ord=P.ord,Quants.norm = norm.spec$quants.norm,Corr.norm = norm.spec$corr.norm, # method 1
                                                      Quants.norm.B = norm.spec.B$quants.norm,Corr.norm.B = norm.spec.B$corr.norm, # method 2
                                                      Coef.bin=coef.chain, # method 3
                                                      simple.coef.cens=simple.coef.cens, simple.scale.cens=simple.scale.cens),
                                      "simplebump" = list(coef.XonZ=coef.XonZ, # method 1+2+3
                                                          Coef.cat=coef.AonB, # method 2+3
                                                          P.ord=P.ord,Quants.norm = norm.spec$quants.norm,Corr.norm = norm.spec$corr.norm, # method 1
                                                          Quants.norm.B = norm.spec.B$quants.norm,Corr.norm.B = norm.spec.B$corr.norm, # method 2
                                                          Coef.bin=coef.chain, # method 3
                                                          simplebump.coef.cens=simplebump.coef.cens, simplebump.scale.cens=simplebump.scale.cens),
                                      "cov" = list(coef.XonZ=coef.XonZ, # method 1+2+3
                                                   Coef.cat=coef.AonB, # method 2+3
                                                   P.ord=P.ord,Quants.norm = norm.spec$quants.norm,Corr.norm = norm.spec$corr.norm, # method 1
                                                   Quants.norm.B = norm.spec.B$quants.norm,Corr.norm.B = norm.spec.B$corr.norm, # method 2
                                                   Coef.bin=coef.chain, # method 3
                                                   cov.coef.cens=cov.coef.cens, cov.scale.cens=cov.scale.cens),
                                      "covbump" = list(coef.XonZ=coef.XonZ, # method 1+2+3
                                                       Coef.cat=coef.AonB, # method 2+3
                                                       P.ord=P.ord,Quants.norm = norm.spec$quants.norm,Corr.norm = norm.spec$corr.norm, # method 1
                                                       Quants.norm.B = norm.spec.B$quants.norm,Corr.norm.B = norm.spec.B$corr.norm, # method 2
                                                       Coef.bin=coef.chain, # method 3
                                                       covbump.coef.cens=covbump.coef.cens, covbump.scale.cens=covbump.scale.cens)),
                       "1" = switch(censtype,
                                    "simple" = list(coef.XonZ=coef.XonZ, # method 1+2+3
                                                    P.ord=P.ord,Quants.norm = norm.spec$quants.norm,Corr.norm = norm.spec$corr.norm, # method 1
                                                    simple.coef.cens=simple.coef.cens, simple.scale.cens=simple.scale.cens),
                                    "simplebump" = list(coef.XonZ=coef.XonZ, # method 1+2+3
                                                        P.ord=P.ord,Quants.norm = norm.spec$quants.norm,Corr.norm = norm.spec$corr.norm, # method 1
                                                        simplebump.coef.cens=simplebump.coef.cens, simplebump.scale.cens=simplebump.scale.cens),
                                    "cov" = list(coef.XonZ=coef.XonZ, # method 1+2+3
                                                 P.ord=P.ord,Quants.norm = norm.spec$quants.norm,Corr.norm = norm.spec$corr.norm, # method 1
                                                 cov.coef.cens=cov.coef.cens, cov.scale.cens=cov.scale.cens),
                                    "covbump" = list(coef.XonZ=coef.XonZ, # method 1+2+3
                                                     P.ord=P.ord,Quants.norm = norm.spec$quants.norm,Corr.norm = norm.spec$corr.norm, # method 1
                                                     covbump.coef.cens=covbump.coef.cens, covbump.scale.cens=covbump.scale.cens)),
                       "2" = switch(censtype,
                                    "simple" = list(coef.XonZ=coef.XonZ, # method 1+2+3
                                                    Coef.cat=coef.AonB, # method 2+3
                                                    Quants.norm.B = norm.spec.B$quants.norm,Corr.norm.B = norm.spec.B$corr.norm,# method 2
                                                    simple.coef.cens=simple.coef.cens, simple.scale.cens=simple.scale.cens),
                                    "simplebump" = list(coef.XonZ=coef.XonZ, # method 1+2+3
                                                        Coef.cat=coef.AonB, # method 2+3
                                                        Quants.norm.B = norm.spec.B$quants.norm,Corr.norm.B = norm.spec.B$corr.norm,# method 2
                                                        simplebump.coef.cens=simplebump.coef.cens, simplebump.scale.cens=simplebump.scale.cens),
                                    "cov" = list(coef.XonZ=coef.XonZ, # method 1+2+3
                                                 Coef.cat=coef.AonB, # method 2+3
                                                 Quants.norm.B = norm.spec.B$quants.norm,Corr.norm.B = norm.spec.B$corr.norm,# method 2
                                                 cov.coef.cens=cov.coef.cens, cov.scale.cens=cov.scale.cens),
                                    "covbump" = list(coef.XonZ=coef.XonZ, # method 1+2+3
                                                     Coef.cat=coef.AonB, # method 2+3
                                                     Quants.norm.B = norm.spec.B$quants.norm,Corr.norm.B = norm.spec.B$corr.norm,# method 2
                                                     covbump.coef.cens=covbump.coef.cens, covbump.scale.cens=covbump.scale.cens)),
                       "3" = switch(censtype,
                                    "simple" = list(coef.XonZ=coef.XonZ, # method 1+2+3
                                                    Coef.cat=coef.AonB, # method 2+3
                                                    Coef.bin=coef.chain, # method 3
                                                    simple.coef.cens=simple.coef.cens, simple.scale.cens=simple.scale.cens),
                                    "simplebump" = list(coef.XonZ=coef.XonZ, # method 1+2+3
                                                        Coef.cat=coef.AonB, # method 2+3
                                                        Coef.bin=coef.chain, # method 3
                                                        simplebump.coef.cens=simplebump.coef.cens, simplebump.scale.cens=simplebump.scale.cens),
                                    "cov" = list(coef.XonZ=coef.XonZ, # method 1+2+3
                                                 Coef.cat=coef.AonB, # method 2+3
                                                 Coef.bin=coef.chain, # method 3
                                                 cov.coef.cens=cov.coef.cens, cov.scale.cens=cov.scale.cens),
                                    "covbump" = list(coef.XonZ=coef.XonZ, # method 1+2+3
                                                     Coef.cat=coef.AonB, # method 2+3
                                                     Coef.bin=coef.chain, # method 3
                                                     covbump.coef.cens=covbump.coef.cens, covbump.scale.cens=covbump.scale.cens)))


  return (c(list(n=length(Y),
            P.presc.topK=P.presc.topK,
            prescription.mode.topK=prescription.mode.topK,
            adj.coef.event=adj.coef.event,adj.scale.event=adj.scale.event,
            logHR.X=cox.coef.adjusted[[1]]),output_vars))
}


#' Title
#'
#' @param Y a binary vector of outcomes
#' @param X a binary vector of exposures
#' @param B a matrix of binary covariates, with each column representing a covariate
#' @param A a matrix of categorical covariates, with each column representing a covariate
#' @param method specify the type of summary statistics abstracted corresponding to the data generation method,
#' with options including 'all', 1, 2, 3. 'all' will return all summary statistics can be used for generate data with three methods,
#' while 1, 2, 3 corresponding to the multivariate normal thresholding for A and B, the multivariate normal thresholding only for B,
#' and the chain of regression methods, respectively. See more details in the manuscript.
#'
#' @return a list of summary statistics
#' @export
#'
#' @examples
#' data(example_data)
#' E <- example_data$E
#' Y <- example_data$Y
#' X <- example_data$X
#' B <- as.matrix(example_data[, c("B.1", "B.2", "B.3", "B.4")])
#' A <- as.matrix(example_data[, c("A1", "A2")])
#' summstat.binary <- get.summstat.binary(Y, X, B, A, method = "all")
#' summstat.binary
get.summstat.binary <- function(Y,X,B,A,method="all"){
  # #test
  # A=C
  # method=1

  glm.ctrl <- glm.control(epsilon = 1e-8, maxit = 25, trace = FALSE)

  ### (Correlated) binary covariates B
  B <- as.matrix(B)
  n.B <- ncol(B)

  #### Marginal mean, i.e. prevalence of each binary covariate
  P <- apply(as.matrix(B),2,function(XX){as.numeric(table(XX)[[2]]/sum(table(XX)))})
  #### Common probability, i.e. number of one's shared by each pair of binary variables divided by sample size
  Common.P <- t(as.matrix(B))%*%as.matrix(B)/nrow(B)
  Corr.B <- cor(B)

  ## Categorical variables A
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
  coef.chain <- vector("list", n.B)
  names(coef.chain) <- colnames(B)
  coef.chain[[1]] <- mean(B[,1])
  for (i in 2:n.B) {
    coef.chain[[i]] <- tryCatch({
      coef(glm(B[,i] ~ B[, 1:(i-1)], family = "binomial"))
    }, error = function(e) {
      message("Chain regression failed for B[", i, "]: ", e$message)
      NA
    })
  }

  #### Coef for relationship between binary and categorical confounders

  coef.AonB <- vector("list", n.A)
  names(coef.AonB) <- colnames(A)

  # For the first categorical variable, regress on all B variables
  coef.AonB[[1]] <- tryCatch({
    coef(nnet::multinom(A[, 1] ~ B, trace = FALSE))
  }, error = function(e) {
    message("Model fitting failed for A[1]: ", e$message)
    NA
  })

  # For subsequent categorical variables, regress on the indicators of previous categorical variables and B
  if (n.A > 1) {
    for (i in 2:n.A) {
      # Get the corresponding indicator columns for all previous categorical variables
      start_index <- 1
      end_index <- sum(sapply(1:(i - 1), function(k) length(unique(A[, k])) - 1))

      # Extract the relevant indicator columns
      previous_A_indicators <- A.indicator[, start_index:end_index, drop = FALSE]

      # Create a formula using previous indicators and B as predictors
      formula <- as.formula(paste("A[, i] ~ . + B"))  # Use . to represent all predictors

      # Fit the multinomial logistic regression model, catch errors due to mismatched rows or other issues
      coef.AonB[[i]] <- tryCatch({
        model_data <- data.frame(A = A[, i], previous_A_indicators, B)
        model <- nnet::multinom(A ~ ., data = model_data, trace = FALSE)
        coef(model)
      }, warning = function(w) {
        message("Warning during model fitting for A[", i, "]: ", w$message)
        NA
      }, error = function(e) {
        message("Model fitting failed for A[", i, "]: ", e$message)
        NA
      })
    }
  }


  ####################################

  ### Exposure variable X|B,A (Propensity Score Model Coefficients)
  #### Intercept and coefficients from logistic regression
  ps.fit <- glm.fit(cbind(1,B,A.indicator), X, control=glm.ctrl, family=binomial())
  class(ps.fit) <- 'glm'
  coef.XonZ <- coef(ps.fit)
  ps.by.x <- cbind(X, fitted(ps.fit))
  ps.vcov <- vcov(ps.fit)
  probs <- predict(ps.fit, type = "response")
  ps.c <- mean(sample(probs[X == 1L], 1000000L, TRUE) > sample(probs[X == 0L], 1000000L, TRUE))

  #### Matthew added: Binary outcome
  #### Logistic regression for the event distribution

  logit.fit <- glm.fit(cbind(1,X,B,A.indicator), Y, control=glm.ctrl, family=binomial())
  class(logit.fit) <- 'glm'
  coef.Yon1 <- coef(logit.fit)[1]
  coef.YonX <- coef(logit.fit)[2]
  coef.YonZ <- coef(logit.fit)[-c(1,2)]

  control.events <- sum(Y[X==0])
  compare.events <- sum(Y[X==1])
  N.X <- table(X)

  # summary statistics for the ord method
  norm.spec <- .ordtonorm(probs=P.ord, Cor=Corr.ord)

  norm.spec.B <- .ordtonorm(probs=lapply(P,FUN=function(x){c(1-x,x)}), Cor=Corr.B)


  # return binary outcome version
  if (method=="all"){
    return(list(n=length(Y),
                coef.Yon1=coef.Yon1,coef.YonX=coef.YonX,coef.YonZ=coef.YonZ,
                coef.XonZ=coef.XonZ, # method 1+2+3
                Coef.cat=coef.AonB, # method 2+3
                P.ord=P.ord,Quants.norm = norm.spec$quants.norm,Corr.norm = norm.spec$corr.norm, # method 1
                Quants.norm.B = norm.spec.B$quants.norm,Corr.norm.B = norm.spec.B$corr.norm, # method 2
                Coef.bin=coef.chain # method 3
                )
           )
  } else if (method==1){
    return(list(n=length(Y),
                coef.Yon1=coef.Yon1,coef.YonX=coef.YonX,coef.YonZ=coef.YonZ,
                coef.XonZ=coef.XonZ, # method 1+2+3
                P.ord=P.ord,Quants.norm = norm.spec$quants.norm,Corr.norm = norm.spec$corr.norm # method 1
                )
           )
  } else if (method==2){
    return(list(n=length(Y),
                coef.Yon1=coef.Yon1,coef.YonX=coef.YonX,coef.YonZ=coef.YonZ,
                coef.XonZ=coef.XonZ, # method 1+2+3
                Coef.cat=coef.AonB, # method 2+3
                Quants.norm.B = norm.spec.B$quants.norm,Corr.norm.B = norm.spec.B$corr.norm # method 2
                )
           )
  } else if (method==3){
    return(list(n=length(Y),
                coef.Yon1=coef.Yon1,coef.YonX=coef.YonX,coef.YonZ=coef.YonZ,
                coef.XonZ=coef.XonZ, # method 1+2+3
                Coef.cat=coef.AonB, # method 2+3
                Coef.bin=coef.chain # method 3
                )
           )
  } else { stop("Invalid method specified")}
}
#####


# Hidden functions used to generate data ----
.ordgendata <- function(n, sigma, quants.norm){
  retval = mvtnorm::rmvnorm(n = n, sigma = sigma)
  for (i in 1:ncol(sigma)) {
    retval[, i] = cut(x = retval[, i], breaks = c(-1/0, quants.norm[[i]]),
                      right = FALSE)
  }
  retval - 1
}


.gencov.ord <- function(n, P.ord, Quant.norm, Corr.norm, coef.XonZ){
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


.gencov <- function(n, Corr.norm.B, Quant.norm.B, coef.AonB, coef.XonZ){
  B <- .ordgendata(n, sigma=Corr.norm.B, quants.norm=Quant.norm.B)
  # B <- .rmvbin(n, margprob=P, commonprob=Common.P) # TODO: replace

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
      # P.A2 <- exp( cbind(rep(1,n),B)%*%t(coef.AonB[[i]]) )
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


.gencov.chain <- function(n, coef.chain, coef.AonB, coef.XonZ,names.cat=NULL){
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
      P.A2 <- exp( cbind(rep(1,n),A.indicator,B)%*%t(coef.AonB[[i]]) ) # Matthew: Q: exp( cbind(rep(1,n),B)%*%t(coef.AonB[[i]]) ) Question: want to regress A[2] on A[1]?
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
} # TODO: intuitive, capture more info? but no theory; different property compares to other method?


.gendata.survival <- function(logHR.X.site=NULL, # optional
                              P.ord=NULL, Quant.norm=NULL, Corr.norm=NULL, # method 1 .gencov.ord
                              Corr.norm.B=NULL, Quant.norm.B=NULL, # method 2 .gencov
                              coef.chain=NULL, # method 3 .gencov.chain
                              user.data=NULL,noX=FALSE, # method 4 user data
                              n, coef.XonZ, coef.AonB=NULL, # required
                              coef.event, scale.event, coef.cens, scale.cens,
                              censtype="simple", trunc=365,
                              method=1,
                              strat.var.cens=NULL,strat.var.event=NULL,
                              P.presc.topK=NULL, prescription.mode.topK=NULL # censtype%in%c("simplebump","covbump")
                              )
  {
  # # test 10.20
  # noX=FALSE


  # if (is.null(coef.cens) || is.null(scale.cens) || is.null(coef.event) || is.null(scale.event)){
  #   stop("For survival outcome, coef.cens, scale.cens, coef.event, and scale.event must be specified")
  # }
  if (method == 4 && is.null(user.data)) stop("User data must be provided for method 4")

  ### set up specified HR of X
  if (!noX && !is.null(logHR.X.site)) {
    coef.event[grep("X",names(coef.event))] <- -logHR.X.site*scale.event
  }

  ### Generate covariates via specified method
  exppluscovs <- switch(method,
                        `1` = .gencov.ord(n, P.ord, Quant.norm, Corr.norm, coef.XonZ),
                        `2` = .gencov(n, Corr.norm.B, Quant.norm.B, coef.AonB, coef.XonZ),
                        `3` = .gencov.chain(n, coef.chain, coef.AonB, coef.XonZ),
                        `4` = {
                          exppluscovs <- list(B=NULL, A.indicator=NULL)
                          exppluscovs$A.indicator <- user.data[,-1]
                          exppluscovs$X <- if (noX) NULL else user.data[, 1]
                          exppluscovs
                        },
                        stop("Invalid method specified"))

  ### Confounders Z
  Z.model.data <- cbind(exppluscovs$B, exppluscovs$A.indicator)

  ### Exposure variable X
  X <- exppluscovs$X

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
      censorT[strat.var.cens==0] <- ceiling((-log(u[strat.var.cens==0])*exp( # TODO: strat.var.cens missing
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

  } else
  {
    E <- apply(cbind(eventT,trunc,censorT),1,min,na.rm=T)
    Y <- as.numeric(ifelse(E == eventT,1,0))

    E.1 <- apply(cbind(eventT.1,trunc,censorT.1),1,min,na.rm=T)
    E.0 <- apply(cbind(eventT.0,trunc,censorT.0),1,min,na.rm=T)

    Y.1 <- as.numeric(ifelse(E.1 == eventT.1,1,0))
    Y.0 <- as.numeric(ifelse(E.0 == eventT.0,1,0))
  }


  ## Marginal sample needs to be same size as original
  marg.x1 <- sample(1:n,size=sum(X))

  return(list(marginal.dat=rbind(cbind(x=1, y=Y.1, obst=E.1)[marg.x1,],cbind(x=0, y=Y.0, obst=E.0)[-marg.x1,]),
              B=exppluscovs$B,A.indicator=exppluscovs$A.indicator,
              X=exppluscovs$X,Y=Y,E=E)
    )

}


.gendata.binary <- function(logOR.X=NULL, # optional
                            P.ord=NULL, Corr.norm=NULL, Quant.norm=NULL,  # method 1 .gencov.ord
                            Corr.norm.B=NULL, Quant.norm.B=NULL, # method 2 .gencov
                            coef.chain=NULL, # method 3 .gencov.chain
                            user.data=NULL,noX=FALSE, # method 4 user data
                            n, coef.AonB=NULL, coef.XonZ, # required
                            coef.Yon1, coef.YonX, coef.YonZ,
                            method=1)
{
  # # test 10.8
  # n=SS$n
  # P=SS$P
  # Common.P=SS$Common.P
  # coef.XonZ=SS$coef.XonZ
  # coef.AonB=SS$Coef.cat
  # method=2
  # coef.Yon1=SS$coef.Yon1
  # coef.YonX=SS$coef.YonX
  # coef.YonZ=SS$coef.YonZ
  # logOR.X=NULL
  # Quant.norm.B=SS$Quants.norm.B
  # Corr.norm.B=SS$Corr.norm.B


  if (is.null(coef.Yon1) || is.null(coef.YonX) || is.null(coef.YonZ)) {
    stop("For binary outcome, coef.Yon1, coef.YonX, and coef.YonZ must be specified")
  }
  if (method == 4 && is.null(user.data)) stop("User data must be provided for method 4")

  ### Set custom coef.YonX if provided
  if (!is.null(logOR.X)) {
    coef.YonX <- logOR.X
  }

  ### Generate covariates via specified method

  exppluscovs <- switch(method, # 9.24: try this new code
                        `1` = .gencov.ord(n, P.ord, Quant.norm, Corr.norm, coef.XonZ),
                        `2` = .gencov(n, Corr.norm.B, Quant.norm.B, coef.AonB, coef.XonZ),
                        `3` = .gencov.chain(n, coef.chain, coef.AonB, coef.XonZ),
                        `4` = {
                          exppluscovs <- list(B=NULL,A.indicator=NULL)
                          exppluscovs$A.indicator <- user.data[,-c(1)]
                          exppluscovs$X <- if(noX) NULL else user.data[,1]
                          exppluscovs
                        },
                        stop("Invalid method specified")
  )

  ### Confounders Z
  Z.model.data <- cbind(exppluscovs$B, exppluscovs$A.indicator)

  ### Exposure variable X
  X <- exppluscovs$X

  #### Matthew added: Binary outcome
  Y <- rbinom(n,size=1,
              p=1 / (1 + exp(-(coef.Yon1 + X * coef.YonX + Z.model.data %*% coef.YonZ))))
  Y.1 <- rbinom(n,size=1,
                p=1 / (1 + exp(-(coef.Yon1 + coef.YonX + Z.model.data %*% coef.YonZ))))
  Y.0 <- rbinom(n,size=1,
                p=1 / (1 + exp(-(coef.Yon1 + Z.model.data %*% coef.YonZ))))

  ## Marginal sample needs to be same size as original
  marg.x1 <- sample(1:n,size=sum(X))

  return(list(marginal.dat=rbind(cbind(x=1, y=Y.1)[marg.x1,],cbind(x=0, y=Y.0)[-marg.x1,]),
              B=exppluscovs$B,A.indicator=exppluscovs$A.indicator,
              X=exppluscovs$X,Y=Y)
  )
}
#####


# Functions for generating data ----
#' Title
#'
#' @param Summ.Stat a nested list where each element is a list of summary statistics for a specific data site, formatted consistently with the output from the get.summstat.survival function
#' @param n the number of samples to generate, the default is using the sample size from the summary statistics
#' @param censtype specify the type of censoring, with options including 'simple', 'simplebump', 'cov', and 'covbump',
#' see more details in the manuscript
#' @param trunc the truncation time for the survival data
#' @param method specify the type of data generation method, which consistent with the summary statistics input.
#' Options including 1, 2, 3, corresponding to the multivariate normal thresholding for A and B, the multivariate normal thresholding only for B,
#' and the chain of regression methods, respectively. See more details in the manuscript.
#' @param set.logHR.X user can manually specify the log hazard ratio for the exposure variable X. If NULL, it use the log hazard ratio from the summary statistics.
#'
#' @return a dataframe with the generated survival data, where the columns are B, indicator of A, X, Y, E, and site
#' @export
#'
#' @examples
#' data(summstat.survival)
#' summstat.survival = list(summstat.survival)
#' simulated_data_survival <- generate.data.survival(Summ.Stat=summstat.survival,method=1)
#' head(simulated_data_survival)
generate.data.survival <- function(Summ.Stat,n=NULL,censtype="simple", trunc=365,method=1, set.logHR.X=NULL){
  # # test 10.20
  # censtype="simple"
  # trunc=365
  # method=2
  # set.logHR.X=NULL
  # i=1

  n.sites<-length(Summ.Stat)

  ## Data.Site: simulated data across sites
  Data.Simulated <- NULL
  Data.Marg <- NULL

  ## Data.Pool.forPS: pooled data for est pooled PS
  Data.Pool.forPS <- NULL

  ## Loop through sites to generate site specific data
  for(i in 1:n.sites)
  {
    ## read in summary statistics SS
    SS <- Summ.Stat[[i]]
    ## get common parameters
    if(is.null(set.logHR.X)) {logHR.X.site <- SS$logHR.X} else {logHR.X.site <- set.logHR.X}
    if(is.null(n)) {n <- SS$n}
    coef.XonZ <- SS$coef.XonZ
    coef.event <- SS$adj.coef.event
    scale.event <- SS$adj.scale.event

    # set coef.cens and scale.cens for different methods
    if (censtype == "simple") {
      coef.cens <- SS$simple.coef.cens
      scale.cens <- SS$simple.scale.cens
    } else if (censtype == "simplebump") {
      coef.cens <- SS$simplebump.coef.cens
      scale.cens <- SS$simplebump.scale.cens
    } else if (censtype == "cov") {
      coef.cens <- SS$cov.coef.cens
      scale.cens <- SS$cov.scale.cens
    } else if (censtype == "covbump") {
      coef.cens <- SS$covbump.coef.cens
      scale.cens <- SS$covbump.scale.cens
    } else {
      stop("Censoring type is misspecified")
    }

    P.presc.topK <- if (censtype %in% c("simplebump", "covbump")) SS$P.presc.topK else NULL
    prescription.mode.topK <- if (censtype %in% c("simplebump", "covbump")) SS$prescription.mode.topK else NULL

    ## generate site specific data
    if (method==1){
      Corr.norm <- SS$Corr.norm
      Quant.norm <- SS$Quants.norm
      P.ord <- SS$P.ord

      DS <- .gendata.survival(
        logHR.X.site = logHR.X.site,
        n = n, coef.XonZ = coef.XonZ,
        coef.cens = coef.cens, scale.cens = scale.cens,
        coef.event = coef.event, scale.event = scale.event,
        censtype = censtype, trunc = trunc,
        P.presc.topK = P.presc.topK, prescription.mode.topK = prescription.mode.topK,
        method = method, Corr.norm = Corr.norm, Quant.norm = Quant.norm, P.ord = P.ord
      )
    } else if (method==2){
      coef.AonB <- SS$Coef.cat
      Corr.norm.B <- SS$Corr.norm.B
      Quant.norm.B <- SS$Quants.norm.B

      DS <- .gendata.survival(
        logHR.X.site = logHR.X.site,
        n = n, Corr.norm.B=Corr.norm.B, Quant.norm.B=Quant.norm.B, coef.XonZ = coef.XonZ,
        coef.cens = coef.cens, scale.cens = scale.cens,
        coef.event = coef.event, scale.event = scale.event,
        censtype = censtype, trunc = trunc, coef.AonB = coef.AonB,
        P.presc.topK = P.presc.topK, prescription.mode.topK = prescription.mode.topK,
        method = method
      )
    } else if (method==3){
      coef.chain <- SS$Coef.bin
      coef.AonB <- SS$Coef.cat

      DS <- .gendata.survival(
        logHR.X.site = logHR.X.site,
        n = n, coef.XonZ = coef.XonZ,
        coef.cens = coef.cens, scale.cens = scale.cens,
        coef.event = coef.event, scale.event = scale.event,
        censtype = censtype, trunc = trunc, coef.chain = coef.chain, coef.AonB = coef.AonB,
        P.presc.topK = P.presc.topK, prescription.mode.topK = prescription.mode.topK,
        method = method
      )
    } # else if (method==4){}

    if(is.null(DS)) stop("Censoring type is misspecified")


    ## save site specific data
    if (n.sites > 1) {
      Data.Simulated <- rbind(Data.Simulated, data.frame(B=DS$B,A=DS$A.indicator,X=DS$X,Y=DS$Y,E=DS$E,site=i))
      Data.Marg <- rbind(Data.Marg,data.frame(site=i,DS$marginal.dat))
    } else {
      Data.Simulated <- data.frame(B=DS$B, A=DS$A.indicator, X=DS$X,Y=DS$Y,E=DS$E)
      Data.Marg <- data.frame(DS$marginal.dat)
    }

  }

  # return(list(Data.Simulated=Data.Simulated,Data.Marginal=Data.Marg)) # Data.Marginal seems not useful?
  return(Data.Simulated)
}


#' Title
#'
#' @param Summ.Stat a nested list where each element is a list of summary statistics for a specific data site, formatted consistently with the output from the get.summstat.binary function
#' @param n the number of samples to generate, the default is using the sample size from the summary statistics
#' @param method specify the type of data generation method, which consistent with the summary statistics input.
#' Options including 1, 2, 3, corresponding to the multivariate normal thresholding for A and B, the multivariate normal thresholding only for B,
#' and the chain of regression methods, respectively. See more details in the manuscript.
#' @param set.logOR.X user can manually specify the log ODDs ratio for the exposure variable X. If NULL, it use the log ODDs ratio from the summary statistics.
#'
#' @return a dataframe with the generated survival data, where the columns are B, indicator of A, X, Y, and site
#' @export
#'
#' @examples
#' data(summstat.binary)
#' summstat.binary = list(summstat.binary)
#' simulated_data_binary <- generate.data.binary(Summ.Stat=summstat.binary, method=1)
#' head(simulated_data_binary)
generate.data.binary <- function(Summ.Stat,n=NULL,method=1, set.logOR.X=NULL){
  # # test 11.12
  # method=2
  # set.logOR.X=NULL
  # i=1

  n.sites<-length(Summ.Stat)

  ## Data.Site: simulated data across sites
  Data.Simulated <- NULL
  Data.Marg <- NULL

  ## Data.Pool.forPS: pooled data for est pooled PS
  Data.Pool.forPS <- NULL

  ## Loop through sites to generate site specific data
  for(i in 1:n.sites)
  {
    ## read in summary statistics SS
    SS <- Summ.Stat[[i]]
    if(is.null(n)) {n <- SS$n}
    ## generate site specific data

    # TODO: method 4 working?
    if (method==1){
      DS <- .gendata.binary(n=n, coef.XonZ=SS$coef.XonZ,
                            method=method, Corr.norm=SS$Corr.norm, Quant.norm=SS$Quants.norm, P.ord=SS$P.ord,
                            coef.Yon1=SS$coef.Yon1, coef.YonX=SS$coef.YonX, coef.YonZ=SS$coef.YonZ,
                            logOR.X=set.logOR.X)
    } else if (method==2){
      DS <- .gendata.binary(n=n, Corr.norm.B=SS$Corr.norm.B, Quant.norm.B=SS$Quants.norm.B, coef.XonZ=SS$coef.XonZ,
                            coef.AonB=SS$Coef.cat,
                            method=method,
                            coef.Yon1=SS$coef.Yon1, coef.YonX=SS$coef.YonX, coef.YonZ=SS$coef.YonZ,
                            logOR.X=set.logOR.X)
    } else if (method==3){
      DS <- .gendata.binary(n=n, coef.XonZ=SS$coef.XonZ,
                            coef.chain=SS$Coef.bin, coef.AonB=SS$Coef.cat,
                            method=method,
                            coef.Yon1=SS$coef.Yon1, coef.YonX=SS$coef.YonX, coef.YonZ=SS$coef.YonZ,
                            logOR.X=set.logOR.X)
    } # else if (method==4){}
    # DS <- .gendata.binary(n=n, P=SS$P, Common.P=SS$Common.P, coef.XonZ=SS$coef.XonZ,
    #                      coef.chain=SS$Coef.bin, coef.AonB=SS$Coef.cat,
    #                      method=method, Corr.norm=SS$Corr.norm, Quant.norm=SS$Quants.norm, P.ord=SS$P.ord,
    #                      coef.Yon1=SS$coef.Yon1, coef.YonX=SS$coef.YonX, coef.YonZ=SS$coef.YonZ,
    #                      logOR.X=set.logOR.X)

    ## save site specific data
    if (n.sites > 1) {
      Data.Simulated <- rbind(Data.Simulated, data.frame(B=DS$B,A=DS$A.indicator,X=DS$X,Y=DS$Y,site=i))
      Data.Marg <- rbind(Data.Marg,data.frame(site=i,DS$marginal.dat))
    } else {
      Data.Simulated <- data.frame(B=DS$B, A=DS$A.indicator, X=DS$X,Y=DS$Y)
      Data.Marg <- data.frame(DS$marginal.dat)
    }

  }

  # return(list(Data.Simulated=Data.Simulated,Data.Marginal=Data.Marg))
  return(Data.Simulated)
}


# TODO:
# add notations critical for the package
# high dimensional?





