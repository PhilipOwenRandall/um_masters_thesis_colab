# Scrip For Christian - Experiment C2
# C1 now hides from the model randomly distributed variables, the subgroup noise and the sample identifiers
# Tau is not effected by subgroup noise.

# Required Packages - Dependencies Automated Install
#install.packages(c("grf", "doParallel", "doSNOW", "doRNG", "FNN","mvtnorm"))

# Libraries Used
library("mvtnorm")
library(grf)
library(doParallel)
library(doSNOW)
library(doRNG)
library(FNN)

# Part 1: Data Generating Processes
experimentc2 <- function(n, 
                         n.test, 
                         d, 
                         prop,
                         noise)
{
  mean <- rep(1,6)
  # use d=10 as a fixed value so that the dimensions can be used for sigma
  X <- matrix(rnorm(n * 6, mean, 1), n, 6)
  X.test = matrix(rnorm(n.test * 6, mean, 1), n.test, 6)
  W <- rbinom(n, 1, prop)
  
  # 500 for 500 subgroups, d as the sigma for the normal distribution of the function u(s)
  U <- rnorm(500, 0, d)
  U.r <- rep(U,n/500)
  S <- rep(seq(1,500,1), n/500)
  X <- cbind(X,S,U.r)
  colnames(X)[colnames(X) == 'U.r'] <- 'U'
  
  S.test <- rep(seq(1,500,1), n.test/500)
  U.test <- rep(U,n.test/500)
  
  X.test <- cbind(X.test,S.test,U.test)
  colnames(X.test)[colnames(X.test) == 'U.test'] <- 'U'
  colnames(X.test)[colnames(X.test) == 'S.test'] <- 'S'
  
  tau <- 2*X.test[,1] + X.test[,2] + U.test
  Y <- (2*X[,1] + X[,2])* W + X[,3]+ U.r + rnorm(n, 0, noise)
  
  return(list(X=X,
              X.test=X.test,
              W=W,
              Y=Y,
              tau=tau))
}

# Part 2: Performance Evaluation

mse <- function(predictions, true) {
  return(mean((predictions-true)^2))
}



bias <- function(predictions, true){
  return(abs(mean(predictions-true)))
}

covered <- function(predictions, true, sigma) {
  return(mean(abs(predictions-true) / sigma <= 1.96))
}

evaluate <- function(predictions, true, sigma) {
  return(list(mse = mse(predictions,true),
              bias = bias(predictions,true),
              coverage = covered(predictions,true,sigma)))
}

# Part 3: Estimators

# Causal Forest with user-specific parameters
# Honesty can be disabled by setting 'honest = FALSE'
# Local centering can be disabled by setting Y.hat and W.hat to vectors of zeros

CF_estimator <- function(X,
                         X.test,
                         Y,
                         W,
                         Y.hat = NULL,
                         W.hat = NULL,
                         num.trees = 1000,
                         clusters = NULL,
                         honesty = TRUE,
                         honesty.fraction = 0.5,
                         honesty.prune.leaves = TRUE,
                         min.node.size = 5) {
  CF <- causal_forest(X,
                      Y,
                      W,
                      Y.hat = Y.hat,
                      W.hat = W.hat,
                      num.trees = num.trees,
                      clusters = clusters,
                      min.node.size = min.node.size,
                      honesty = honesty,
                      honesty.fraction = honesty.fraction,
                      honesty.prune.leaves = honesty.prune.leaves,
                      seed = 1)
  
  estimates <- predict(CF,
                       X.test,
                       estimate.variance = TRUE)
  
  ate <- average_treatment_effect(CF, target.sample = "all")
  
  test_cal <- test_calibration(CF) 
  
  return(list(predictions = estimates$predictions,
              sigma = sqrt(estimates$variance.estimates),
              ate = ate,
              var.imp = variable_importance(CF),
              mean.pred = test_cal[1,1],
              differential.pred = test_cal[2,1]
  ))
}

# Part 4: Simulations
simulation_procedure <- function(d) {
  n <- 4000
  n.test <- 1000
  noise <- 0.5
  data <- experimentc2(n,          # The name of this function needs to be changed.
                       n.test,
                       d=d,
                       prop = 0.5,
                       noise)
  
  # causal forest  
  cf <- CF_estimator(data$X,
                     data$X.test,
                     data$Y,
                     data$W,
                     num.trees = 1000)
  
  cf.evaluation <- evaluate(cf$predictions, data$tau, cf$sigma)
  
  return (data.frame(n = n,
                     d = d,
                     cf.mse = cf.evaluation$mse, 
                     cf.bias = cf.evaluation$bias,
                     cf.sigma = mean(cf$sigma),
                     cf.coverage = cf.evaluation$coverage,
                     cf.mean.pred = cf$mean.pred,
                     cf.differential.pred = cf$differential.pred))
  
}

# Part 5: Running Scrip with parallelisation

n.simulation <- 1000
parameter.values <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.75, 1, 2)
columns = c("n",
            "d",
            "cf.mse",
            "cf.bias",
            "cf.sigma",
            "cf.coverage",
            "mean.pred",
            "differential.pred")


cores=detectCores()
cl <- makeCluster(cores[1])
registerDoSNOW(cl)

# progess bar
pb <- txtProgressBar(max=n.simulation, style=3)
progress <- function(n) 
{
  setTxtProgressBar(pb, n)
}
opts <- list(progress=progress)

# initialize dataframe
output <- setNames(data.frame(matrix(ncol = length(columns), nrow = 0)),
                   columns)

# run simulation
set.seed(1)

for(parameter in parameter.values){
  print(paste("Running ",parameter, ":"))
  results = foreach(i=1:n.simulation,
                    .combine=rbind,
                    .options.snow=opts,
                    .packages=c('grf', 'FNN','mvtnorm')) %dopar% {
                      simulation_procedure(parameter)
                    }
  output <- rbind.data.frame(output, results)
}
output <- setNames(output, columns)
save.image(paste("Experiment_C2_RWTH.RData",sep=""))
