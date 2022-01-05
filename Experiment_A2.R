# Scrip For Christian - Experiment A2
# Experiment A2 has strong heterogeneity (normal dist of data)


library("mvtnorm")
library(grf)
library(doParallel)
library(doSNOW)
library(doRNG)
library(FNN)

# Part 1: Data Generating Processes

# Experiment A2: Strong Treatment Effect Heterogeneity (Norm Dist)
experimenta2 <- function(n, 
                         n.test, 
                         d, 
                         prop,
                         noise)
{
  mean <- 2
  X <- matrix(rnorm(n * d, mean, 1), n, d)
  X.test <- matrix(rnorm(n.test * d, mean, 1), n.test, d)
  W <- rbinom(n, 1, prop)
  
  tau <- (2*X.test[,1]) + X.test[,2]
  Y <- ((2*X[,1]) + X[,2])* W + X[,3] + rnorm(n, 0, noise)
  
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
              differential.pred = test_cal[2,1]))
}

# Part 4: Simulation

simulation_procedure <- function(d) {
  n <- 4000
  n.test <- 1000
  noise <- 0.5
  data <- experimenta2(n,
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

# Part 5: Running Script

n.simulation <- 1000
parameter.values <- c(6,8,12,16,20)
form <- "linear"
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
# using the doRNG package to guarantee reproducible results
set.seed(1)

for(parameter in parameter.values){
  print(paste("Running ",parameter, ":"))
  results = foreach(i=1:n.simulation,
                    .combine=rbind,
                    .options.snow=opts,
                    .packages=c('grf', 'FNN')) %dopar% {
                      simulation_procedure(parameter)
                    }
  output <- rbind.data.frame(output, results)
}
output <- setNames(output, columns)
save.image(paste("Experiment_A2.RData",sep=""))