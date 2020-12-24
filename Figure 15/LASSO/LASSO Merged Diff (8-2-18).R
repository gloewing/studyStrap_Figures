 #number of parallel foreach loops, Can be number of electrodes running for simultaneous studystrap sampling
samps <- 15

num.elec <- 15 #number of electrode datasets to sample from 
straps <- 500 #study replicates to include
library(doParallel)
library(glmnet)
set.seed(1)

clust <- makeCluster(samps, outfile="")
registerDoParallel(clust)


results <- foreach(y=1:samps, .combine=list, .multicombine=TRUE) %dopar%{
  
  return.results <- function(y){
    library(glmnet) #library for model fitting
    set.seed(y) # set seed
    num.elec <- 15 # total number of electrodes

    # Load Data
    full <- as.data.frame(read.csv("combined"))
    test.elec <- full[full$Electrode == y, (ncol(full) - 1000):ncol(full)] # set test electrode
    full <- as.matrix(full[full$Electrode != y, (ncol(full) - 1000):ncol(full)]) # take full dataset out of test electrode
    
    # Derivative
    test.elec <- cbind(test.elec[,1], t(diff(t(test.elec[,-1]))))
    colnames(test.elec) <- c("DA", paste0("V",1:(ncol(test.elec) - 1)))
    
    # Derivative
    full <- cbind(full[,1], t(diff(t(full[,-1]))))
    colnames(full) <- c("DA", paste0("V",1:(ncol(full) - 1)))
    
    
    
    best.comp <- 0.25
    cvFit <- glmnet(x = as.matrix(full[,-1]),
                    y = as.matrix(full[,1]),
                    family = "gaussian", 
                    lambda = best.comp)

    preds <- predict(cvFit, as.matrix(test.elec[, -1])) #make predictions
    results.vec <- sqrt(mean((test.elec[, 1] - preds)^2))
    print(paste0("LASSO Merged Error: Elec_", y, ": ", results.vec))
  
    return(results.vec)
  }
  
  
  
  return.results(y)
}
