 #number of parallel foreach loops, Can be number of electrodes running for simultaneous studystrap sampling
samp.list <- 11:15
samps <- length(samp.list)
num.elec <- 15 #number of electrode datasets to sample from 
straps <- 500 #study replicates to include
library(doParallel)
best.comp <- 28
tune.metric <- 5
tune.fraction <- 1
seed.type <- "y"

samp.size <- "full"
tune.grid <- c(0.0000001, 0.0000005,0.000001, 0.000005, 0.00001, 0.00005, 0.0001, 
               0.0005, 0.001, 0.005, 0.01, 0.1, 0.25, 0.5, 1, 2,5, 10)
output.filename <- paste0("Merged LASSO 0 Trunc_seed_", 
                          seed.type, "_samps_ ", min(samp.list),":",max(samp.list), "_ncomp_", best.comp, "_tune.met_", tune.metric, 
                          "_tune.frac_", tune.fraction, "_sampSize_", samp.size)
set.seed(1)
print(output.filename)
source("Study Strap Functions.R")

clust <- makeCluster(samps, outfile="")
registerDoParallel(clust)

final.results <- matrix(NA, ncol = 2, nrow = num.elec)
colnames(final.results) <- c("Elec", "RMSE")

results <- foreach(y=samp.list, .combine=list, .multicombine=TRUE) %dopar%{
  
  return.results <- function(y){
    library(glmnet) #library for model fitting
    set.seed(y) # set seed
    num.elec <- 15 # total number of electrodes

    
    results.vec <- c(y, NA)
    
    # Load Data
    full <- read.csv("combined")

    test.elec <- full[full$Electrode == y, (ncol(full) - 1000):ncol(full)] # set test electrode
    full <- full[full$Electrode != y, ] # take full dataset out of test electrode
 
    #########
    # CV
    #########
    # make CV indx sets
    train.elecs <- unique(full$Electrode)
    indx <- list(length = length(train.elecs))
    
    test.indx <- list(length = length(train.elecs))
    for ( x in 1:length(train.elecs) ){
      test.indx[[x]] <- which(full$Electrode == train.elecs[x])
      rows.indx <- which(full$Electrode != train.elecs[x])
      indx[[x]] <- sample(rows.indx, length(test.indx[[x]]) * tune.fraction, replace = FALSE)
    }
    rm(rows.indx)
    
    # Derivative
    colnames(test.elec) <- c("DA", paste0("V",1:(ncol(test.elec) - 1)))
    
    # Derivative
    full <- full[, (ncol(full) - 1000):ncol(full)] # set test electrode
    colnames(full) <- c("DA", paste0("V",1:(ncol(full) - 1)))
    
    
    ############
    # Tuning
    ############
    if (tune.metric == 5){
      print("tuning 5")
      # keep CV data
      MSE.mat <- matrix(NA, ncol = length(tune.grid), nrow = length(train.elecs)) # matrix of test errors associated with only metric metric used
      colnames(MSE.mat) <- c(paste0("MSE-Parameter: ", tune.grid))
      
      # iterate through tuning parameters

          #iterate through electrodes
          for ( x in 1:length(train.elecs) ){
            # train classifier on one electrode in training elecs

            tune.fit <- glmnet(y = as.matrix(full[indx[[x]], 1]), x = as.matrix(full[indx[[x]], -1]), 
                               lambda = tune.grid, family = "gaussian", alpha = 1)

            preds <- predict(tune.fit, as.matrix(full[test.indx[[x]], -1]))
            
            #iterate through lambda values and save RMSEs
            for ( i in 1:length(tune.grid)){
                MSE.mat[x, i] <- sqrt(mean((full$DA[ test.indx[[x]] ] - preds[,i])^2)) # save MSEs of each held out electrode for each tuning parameter
            }
            rm(preds)
                        
          }
      #}
      print(MSE.mat)
      
      MSE.mean <- colMeans(MSE.mat)
      print(MSE.mean)
      best.comp <- tune.grid[which.min(MSE.mean)] # choose best tuning parameter
      print(best.comp)
      rm(tune.fit, MSE.mat)
    
    }
    print("primary model")
    
    pcr.fit <- glmnet(y = as.matrix(full[, 1]), x = as.matrix(full[, -1]), 
                       lambda = best.comp, family = "gaussian", alpha = 1)
    
    preds <- predict(pcr.fit, as.matrix(test.elec[, -1])) #make predictions
    preds <- preds * I(preds > 0) # truncate
    results.vec[2] <- sqrt(mean((test.elec[,1] - preds)^2))
    print(paste0("LASSO Merged 0 Truncate Error: Elec_", y, ": ", results.vec))
    print(paste0("LASSO Merged 0 Truncate Lambda: Elec_", y, ": ", best.comp))
    return(results.vec)
  }
  
  return.results(y)
}


for(i in 1:samps){
  final.results[i,] <- results[[i]]
}

write.csv(final.results, output.filename)

