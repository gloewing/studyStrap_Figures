 #number of parallel foreach loops, Can be number of electrodes running for simultaneous studystrap sampling
samp.list <- 1:15
samps <- length(samp.list)
num.elec <- 15 #number of electrode datasets to sample from 
library(doParallel)
best.comp <- 28
tune.metric <- 5
tune.fraction <- 1
seed.type <- "y"
comp.vec <- c(29, 22, 20, 29, 17, 21, 12, 30, 21, 21, 21, 21, 27, 13, 29)
samp.size <- 2500 #
tune.grid <- seq(15,45, by = 2)
output.filename <- paste0("Merged PCR Expedited Tuning 0 Trunc_seed_", 
                          seed.type, "_samps_ ", min(samp.list),":",max(samp.list), 
                          "_ncomp_", best.comp, "_tune.met_", tune.metric, 
                          "_tune.frac_", tune.fraction, "_samp_size_", samp.size)
max.comp <- max(tune.grid)
set.seed(1)
print(output.filename)
source("Study Strap Functions.R")
clust <- makeCluster(samps, outfile="")
registerDoParallel(clust)
final.results <- matrix(NA, ncol = 3, nrow = num.elec)
colnames(final.results) <- c("Elec", "RMSE", "ncomp")

results <- foreach(y=samp.list, .combine=list, .multicombine=TRUE) %dopar%{
  
  return.results <- function(y){
    library(pls) #library for model fitting
    set.seed(y) # set seed
    num.elec <- 15 # total number of electrodes

    best.comp <- comp.vec[y]
    results.vec <- c(y, NA, NA)
    
    # Load Data
    full <- read.csv("combined")
    #full <- read.csv(paste0("sub_samp", samp.size))
    X.colmns <- (ncol(full) - 999):ncol(full)
    
    test.elec <- full[full$Electrode == y, ] # set test electrode
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
    
    test.elec <- test.elec[, (ncol(test.elec) - 1000):ncol(test.elec)] # set test electrode
    
    colnames(test.elec) <- c("DA", paste0("V",1:(ncol(test.elec) - 1)))

    full <- full[, (ncol(full) - 1000):ncol(full)] # set test electrode
    colnames(full) <- c("DA", paste0("V",1:(ncol(full) - 1)))
    
    
    ############
    # Tuning
    ############
    if (tune.metric == 5){
      print("tuning 5")
      # keep CV data
      MSE.mat <- matrix(NA, ncol = max.comp, nrow = length(train.elecs))
      colnames(MSE.mat) <- c(paste0("MSE-Parameter: ", 1:max.comp))
      
               #iterate through electrodes
          for ( x in 1:length(train.elecs) ){
                tune.fit <- pcr(DA ~., data = full[indx[[x]], ], ncomp = max.comp,
                model = FALSE, x = FALSE, y = FALSE) # fit classifier with current tuning parameter
            tune.fit <- fatTrim(tune.fit)
            
            for(i in 1:max.comp){
                preds <- predict(tune.fit, full[test.indx[[x]], -1], ncomp = i)
                
                MSE.mat[x, i] <- sqrt(mean((full$DA[ test.indx[[x]] ] - preds)^2)) # save MSEs of each held out electrode for each tuning parameter
            
            }
      }
      print(MSE.mat)
      
      MSE.mean <- colMeans(MSE.mat)
      print(MSE.mean)
      best.comp <- seq(1, max.comp)[which.min(MSE.mean)] # choose best tuning parameter
      
      print(paste0("elec_", y, "_best.comp_", best.comp))
      rm(tune.fit, preds, MSE.mat)
    
    }
    
    pcr.fit <- pcr(DA ~., data = as.data.frame(full), ncomp = best.comp, 
                   model = FALSE, x = FALSE, y = FALSE) # fit classifier
    pcr.fit <- fatTrim(pcr.fit)
    preds <- predict(pcr.fit, test.elec[, -1], ncomp = best.comp) #make predictions
    preds <- preds * I(preds > 0) # truncate
    results.vec[2] <- sqrt(mean((test.elec[,1] - preds)^2))
    results.vec[3] <- best.comp
    print(paste0("PCR Merged 0 Truncate Error: Elec_", y, ": ", results.vec))
  
    return(results.vec)
  }
  
  return.results(y)
}


for(i in 1:samps){
  final.results[i,] <- results[[i]]
}

write.csv(final.results, output.filename)

