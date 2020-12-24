# Gabe Loewinger
# 8/29/18
# Study Straps
samps <- 15
setwd("/n/home12/gloewinger") # set so can read file
bag.size <- 14 #num.elec - 1 # number of random selection with replacement samples to put into a bag 
num.elec <- 15 #number of electrode datasets to sample from 
library(doParallel)
set.seed(1)
seed.type <- "y"
samp.size <- "full"
tune.fraction <- 1
straps <- 500
comp.vec <- c(29, 22, 20, 29, 17, 21, 12, 30, 21, 21, 21, 21, 27, 13, 29)
output.filename <- paste0("Study Strap PCR Opt Res_seed_", 
                          seed.type, "_samp_ ", samps, 
                          "_strap_", straps, "_bgSz_", bag.size, "sampSz", samp.size)
print(output.filename)

#------------------------------------------------------
# Wiggle Room 
#------------------------------------------------------
# this variation has no base 
point.finder <- function(train.mat, test.mat, win = 2, wiggle = 10){
  #train.mat -- the matrix that forms the base points (can be the same as test.mat)
  # test.mat -- what is analyzed (can be the same as train.mat)
  # win      -- the window around the points that are included in the output in addition to the base points/wiggle points
  # wiggle   -- the window around which inflection points are detected
  # takes in a matrix of CVs: n x 1000
  # averages the CVs to use as a base CV for comparison
  # generates a window of 5 points to the left and right 
  # and if a new CVs inflection points fall in that range it uses either those inflection points
  # or the inflection points of the base CV otherwise use 0
  # 
  train.mat <- as.matrix(train.mat) #necessary for diff function
  test.mat <- as.matrix(test.mat) #necessary for diff function
  base.CV <- colMeans(train.mat) #average CV
  infl <- c(FALSE, diff(diff(base.CV)>0)!=0) #inflection points of average CV ---ORIGINAL WAY
  #infl <- diff(sign(diff(base.CV)))
  #base.points <- which(infl == -2) #inflection points of average CV
  base.points <- which(infl == TRUE) #inflection points of average CV
  #-------------------------
  # Eliminate similar points within window of 30 # added this in
  #-------------------------
  inf.points <- c()
  
  for(y in 1:length(base.points)){
    if (any(abs(base.points[y] - base.points) <= 30)){
      if (any(abs(base.points[y] - inf.points) <= 30)){ #make sure there arent already these values in the new list
      }
      else{
        inf.points <- c(inf.points, mean(base.points[which(abs(base.points[y] - base.points) <= 30)], base.points[y])) #if any are within a window of 30 take the mean
      }
    }
    else{
      inf.points <- c(inf.points, base.points[y])
    }
  }
  base.points <- inf.points #change it back for consistency below
  #-------------------------
  
  window0 <- matrix(NA, nrow = nrow(test.mat), ncol = length(base.points)) #- if a point falls in a window,  use the current at that point, and 0 if it falls outside the window
  window.base <- matrix(NA, nrow = nrow(test.mat), ncol = length(base.points)) #- if a point falls in a window, and use the current at the base point if it falls outside the window
  base0 <- matrix(NA, nrow = nrow(test.mat), ncol = length(base.points)) #- use the base point if any inflection point falls in a base point window, 0s otherwise
  base.base <- matrix(NA, nrow = nrow(test.mat), ncol = length(base.points)) #- use base points regardless
  
  #gives window of 10 voltage points around modes
  window10 <- matrix(NA, nrow = nrow(test.mat), ncol = length(base.points) * (win * 2 + 1)) #- if a point falls in a window,  use the current at that point, and 0 if it falls outside the window
  window.base10 <- matrix(NA, nrow = nrow(test.mat), ncol = length(base.points) * (win * 2 + 1)) #- if a point falls in a window, and use the current at the base point if it falls outside the window
  base10 <- matrix(NA, nrow = nrow(test.mat), ncol = length(base.points) * (win * 2 + 1)) #- use the base point if any inflection point falls in a base point window, 0s otherwise
  base.base10 <- matrix(NA, nrow = nrow(test.mat), ncol = length(base.points) * (win * 2 + 1)) #- use base points regardless
  
  
  window.points <- c()
  for (x in 1:length(base.points)){
    window.points <- c(window.points, (base.points[x] - wiggle):(base.points[x] + wiggle)) #create window of points where
    #inflection points can fall
  }
  d <- 0
  skips <- 0 # counts the number of times there are no inflection points in the curve and use previous CV's inflection points
  infl.vec <-c()
  for (i in 1:nrow(test.mat)){ #NEED TO FIGURE OUT WHAT TO DO IF INFL is SHORTER THAN BASE.POINTS
     if (length(which(c(FALSE, diff(diff(test.mat[i,])>0)!=0)==TRUE)) == 0){
      skips <- skips + 1 #if there are no inflection points then use the inflection points from the previous CV
    }
    else{
      infl <- c(FALSE, diff(diff(test.mat[i,])>0)!=0) #inflection points for each CV ---ORIGINAL WAY
      infl <- which(infl == TRUE) #add inflection point point numbers (voltage potential number) to vector ---ORIGINAL WAY
      #-------------------------
      # Eliminate similar points within window of 30 # added this in
      #-------------------------
      inf.points <- c()
      
      for(y in 1:length(infl)){
        if (any(abs(infl[y] - infl) <= 30)){
          if (any(abs(infl[y] - inf.points) <= 30)){ #make sure there arent already these values in the new list
          }
          else{
            inf.points <- c(inf.points, mean(infl[which(abs(infl[y] - infl) <= 30)], infl[y])) #if any are within a window of 30 take the mean
          }
        }
        else{
          inf.points <- c(inf.points, infl[y])
        }
      }
      infl <- inf.points #change it back for consistency below
      
    }
    

    #add conditional to determine if base.points and infl are of different lengths and if they arent then use either 0 or basepoints if its shorter (maybe use furthest point if its basepoints)
    w10.vec <- c()
    wb10.vec <- c()
    b10.vec <- c()
    bb10.vec <- c()
    for (x in 1:length(base.points)){
      #print(x)
      if(abs(infl[which.min(abs(infl - base.points[x]))] - base.points[x]) <= wiggle){ # if the biggest difference between a point another point is <=10, keep it

        window0[i,x] <- test.mat[i,infl[which.min(abs(infl - base.points[x]))]] #if the point fits then use the inflection point specific to each CV
        w10.vec <- c(w10.vec, test.mat[i,(infl[which.min(abs(infl - base.points[x]))] - win):(infl[which.min(abs(infl - base.points[x]))] + win)]) #adds values of a window around the inflection point
        window.base[i,x] <- test.mat[i,infl[which.min(abs(infl - base.points[x]))]] #if the point fits then use the inflection point specific to each CV
        wb10.vec <- c(wb10.vec, test.mat[i,(infl[which.min(abs(infl - base.points[x]))] - win):(infl[which.min(abs(infl - base.points[x]))] + win)])
        base0[i,x] <- test.mat[i, base.points[x]] #if the point fits, then use the closest base point
        b10.vec <- c(b10.vec, test.mat[i, (base.points[x] - win):(base.points[x] + win)])
        
        # print(4)
        base.base[i,x] <- test.mat[i, base.points[x]] #either way use base inflection points
        bb10.vec <- c(bb10.vec, test.mat[i, (base.points[x] - win):(base.points[x] + win)])
      }
      
      else{
        d <- d + 1
        window0[i,x] <- 0 #0s if inflection points don't fall into window
        w10.vec <- c(w10.vec, rep(0, win * 2 + 1))
        window.base[i,x] <- test.mat[i, base.points[x]] #if the point doesn't fit then use the base inflection points
        wb10.vec <- c(wb10.vec, test.mat[i, (base.points[x] - win):(base.points[x] + win)])
        
        
        base0[i,x] <- 0 #0s if inflection points don't fall into window, use 0
        b10.vec <- c(b10.vec, rep(0,win * 2 + 1))
        base.base[i,x] <- test.mat[i, base.points[x]] #either way use the base inflection points
        bb10.vec <- c(bb10.vec, test.mat[i, (base.points[x] - win):(base.points[x] + win)])
        
      }
    }
    
    # window points added to matrix after done iterating through base points
    window10[i, ] <- w10.vec
    window.base10[i,] <- wb10.vec
    base10[i,] <- b10.vec
    base.base10[i,] <- bb10.vec
  }
  
  print("Point counter")
  print(d)
  print(paste0("skip counter: ", skips))
  return(list(window0, window.base, base0, base.base, window10, window.base10, base10, base.base10, c(base.points, base.CV[base.points]), base.CV))
}
print("wiggle loaded")
#############################################################
clust <- makeCluster(samps, outfile="")
registerDoParallel(clust)
# results matrix
final.results <- matrix(NA, ncol = 14, nrow = samps)
colnames(final.results) <- c("Mean", "Weighted", "Weighted Distance", "CV Distance", "Volt Dist", 
                             "Within Bag Feats", "Within Bag RMSE", "Weighted Bag RMSE", "Within Strap Feats", 
                             "Withinstrap test error", "Indiv Strap vs Test Elec Error", "OOB test error", "OOB Weighted",
                             "Stacking")

results <- foreach(y=1:samps, .combine=list, .multicombine=TRUE) %dopar%{
  
  return.results <- function(y){
    
    # set parameters
    library(pls) #library for model fitting
    library(nnls)
    set.seed(y) # set seed
    i <- 4 # wiggle function parameter
    best.comp <- comp.vec[y] # components optimized elsewhere
    
    # results matrix
    final.results <- matrix(NA, ncol = 14, nrow = samps)
    colnames(final.results) <- c("Mean", "Weighted", "Weighted Distance", "CV Distance", "Volt Dist", 
                                 "Within Bag Feats", "Within Bag RMSE", "Weighted Bag RMSE", "Within Strap Feats", 
                                 "Withinstrap test error", "Indiv Strap vs Test Elec Error", "OOB test error", "OOB Weighted",
                                 "Stacking")
    print(paste0("Electrode ", y))
    
    # Load Data
     full <- as.data.frame(read.csv("combined"))
    #full <- read.csv(paste0("sub_samp", samp.size)) # took out dataframe
    test.elec <- as.data.frame(full[full$Electrode == y, (ncol(full) - 1000):ncol(full)]) # set test electrode
    full <- as.data.frame(full[full$Electrode != y, ]) # take full dataset out of test electrode
    colmns <- (ncol(full)-1000):ncol(full) # columns of full including DA (outcome) and design matrix
    X.colmns <- (ncol(full)-999):ncol(full) # columns of full just the design matrix
    
    # similarity metrics
    test.similarity <- point.finder(test.elec[,-1], test.elec[,-1])[[9]] #for comparing with other electrodes
    avg.CV <- colMeans(test.elec[,-1])
    print("test.simialrity")
 
    # matricies for storing data
    RF.preds <- matrix(NA, nrow = nrow(test.elec), ncol = straps) # column for each study strap replicate
    weighted.predictions <- matrix(NA, nrow = nrow(test.elec), ncol = straps) # column for each study strap replicate
    similarity.matrix <- matrix(NA, ncol = 8, nrow = straps) # row for each study strap replicate
    avg.CV.matrix <- matrix(NA, ncol = 1000, nrow = straps) #matrix of average electrode CVs, study strap replicate
    r <- 0 #counter
    
    # stacking subset
    stackIndx <- c()
    for(elecIndx in seq(1,15)[-y]){
      # electrodes
      rowIndx <- which(full$Electrode == elecIndx)
      
      # unique DA concentrations
      for(das in unique(full$DA[rowIndx])){
        daIndx <- which(full$DA == das & full$Electrode == elecIndx)
        
        stackIndx <- c(stackIndx,
                       sample(daIndx, round(length(daIndx) / 10), replace = FALSE)
        )
        
      }
    }
    
    rm(daIndx, rowIndx)
    stackMat <- matrix(nrow = length(stackIndx), ncol = straps) # matrix to store predictions for stacking

    colnames(test.elec)[1] <- "DA"
    
    ################################################################## 
    # PCR individual - Test Electrode Within Electrode Test Error
    ################################################################## 
    train.rows <- sample(1:nrow(test.elec), nrow(test.elec)/2, replace = FALSE) #training rows
    pcr.fit <- pcr(DA ~., data = as.data.frame(test.elec[train.rows,]), ncomp = best.comp,
                   model = FALSE, y = FALSE, x = FALSE)
    final.results[y,6] <- best.comp #save the number of features used
    
    pcr.pred <- predict(pcr.fit, test.elec[-train.rows, -1], ncomp = best.comp) #make predictions
    print(paste0(" within test elec test error ", y))
    final.results[y,7] <- sqrt(mean((test.elec[-train.rows,1] - pcr.pred)^2)) #within electrode test error on held out set
    rm(pcr.fit)
    print(final.results[y,7])
    ################################################################## 
    
    strap.results <- matrix(NA, ncol = 4, nrow = straps) # matrix to store individual results for the for loop below
    colnames(strap.results) <- c("Within Strap Feats", "Withinstrap test error", "Indiv Strap vs Test Elec Error", "OOB test error")
    
    # for loop for the Study Strap sampling
    
    for (z in 1:straps){
      print(paste0("electrode: ", y, "strap #: ", z))
      start.time1 <- proc.time()
      elecs <- seq(1, num.elec)[-y] #take out the electrode that is currently being tested (y)
      r <- r + 1
      
            
      ########################
      ######## Boostraping
      ########################
      
      print("straps compiler")
      print(y)
      
      bag <- sample(elecs, bag.size, replace = TRUE)
      strap.table <- as.data.frame(table(bag)) #sample with replacement and put in data table format
      # OOB <- setdiff(elecs, unique(bag)) # electrodes that are OOB
      counter <- 0
      # make study strap matrix
      ########################################
      # remake pseudo electrode with data
      ########################################
      indx <- c() # vector of indices corresponding to rows being sub-sampled
      indx.tune <- c() # vector of indices corresponding to rows being sub-sampled for tuning
      for(i in 1:nrow(strap.table)){
          
          sub.elec <- as.numeric(as.character(strap.table[i,1])) # current electrode being sampled from
          elec.indx <- which(full$Electrode == sub.elec) # rows corresponding to electrode sampled in current bag
          num.obs <- round( (length(elec.indx) / bag.size) * as.numeric(as.character(strap.table[i,2])) )
          elec.indx <- elec.indx[sample(1:length(elec.indx), num.obs, replace = FALSE)] #sub-sample rows
          elec.tune.indx <- sample(elec.indx, round(length(elec.indx) * tune.fraction), replace = FALSE) # sub sample rows to train classifier on
          indx <- c(indx, elec.indx ) # sample as many rows as indicated in elec.obs and add to indices vector
          indx.tune <- c(indx.tune, elec.tune.indx)
      }
      print(paste0(length(indx), " --number of rows in strap for classifier for elec_ ", y, ",_ strap: ", z))
      
      rm(sub.elec)
      rm(elec.indx)
      rm(elec.tune.indx)
      ########################
      print("strap complete")
      avg.CV.matrix[z,] <- colMeans(full[indx, X.colmns]) #average CVs
      similarity.matrix[z, ] <- point.finder(full[indx, X.colmns], full[indx, X.colmns])[[9]] #similarity
      ######## 
      
      # wiggle processing
      print("study strap")
      ################################################################## 
      # PCR for Ensembling and Within Bag
      ################################################################## 
      
      print("RF within")
      strap.results[z,1] <- best.comp # save the number of features used
      
      # remove to avoid problems
      #################################
      # Within Bag Test Error
      #################################
      strap.results[z, 2] <- 1 # just avoid any problems
      print(paste0(y, " within test error", strap.results[z,2]))
      
      # Full Classifier for Study Strap
      pcr.fit <- pcr(DA ~., data = full[indx, colmns], 
                     ncomp = best.comp, model = FALSE, y = FALSE, x = FALSE) # fit classifier
      
      RF.preds[,r] <- predict(pcr.fit, test.elec[,-1], ncomp = best.comp) # make predictions on test electrode
      
      # stacking
      stackMat[,r] <- predict(pcr.fit, full[stackIndx,-1], ncomp = best.comp) # make predictions for stacking
      stackMat[,r] <- stackMat[,r] * I(stackMat[,r] > 0) # truncate stacking
      rm(pcr.fit)

      # print test error of all previous straps
      print(paste0("Mean Test Error Ensemble (Average all previous predictions) on elec_", y))
      if(r > 1){
        print(sqrt(mean((test.elec[, 1] - rowMeans(RF.preds[,1:r]))^2)))
      }
      
      ##############################
      # Write weighted Predictions
      ##############################
      print("similarity weighting")
      similarity.vec <- vector(length = r)
      for (x in 1:r){
        diff.vec <- similarity.matrix[x,] - test.similarity
        similarity.vec[x] <- sum(diff.vec^2)
      }
      print("similarity weights complete")
      
      #mean weights
      weights.sim <- (1/similarity.vec) / sum(1/similarity.vec)
      print(paste0("length weights sim",length(weights.sim)))
      print(paste0("dim RF preds", dim(RF.preds[,1:r])))
      
      # check to see if r ==1. if it is then dont do matrix multiplication (unnecessary and dimensions don't conform)
      if (r == 1){
        weighted.preds <- RF.preds[,r]
        print("similarity weighted completed")
        weighted.predictions[,r] <- weighted.preds # write csv
       
     }else{
      weighted.preds <- RF.preds[,1:r] %*% weights.sim
      print("similarity weighted completed")
      weighted.predictions[,r] <- weighted.preds # write csv
           }
      
      # print test error of all previous weighted straps
      print(paste0("Mean Test Error WEIGHTED Ensemble (Average all previous predictions) on elec_", y))
      if(r > 1){
        print(sqrt(mean((test.elec[, 1] - weighted.predictions[,r])^2)))
      }
      
      #------------------------------------------------------------------------------------
      
      strap.results[z, 3] <- sqrt(mean((test.elec[, 1] - RF.preds[,r])^2)) # test electrode test error on this model alone
      print("main PCR model and predictions; test error test elec against one model")
      print(strap.results[z, 3]) # print
      
      ############################
      # OOB Test Error
      ############################
      print("OOB")
      oob.error <- c(1) # test error from OOB samples

      ############################
      # print(oob.error)
      strap.results[z, 4] <- mean(oob.error) # mean OOB test error
      
      print("mean OOB error")
      print(mean(oob.error))
      
      end.time1 <- proc.time()
      print("one study strap replicate time")
      print(end.time1 - start.time1)

    }
    ########################################################
    
    ##############
    # Compile results from study strap sampling and weight predictions
    ##############

    ##############
    # stacking
    ##############
    stackCoef <- coef( nnls(A = as.matrix(cbind(1, stackMat)), b = as.matrix(full$DA[stackIndx])) )
    rm(stackMat, full)
    
    #save average results from one electrode's study strap 
    print("stats: stat results")
    print(strap.results)
    strap.results.mean <- colMeans(strap.results)
    final.results[y, 9] <- strap.results.mean[1] # mean number of within bag features/principal components used
    final.results[y, 10] <- strap.results.mean[2]
    final.results[y, 11] <- strap.results.mean[3]
    final.results[y, 12] <- strap.results.mean[4]
    print(c(final.results[y, 10], final.results[y, 11], final.results[y, 12], strap.results.mean))
    
    mean.preds <- rowMeans(RF.preds)
    final.results[y,1] <- sqrt(mean((test.elec[,1] - mean.preds)^2))
    
    print(paste0("weights_", y, "main test error on average results: ", final.results[y,1]))
    sim.vec <- vector(length = nrow(similarity.matrix))
    sim.vec2 <- vector(length = nrow(similarity.matrix))
    sim.vec3 <- vector(length = nrow(similarity.matrix))
    sim.vec4 <- vector(length = nrow(similarity.matrix))
    sim.vec5 <- vector(length = nrow(similarity.matrix)) # within strap test error
    sim.vec6 <- vector(length = nrow(similarity.matrix)) # OOB test error
    
    for (x in 1:nrow(similarity.matrix)){
      diff <- similarity.matrix[x,] - test.similarity
      print(paste0("diff", diff))
      sim.vec[x] <- mean((sum(diff[1]^2 + diff[5]^2)) + (sum(diff[2]^2 + diff[6]^2)) + (sum(diff[3]^2 + diff[7]^2)) + (sum(diff[4]^2 + diff[8]^2)))
      sim.vec2[x] <- mean(sqrt(sum(diff[1]^2 + diff[5]^2)) + sqrt(sum(diff[2]^2 + diff[6]^2)) + sqrt(sum(diff[3]^2 + diff[7]^2)) + sqrt(sum(diff[4]^2 + diff[8]^2)))
      sim.vec3[x] <- sum(abs(avg.CV - avg.CV.matrix[x]))
      sim.vec4[x] <- mean((sum(diff[1:4]^2)))
    }
      sim.vec5 <- strap.results[,2] # within strap test error
      sim.vec6 <- strap.results[,4] # OOB test error
    
      
      print("sim.vec")
      print(rbind(sim.vec, sim.vec2, sim.vec3, sim.vec4))
      
    #mean weights
    weights <- (1/sim.vec) / sum(1/sim.vec)
    weighted.preds <- RF.preds %*% weights
    final.results[y,2] <- sqrt((mean((test.elec[,1] - weighted.preds)^2)))
    print("weights")
    
    # actual distance
    weights <- (1/sim.vec2) / sum(1/sim.vec2)
    weighted.preds <- RF.preds %*% weights
    final.results[y,3] <- sqrt((mean((test.elec[,1] - weighted.preds)^2)))
    
    # distance between each point in CV
    weights <- (1/sim.vec3) / sum(1/sim.vec3)
    weighted.preds <- RF.preds %*% weights
    final.results[y,4] <- sqrt((mean((test.elec[,1] - weighted.preds)^2)))
    
    # distance between each point in CV squared
    weights <- (1/sim.vec4) / sum(1/sim.vec4)
    weighted.preds <- RF.preds %*% weights
    final.results[y,5] <- sqrt((mean((test.elec[,1] - weighted.preds)^2)))
    
    # within bag test error weighted
    weights <- (1/sim.vec5) / sum(1/sim.vec5)
    weighted.preds <- RF.preds %*% weights
    final.results[y,8] <- sqrt((mean((test.elec[,1] - weighted.preds)^2)))
    
    # OOB test error weighted
    weights <- (1/sim.vec6) / sum(1/sim.vec6)
    weighted.preds <- RF.preds %*% weights
    final.results[y,13] <- sqrt((mean((test.elec[,1] - weighted.preds)^2)))
    
    # stacking weights
    weights <- stackCoef
    weighted.preds <- cbind(1, RF.preds) %*% weights
    final.results[y,14] <- sqrt((mean((test.elec[,1] - weighted.preds)^2)))
    
    print(paste0("final results, electrode: ", y))
    print(final.results[y,])
    return(final.results[y,])
  }
  return.results(y)
  

}
#iterate through and collect items from parallel list results
for(i in 1:samps){
  final.results[i,] <- results[[i]]
}

setwd("/n/home12/gloewinger/bag_size")
write.csv(final.results, output.filename)

