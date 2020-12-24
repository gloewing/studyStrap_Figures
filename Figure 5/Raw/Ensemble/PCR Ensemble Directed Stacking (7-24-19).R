library(doParallel)
library(abind)
library(pls)
library(nnls)
set.seed(1)
samp.size <- 2500
seed.type <- "y"
tune.metric <- 6
sub.samp.ind <- 0 # 1 if we subsample
tune.fraction <- 1
output.filename <- paste0("Directed Stacked PCR Results 0 Truncate (7-27-19)", 
                          seed.type, "_tuneMet_", tune.metric, 
                          "_smpSz_", samp.size, "_subSampInd_", sub.samp.ind)
print(output.filename)
bag.size <- 14
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
    #print(i)
    ## volt.mat[i,] <- voltammograms.m[seq(1,17001, by = 1000)[i]:seq(1000,18000, by = 1000)[i],1]
    #infl <- c(FALSE, diff(diff(test.mat[i,])>0)!=0) #inflection points for each CV ---ORIGINAL WAY
    #infl <- which(infl == TRUE) #add inflection point point numbers (voltage potential number) to vector ---ORIGINAL WAY
    #if (sum(diff(sign(diff(test.mat[i,]))) == -2) == 0){
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
      
      
      #infl <- c(FALSE, diff(sign(diff(test.mat[i,]))))
      #infl <- which(infl == -2)
    }
    
    #points.keep <- c()
    #print(i)
    
    #add conditional to determine if base.points and infl are of different lengths and if they arent then use either 0 or basepoints if its shorter (maybe use furthest point if its basepoints)
    w10.vec <- c()
    wb10.vec <- c()
    b10.vec <- c()
    bb10.vec <- c()
    for (x in 1:length(base.points)){
      #print(x)
      if(abs(infl[which.min(abs(infl - base.points[x]))] - base.points[x]) <= wiggle){ # if the biggest difference between a point another point is <=10, keep it
        #print(1)
        #points.keep <- c(points.keep, infl[which.min(abs(infl - y))])
        window0[i,x] <- test.mat[i,infl[which.min(abs(infl - base.points[x]))]] #if the point fits then use the inflection point specific to each CV
        w10.vec <- c(w10.vec, test.mat[i,(infl[which.min(abs(infl - base.points[x]))] - win):(infl[which.min(abs(infl - base.points[x]))] + win)]) #adds values of a window around the inflection point
        # w10.vec <- c(w10.vec, test.mat[i,(infl[which.min(abs(infl - base.points[x])) - window):(infl[which.min(abs(infl - base.points[x]))] + window)])
        
        # print(2)
        window.base[i,x] <- test.mat[i,infl[which.min(abs(infl - base.points[x]))]] #if the point fits then use the inflection point specific to each CV
        wb10.vec <- c(wb10.vec, test.mat[i,(infl[which.min(abs(infl - base.points[x]))] - win):(infl[which.min(abs(infl - base.points[x]))] + win)])
        
        # print(3)
        # base0[i,x] <- test.mat[i, base.points[which.min(base.points - infl[x])]] #if the point fits, then use the closest base point
        # base.base[i,x] <- test.mat[i, base.points[which.min(base.points - infl[x])]]
        base0[i,x] <- test.mat[i, base.points[x]] #if the point fits, then use the closest base point
        b10.vec <- c(b10.vec, test.mat[i, (base.points[x] - win):(base.points[x] + win)])
        
        # print(4)
        base.base[i,x] <- test.mat[i, base.points[x]] #either way use base inflection points
        bb10.vec <- c(bb10.vec, test.mat[i, (base.points[x] - win):(base.points[x] + win)])
      }
      
      else{
        #print("else")
        d <- d + 1
        window0[i,x] <- 0 #0s if inflection points don't fall into window
        w10.vec <- c(w10.vec, rep(0, win * 2 + 1))
        
        #window.base[i,x] <- test.mat[i, base.points[which.min(base.points - infl[x])]] #if the point doesn't fit then use the base inflection points
        window.base[i,x] <- test.mat[i, base.points[x]] #if the point doesn't fit then use the base inflection points
        wb10.vec <- c(wb10.vec, test.mat[i, (base.points[x] - win):(base.points[x] + win)])
        
        
        base0[i,x] <- 0 #0s if inflection points don't fall into window, use 0
        b10.vec <- c(b10.vec, rep(0,win * 2 + 1))
        
        # base.base[i,x] <- test.mat[i, base.points[which.min(base.points - infl[x])]] #if the point fits then use the base inflection points
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



num.elec <-15

# # Parallelization
clust <- makeCluster(num.elec, outfile = "")
registerDoParallel(clust)

final.results <- matrix(NA, ncol = 8, nrow = num.elec)
colnames(final.results) <- c("Mean", "Weighted", "Weighted Distance", "CV Distance", 
                             "CV Distance Squared", "Individual Feats", "Individual RMSE", "Weighted Bag RMSE")
comp.vec <- c(29, 22, 20, 29, 17, 21, 12, 30, 21, 21, 21, 21, 27, 13, 29)

# stacking
stack.results <- foreach(y=1:num.elec, .combine=list, .multicombine=TRUE, .packages = c("pls", "nnls")) %dopar%{
  stack.fn <- function(y){ 
    library(pls)
    library(nnls)
    set.seed(y)
    i <- 4 # wiggle function parameter
    num.elec <- 15
    best.comp <- comp.vec[y]
    #function so output of parallel foreach loop is easy
    if (sub.samp.ind == 1 ){
        current.elec <- read.csv(paste0("sub_samp", samp.size)) # took out dataframe
        
        colmns <- (ncol(current.elec)-1000):ncol(current.elec) # columns of full including DA (outcome) and design matrix
        X.colmns <- (ncol(current.elec)-999):ncol(current.elec) # columns of full just the design matrix
        
        current.elec <- current.elec[current.elec$Electrode == y, colmns]
               
    }else{
        current.elec <- read.csv(paste0("E", y, "_Final"))
        current.elec <- current.elec[, 5:1005] #cut out everything but DA labels and design matrix
        
    }
        #current.elec <- as.data.frame(cbind(current.elec[,1], point.finder(current.elec[,-1], current.elec[,-1])[[i]])) #process data
    train.rows <- sample(1:nrow(current.elec), nrow(current.elec)/2, replace = FALSE) #make training set rows
    colnames(current.elec) <- c("DA", paste0("V",1:(ncol(current.elec)-1)))
    feat.num <- 4 #no opt
    pcr.fit <- 1
    elec.RMSE <- 1# sqrt(mean((current.elec[-train.rows,]$DA - predict(pcr.fit, current.elec[-train.rows, -1], ncomp = best.comp))^2))
   
    mod_RF <- pcr(DA ~., data = as.data.frame(current.elec), ncomp = best.comp, model = FALSE) #fit RF on full dataset
    stack.preds <- c() #vector (will be turned into a matrix): first column are labels and second column are predictions made on covariates of different electrodes on classifier of current.elec
    for (x in 1:num.elec){
      if(x == y){ #if same electrode classifier and design matrix
        stack.preds <- rbind(stack.preds, cbind(x, current.elec[,1], predict(mod_RF, current.elec[,-1], ncomp = best.comp))) #attach labels and predictions made on classifier trained on current elec data and predicted on current.elec because x == y
      }
      else{
          
          #import and process elec data
          if (sub.samp.ind == 1 ){
              predict.elec <- read.csv(paste0("sub_samp", samp.size)) # took out dataframe
              colmns <- (ncol(predict.elec)-1000):ncol(predict.elec) # columns of full including DA (outcome) and design matrix
              X.colmns <- (ncol(predict.elec)-999):ncol(predict.elec) # columns of full just the design matrix
              predict.elec <- predict.elec[predict.elec$Electrode == x, colmns]
          }else{
              predict.elec <- read.csv(paste0("E", x, "_Final"))
              predict.elec <- predict.elec[, 5:1005] #cut out everything but DA labels and design matrix
              }
          
        colnames(predict.elec) <- c("DA", paste0("V",1:(ncol(predict.elec)-1)))
          
        stack.preds <- rbind(stack.preds, cbind(x, predict.elec[,1], predict(mod_RF, predict.elec[,-1], ncomp = best.comp))) #attach labels and predictions made on classifiers trained on current elec data and predicted on predict.elec's design matrix
        }
    }
    
    #Truncate -- Truncate predictions so bigger than or equal to 0
    stack.preds <- stack.preds * I(stack.preds > 0)
    
    return(list(elec.RMSE, stack.preds))
  }
  stack.r <- stack.fn(y) #call function and save it to variable so it is returned from foreach
}


######################################################
## combine matricies from foreach loop
######################################################
#foreach loop above for each electrode where each row is a different electrode and the first column corresponds to the within electrode test error and the remaining columns are the coefficients from a nnls
for(elec in 1:num.elec){
  if (elec == 1){
    stack.preds <- stack.results[[elec]][[2]] #add both the DA labels and the preds columns
    elec.RMSE <- stack.results[[elec]][[1]] #add RMSE values
  }
  else{
    stack.preds <- cbind(stack.preds, stack.results[[elec]][[2]][,3]) #add just the preds columns
    elec.RMSE <- c(elec.RMSE, stack.results[[elec]][[1]]) #add RMSE values
  }
}
print("dim stacked preds")
print(dim(stack.preds))
print(elec.RMSE)

######################################################
## NNLS for Stacking Weights
######################################################
#list of results: 1) RMSEs for within electrode RF classifier; 2) NNLS coefficients without intercept; 3) NNLS coefficients with intercept
nnls.coefs <- c()
nnls.coefs.int <- c()
for(y in 1:num.elec){
  if(y == 1){
    elec.preds <- stack.preds[stack.preds[,1] != y, -c(1,y + 2)]# take out predictions made on model trained on electrode y (column y) and labels from electrode y (rows where first column (electrode) ==y) take out first column because it is just electrode number 
    nnls.coefs <- coef(nnls(A = as.matrix(elec.preds[,-1]), b = as.matrix(elec.preds[,1])))
    nnls.coefs.int <- coef(nnls(A = as.matrix(cbind(1,elec.preds[,-1])), b = as.matrix(elec.preds[,1])))
  }
  else{
  print(paste0("y: ", y))
  elec.preds <- stack.preds[stack.preds[,1] != y, -c(1,y + 2)]# take out predictions made on model trained on electrode y (column y) and labels from electrode y (rows where first column (electrode) ==y) take out first column because it is just electrode number 
  nnls.coefs <- rbind(nnls.coefs, coef(nnls(A = as.matrix(elec.preds[,-1]), b = as.matrix(elec.preds[,1]))))
  nnls.coefs.int <- rbind(nnls.coefs.int, coef(nnls(A = as.matrix(cbind(1,elec.preds[,-1])), b = as.matrix(elec.preds[,1]))))
  }
}

stacked.results <- list(elec.RMSE, nnls.coefs, nnls.coefs.int) #attach the coefficients of nnls to results vec and return that
rm(stack.preds) #save memory so delete matrix of predictions

print("stacked results")
print(stacked.results)
  # put into function that outputs a list where first element is within electrode test errors and second element is the samples x 2 matrix of labels and predictions for each electrode
  # load electrode y
  # build classifier on half and calculate within elec test error, save to vector
  # build full data classifier
  # for loop
  # make preds based on design matrix of each elec and put into a matrix where first column is DA labels and second column is preds
  # rbind this to stack.preds matrix
  # nnls with and without intercept constrained to be 0
  # extract coefficients
  # output list of length 2 (matrix and vector of within elec RMSEs) (or combine into a single vector)
  # end foreach loop
  

# Parallelization
clust <- makeCluster(num.elec, outfile = "")
registerDoParallel(clust)
comp.vec <- c(29, 22, 20, 29, 17, 21, 12, 30, 21, 21, 21, 21, 27, 13, 29)

######################################################
# main parallel foreach loop
######################################################

results <- foreach(y=1:num.elec, .combine=list, .multicombine=TRUE) %dopar%{

    return.results <- function(y, stacked){ #takes in stacked matrix
      set.seed(y)
      source("Study Strap Functions.R") # load functions
      i <- 4
      num.elec <- 15
      best.comp <- comp.vec[y]
      library(pls)
      set.seed(y)
      final.results <- matrix(NA, ncol = 18, nrow = num.elec)
      colnames(final.results) <- c("Mean", "Wiggle Distance", "Wiggle and Stacking", "Stacked no int", 
                                 "Stacked Int", "NA", "Within Elec", "NA", "Stack Standard", "Stack Int Standard",
                                 "Stacked no int Directed", "Stacked Int Directed","Stack Standard Direct", "Stack Int Standard Direct",
                                 "Stacked no int AR", "Stacked Int AR","Stack Standard AR", "Stack Int Standard AR")
      print(paste0("Electrode ", y))
      test.num <- y
      
      # Load Data
      full <- read.csv(paste0("sub_samp", samp.size))
      colmns <- (ncol(full)-1000):ncol(full) # columns of full including DA (outcome) and design matrix
      X.colmns <- (ncol(full)-999):ncol(full) # columns of full just the design matrix
      
      # make average CV matrix where first column is electrode number and each row is an average electrode
      CVs <- matrix(NA, ncol = 1002, nrow = num.elec)
      print("CVs matrix")
      for(e in seq(1,num.elec)){ # do not average the CV being tested
          CVs[e, ] <- c(e, colMeans(full[full$Electrode == e, (ncol(full) - 1000):ncol(full)]))
      }
      full.CVs <- CVs # average CVs including test elec
      
      
      if (sub.samp.ind == 1 ){
          test.elec <- full[full$Electrode == test.num, colmns]
          full <- full[full$Electrode != test.num,]
      }else{
          rm(full)
          test.elec <- as.data.frame(read.csv(paste0("E", test.num, "_Final")))
          test.elec <- test.elec[, 5:1005]
          }
      
      avg.CV <- colMeans(test.elec[,-1])
      test.similarity <- point.finder(test.elec[,-1], test.elec[,-1])[[9]]
        ################################################################## 
        # individual
        ################################################################## 
        
        # within electrode RMSE
        final.results[,7] <- stacked[[1]] #the within electrode RMSE was caclculated above. adds all electrodes data into so can use for weighting below
        
      r <- 0 #counter
      RF.preds <- matrix(NA, nrow = nrow(test.elec), ncol = num.elec-1)
      similarity.matrix <- matrix(NA, ncol = 8, nrow = num.elec-1)
      avg.CV.matrix <- matrix(NA, ncol = 1000, nrow = num.elec-1) #matrix of average electrode CVs
      for (z in seq(1, num.elec)[-y]){
        r <- r + 1

        if (sub.samp.ind == 1 ){
            current.elec <- full[full$Electrode == z, colmns]
        }else{
            current.elec <- data.frame(read.csv(paste0("E", z, "_Final"))[, 5:1005]) #for ensembling
        }
        
        avg.CV.matrix[r,] <- colMeans(current.elec[, -1]) #average CVs
        similarity.matrix[r, ] <- point.finder(current.elec[,-1], current.elec[,-1])[[9]] #similarity
        i <- 4
        current.elec <- as.data.frame(current.elec)

        ################################################################## 
        # PCR for Ensembling
        ################################################################## 
        print("PCR")
        pcr.fit.current <- pcr(DA~., data=as.data.frame(current.elec), ncomp = best.comp, model = FALSE)
        pcr.fit.current <- fatTrim(pcr.fit.current)
        assign(paste0("model_", z), pcr.fit.current) # save model for stacking
        rm(current.elec)
        print(paste0("preds on model ", z, "data from test elec ", y))
        RF.preds[,r] <- predict(pcr.fit.current, test.elec[,-1], ncomp = best.comp) #make predictions on test electrode based on current electrode's model
        print(paste0("preds DONE on model ", z, " trained on ", y))
      }
      if (sub.samp.ind == 0 ){
          # use smaller dataset for variations on stacking
            full <- read.csv(paste0("sub_samp", samp.size))
        }
      
      print("Stack nnls")
      library(nnls)
     
      ############################
      # AR Regression Weights
      ############################
      print("AR.immitate")
      AR.bags <- AR.immitate(full, full.CVs, test.elec = y, bag.size = nrow(full.CVs) - 1, paths = 1, 
                             converge.lim = 300000) #100000
      print( paste0("Number of AR bags on elec_", y, "_", length(AR.bags)) )
      AR.bags <- unlist(AR.bags) # each element in list is a vector of rows so just make into a single vector
      
      # set up stacking regression matrix
      print("AR.immitate stack mat")
      stack.mat <- matrix(ncol = num.elec, nrow = length(AR.bags))
      stack.mat[,1] <- full$DA[AR.bags] # add labels
      studies <- seq(1, num.elec)[-y]
      
      for (SSL in studies ){
          study <- which(studies == SSL)
          stack.mat[, study + 1] <- predict( get(paste0("model_", SSL)), 
                                             full[ AR.bags, X.colmns], ncomp = best.comp ) 
      }
      # delete models to save memory:
      rm(list = paste0("model_", studies))
      rm(full) # delete to save memory
      rm(AR.bags)
      stack.mat <- stack.mat * I(stack.mat > 0)
      ###################
      # Stack Regression Weights 
      ###################
      
      # nnls
      print("Stack nnls AR immitate")
      library(nnls)
      nnls.coefs.AR <- coef(nnls(A = as.matrix(stack.mat[,-1]), b = as.matrix(stack.mat[,1])))
      nnls.coefs.int.AR <- coef(nnls(A = as.matrix(cbind(1,stack.mat[,-1])), b = as.matrix(stack.mat[,1])))
      rm(stack.mat)
      
      nnls.coefs.new <- nnls.coefs.AR
      nnls.coefs.int.new <- nnls.coefs.int.AR
      
      #Truncate
      RF.preds <- RF.preds * I(RF.preds > 0)
      
      
      mean.preds <- rowMeans(RF.preds)
      final.results[y,1] <- sqrt(mean((test.elec[,1] - mean.preds)^2))
      
      print(paste0("weights_", y))
      sim.vec <- vector(length = nrow(similarity.matrix))
      sim.vec2 <- vector(length = nrow(similarity.matrix))
      sim.vec3 <- vector(length = nrow(similarity.matrix))
      sim.vec4 <- vector(length = nrow(similarity.matrix))
      sim.vec5 <- vector(length = nrow(similarity.matrix))
        
        for (x in 1:nrow(similarity.matrix)){
            diff <- similarity.matrix[x,] - test.similarity
            sim.vec[x] <- sum(diff^2)
        }
        
      
          sim.vec3 <- stacked[[2]][y,]
          sim.vec4 <- stacked[[3]][y,]
       # weight predictions and calculate RMSE
        ###
        print("dim rf preds")
        print(dim(RF.preds))
        # weight all electrodes predictions equally
        mean.preds <- rowMeans(RF.preds)
        final.results[y,1] <- sqrt(mean((test.elec[,1] - mean.preds)^2))
        
        # euclidean distance between modes squared
        weights <- (1/sim.vec) / sum(1/sim.vec)
        weighted.preds <- RF.preds %*% weights
        final.results[y,2] <- sqrt((mean((test.elec[,1] - weighted.preds)^2)))
        print("weights")
        
        # Combine Stacking and Similarity Weighting
        weights <- ( (1/sim.vec) / sum(1/sim.vec) ) * sim.vec4[-1]
        weighted.preds <- RF.preds %*% weights
        final.results[y,3] <- sqrt((mean((test.elec[,1] - weighted.preds)^2)))
        
        # stacked (NNLS coefs) no intercept
        weights <- sim.vec3
        weighted.preds <- RF.preds %*% weights
        final.results[y,4] <- sqrt((mean((test.elec[,1] - weighted.preds)^2)))
        
        # stacked (NNLS coefs) no intercept
        weights <- sim.vec4[-1]
        weighted.preds <- RF.preds %*% weights
        final.results[y,5] <- sqrt((mean((test.elec[,1] - weighted.preds)^2)))
        
        # stacked (NNLS coefs) no intercept and Wiggle
        weights <- ( (1/sim.vec) / sum(1/sim.vec) ) * sim.vec3
        weighted.preds <- RF.preds %*% weights
        final.results[y,6] <- sqrt((mean((test.elec[,1] - weighted.preds)^2)))
        
        ### standardizeds stacking
        
        # stacked (NNLS coefs) no intercept standardized
        weights <- sim.vec3 / sum(sim.vec3)
        weighted.preds <- RF.preds %*% weights
        final.results[y,9] <- sqrt((mean((test.elec[,1] - weighted.preds)^2)))
        
        # stacked (NNLS coefs) no intercept
        weights <- sim.vec4[-1] / sum(sim.vec4[-1])
        weighted.preds <- RF.preds %*% weights
        final.results[y,10] <- sqrt((mean((test.elec[,1] - weighted.preds)^2)))
        
        print("stacking bagGen base")
        # stacked (NNLS coefs) no intercept
        print("weights 13")
        weights <- nnls.coefs.new 
        weighted.preds <- RF.preds %*% weights
        final.results[y,11] <- sqrt((mean((test.elec[,1] - weighted.preds)^2)))
        
        # stacked (NNLS coefs) intercept
        print("weights 14")
        weights <- nnls.coefs.int.new 
        weighted.preds <- cbind(1, RF.preds) %*% weights
        final.results[y,12] <- sqrt((mean((test.elec[,1] - weighted.preds)^2)))
        
        # stacked (NNLS coefs) no intercept
        print("weights 15")
        weights <- nnls.coefs.new / sum(nnls.coefs.new)
        weighted.preds <- RF.preds %*% weights
        final.results[y,13] <- sqrt((mean((test.elec[,1] - weighted.preds)^2)))
        
        # stacked (NNLS coefs) intercept
        print("weights 16")
        weights <- nnls.coefs.int.new / sum(nnls.coefs.int.new)
        weighted.preds <- cbind(1, RF.preds) %*% weights
        final.results[y,14] <- sqrt((mean((test.elec[,1] - weighted.preds)^2)))
        
        ##### AR Stacking
        print("AR stacking")
        # stacked (NNLS coefs) no intercept
        print("weights 17")
        weights <- nnls.coefs.AR
        weighted.preds <- RF.preds %*% weights
        final.results[y,15] <- sqrt((mean((test.elec[,1] - weighted.preds)^2)))
        
        # stacked (NNLS coefs) intercept
        print("weights 18")
        weights <- nnls.coefs.int.AR 
        weighted.preds <- cbind(1, RF.preds) %*% weights
        final.results[y,16] <- sqrt((mean((test.elec[,1] - weighted.preds)^2)))
        
        # stacked (NNLS coefs) no intercept
        print("weights 19")
        weights <- nnls.coefs.AR / sum(nnls.coefs.AR)
        weighted.preds <- RF.preds %*% weights
        final.results[y,17] <- sqrt((mean((test.elec[,1] - weighted.preds)^2)))
        
        # stacked (NNLS coefs) intercept
        print("weights 20")
        weights <- nnls.coefs.int.AR / sum(nnls.coefs.int.AR)
        weighted.preds <- cbind(1, RF.preds) %*% weights
        final.results[y,18] <- sqrt((mean((test.elec[,1] - weighted.preds)^2)))
        
        print(final.results[y,])
        
        return(final.results[y,])
    }
    return.results(y, stacked = stacked.results) #has stacked matrix
    
    
    #write.csv(final.results, paste0("Weighted Ensemble Wiggle1 Results P(6-24)_", y))
    
}
#iterate through and collect items from parallel list results
results.mat <- matrix(NA, ncol = 18, nrow = num.elec)
for(i in 1:num.elec){
  results.mat[i,] <- results[[i]]
}
colnames(results.mat) <- c("Mean", "Wiggle Distance", "Wiggle and Stacking", "Stacked no int", 
                             "Stacked Int", "NA", "Within Elec", "NA", "Stack Directed Standard", "Stack Int Directed Standard",
                             "Stacked no int Directed", "Stacked Int Directed","Stack AR Direct", 
                           "Stack Int AR",
                           "Stacked no int AR", "Stacked Int AR","Stack Standard AR", "Stack Int Standard AR")

write.csv(results.mat, output.filename)

