# Study Strap AR LR Derivative
# Gabriel Loewinger 8/1/18
# Updated: 3/2/19

start.time <- proc.time()
library(doParallel)
samps <- 15  #number of parallel foreach loops, Can be number of electrodes running for simultaneous studystrap sampling
num.elec <- 15 #number of electrode datasets to sample from 
straps <- 75 #study replicates to include: HERE IS MAX STUDY REPLICATES
converge.lim <- 300000
bag.size <- 75 # number of random selection with replacement samples to put into a bag 
eta <- 1 #learning rate
best.comp <- 20
tune.fraction <- 1 # proportion of rows to train classifer on to optimize tuning parameters
tune.grid <- seq(15,50, by = 2)
samp.size <- 2500#"full"
comp.vec <- c(23, 22, 49, 38, 48, 22, 22, 22, 23, 22, 22, 23, 49, 20, 22)

tune.metric <- 5 # which OOB tuning metric to use (use 1 - 4), 0 corresponds to no tuning and just using pre-specified best.comp

set.seed(100)
seed.type <- "y"
output.filename <- paste0("SS PCR AR LR Derivative Opt Res_seed_", 
                          seed.type, "_samp_ ", samps, "_ncomp_", best.comp, "_LR_", eta, 
                              "_strap_", straps, "_conv_", converge.lim, "_tuneMet_", tune.metric, 
                                  "_tuneFrac_", tune.fraction, "_seqLen_", length(tune.grid), "_smpSz_",
                                        samp.size, "_bgSz_", bag.size)
print(output.filename)

source("Study Strap Functions.R") # load functions
#------------------------------------------------------
print("SS functions loaded")
#############################################################

clust <- makeCluster(samps, outfile="")
registerDoParallel(clust)

final.results <- matrix(NA, ncol = 15, nrow = samps)
colnames(final.results) <- c("Mean", "Weighted", "Weighted Distance", "OOB Error All Rows", "OOB Elecs Merged", 
                             "Within Bag Feats", "Within Bag RMSE", "Weighted Bag RMSE", "Within Strap Feats", 
                             "Withinstrap test error", "Indiv Strap vs Test Elec Error", "OOB test error", 
                             "OOB Weighted", "In bag elecs-OOB rows", "Total Straps")

results <- foreach(y=1:samps, .combine=list, .multicombine=TRUE) %dopar%{
  
  return.results <- function(y){
    best.comp <- comp.vec[y]
    
    library(pls) #library for model fitting
    study.seed <- y
    set.seed(study.seed) # set seed

    # results matrix
    final.results <- matrix(NA, ncol = 15, nrow = samps)
    colnames(final.results) <- c("Mean", "Weighted", "Weighted Distance", "OOB Error All Rows", "OOB Elecs Merged", 
                                 "Within Bag Feats", "Within Bag RMSE", "Weighted Bag RMSE", "Within Strap Feats", 
                                 "Withinstrap test error", "Indiv Strap vs Test Elec Error", "OOB test error", 
                                 "OOB Weighted", "In bag elecs-OOB rows", "Total Straps")
    print(paste0("Electrode ", y))
    
    # Load Data
    full <- read.csv(paste0("sub_samp", samp.size))
    #full <- as.data.frame(read.csv("combined")) # used for full data
    
    colmns <- (ncol(full)-1000):ncol(full) # columns of full including DA (outcome) and design matrix
    X.colmns <- (ncol(full)-999):ncol(full) # columns of full just the design matrix
    
    # make average CV matrix where first column is electrode number and each row is an average electrode
    CVs <- matrix(NA, ncol = 1002, nrow = num.elec)
    print("CVs matrix")
    for(e in seq(1,num.elec)){ # do not average the CV being tested
      CVs[e, ] <- c(e, colMeans(full[full$Electrode == e, (ncol(full) - 1000):ncol(full)]))
    }
    colnames(CVs) <- c("Electrode", "DA", paste0("V_", 1:1000))
    
    CVs[y,] <- NA # make row corresponding to test electrode NAs to ensure it is not used below (but keep the row in for ease of calculations below) 
    
    # set test and full sets
    ## test electrode is the held out dataset
    ## full is the rest of the data from which the study straps are samp,ed
    print("test elec matrix")
    test.elec <- full[full$Electrode == y, (ncol(full) - 1000):ncol(full)] # set test electrode
    print("full matrix")
    full <- full[full$Electrode != y, ] # take full dataset out of test electrode # took out dataframe
      
    # similarity metrics
    test.similarity <- wiggle.fn( colMeans(test.elec[,-1]) ) #for comparing with other electrodes: same as point.finder()[[9]]
    
    avg.CV <- colMeans(test.elec[,-1])
    print("test.simialrity")
 
    # matricies for storing data
    RF.preds <- matrix(NA, nrow = nrow(test.elec), ncol = straps) # column for each study strap replicate
    similarity.matrix <- matrix(NA, ncol = 8, nrow = straps) # row for each study strap replicate
    avg.CV.matrix <- matrix(NA, ncol = 1000, nrow = straps) #matrix of average electrode CVs, study strap replicate
    r <- 0 #counter
    
    
    # Derivative
    test.elec <- as.data.frame(cbind(test.elec[,1], t(diff(t(test.elec[,-1])))))
    colnames(test.elec) <- c("DA", paste0("V_",1:(ncol(test.elec) - 1)))
    
    ################################################################################  
    # Tuning - Metric 5 --- use CV to choose best parameter and use for all straps
    ################################################################################  
    if (tune.metric == 5){
      print("tune metric 5")
      #########
      # CV indx sets
      #########
      train.elecs <- unique(full$Electrode)
      indx <- list(length = length(train.elecs))
      test.indx <- list(length = length(train.elecs))
      for ( x in 1:length(train.elecs) ){
        test.indx[[x]] <- which(full$Electrode == train.elecs[x])
        rows.indx <- which(full$Electrode != train.elecs[x])
        indx[[x]] <- sample(rows.indx, length(test.indx[[x]]) * tune.fraction, replace = FALSE)
      }
      rm(rows.indx)
      
      # CV data
      MSE.mat <- matrix(NA, ncol = max(tune.grid), nrow = length(train.elecs)) # matrix of test errors associated with only metric used
      colnames(MSE.mat) <- c(paste0("MSE-Parameter: ", 1:max(tune.grid)))
      
      # iterate through tuning parameters
      
        
        #iterate through electrodes
    for ( x in 1:length(train.elecs) ){
          # train classifier on one electrode in training elecs
          tune.samp <- data.frame(cbind(DA = full$DA[indx[[x]]], t(diff(t(full[indx[[x]], X.colmns])))))
          
          tune.fit <- pcr(DA ~., data = tune.samp, ncomp = max(tune.grid),
                          model = FALSE, x = FALSE, y = FALSE) # fit classifier with current tuning parameter
          rm(tune.samp)
          
          tune.fit <- fatTrim(tune.fit)
          
          for ( i in 1:max(tune.grid) ){
              preds <- predict(tune.fit, as.data.frame( t(diff(t(full[test.indx[[x]], X.colmns])))), 
                                                       ncomp = i)
              
              MSE.mat[x, i] <- sqrt(mean((full$DA[ test.indx[[x]] ] - preds)^2)) # save MSEs of each held out electrode for each tuning parameter
        }
      }
      print(MSE.mat)
      
      MSE.mean <- colMeans(MSE.mat)
      print(MSE.mean)
      best.comp <- which.min(MSE.mean) # choose best tuning parameter
      print(best.comp)
      rm(tune.fit, preds, MSE.mat, indx, test.indx)
      
    }
    
    
    ################################################################## 
    # PCR individual - Test Electrode Within Electrode Test Error
    ################################################################## 
    train.rows <- sample(1:nrow(test.elec), nrow(test.elec)/2, replace = FALSE) #training rows
    print(paste0(" within test elec PCR ", y))
    pcr.fit <- pcr(DA ~., data = test.elec, ncomp = best.comp, subset = train.rows,
                   model = FALSE, y = FALSE) # took out dataframe
    pcr.fit <- fatTrim(pcr.fit) # save memory
    final.results[y,6] <- best.comp #save the number of features used
    
    pcr.pred <- predict(pcr.fit, test.elec[-train.rows, -1], ncomp = best.comp) #make predictions
    rm(pcr.fit) # save memory
    print(paste0(" within test elec test error ", y))
    final.results[y,7] <- sqrt(mean((test.elec[-train.rows,1] - pcr.pred)^2)) #within electrode test error on held out set
    print(final.results[y,7])
    rm(train.rows) # save memory
    
    ################################################################## 
    
    strap.results <- matrix(NA, ncol = 7, nrow = straps) # matrix to store individual results for the for loop below
    colnames(strap.results) <- c("Within Strap Feats", "Withinstrap test error", "Indiv Strap vs Test Elec Error", 
                                 "OOB test error", "OOB test - All rows", "OOB elecs merged", "In bag elecs-OOB rows")
    total.straps <- 0 # this is the counter for total times bootstrap sampling has occured
    AR.strap.vec <- c() # vector of number of bootstrap samples it takes between acceptances
    counter <- 0 # between acceptance counter
    z <- 1 # number of pseudo electrodes that are accepted -- need for indexing matrix
    # While loop for the Study Strap sampling goes until convergence limit

    # for first run arbitrary values so doesn't run into error for variables not existing
    sim.metric <- 1 # arbitrary so sim.metric/current.sim.metric >= 1
    current.sim.metric <- 1 # arbitrary so sim.metric/current.sim.metric >= 1
    
    while(counter < converge.lim && z <= straps){
      print(paste0("electrode: ", y, "strap #: ", z))
      start.time1 <- proc.time()
      elecs <- seq(1, num.elec)[-y] #take out the electrode that is currently being tested (y)
      r <- r + 1
      
      if( (sim.metric/current.sim.metric) >= 1){  # only reset counter if acceptance was in fact a strictly better (not MH random acceptance)
          counter <- 0 # reset counter at each new sample
      }
      
      ########################
      ######## Boostraping
      ########################
      start.time2 <- proc.time()
      
      print(paste0("straps compiler; elec: ", y, " ; Strap# ", z))
      
      # ensure the bag doesn't have all 14 electrodes for OOB test error below
      sim <- 0 # this is an indicator telling the while loop to break (if = 1)
      
      while (sim == 0 && counter < converge.lim){ #while loop for restrapping: breaks for either acceptance (sim) or accpetance counter limit (converge.lim)
          ind <- 0 # indicator for while loop for bootstrapping 
          # Choose electrodes to be in bag:
          while(ind == 0){
            bag <- sample(elecs, bag.size, replace = TRUE)
            if(length(unique(bag)) < bag.size){
                ind <- 1
            }
          }
          
          strap.table <- as.data.frame(table(bag)) #sample with replacement and put in data table format
          OOB <- setdiff(elecs, unique(bag)) # electrodes that are OOB
          
          #-------------------------------------------------------------------------
          #####################################################################
          # make study strap matrix with average CVs for testing similarity to speed up calculations
          #####################################################################
          # *********** Use average CVs ***********
          # use rows corresponding to studies/electrodes selected in the bag
          study.strap <- CVs[bag, (ncol(CVs)-1000):ncol(CVs)] #cut out everything but DA and voltages
         
      ########################
      # process dataset and add to similarity matrix and check similarity matrix
      similarity.matrix[z, ] <- wiggle.fn( colMeans(study.strap[,-1]) ) #similarity metric (same as point.finder()[[9]])
      rm(study.strap) # remove CVs study strap from above to avoid confusion
      ######## 
      # check to see if current study strap is more similar than previous one
      diff.dist <- (similarity.matrix[z,] - test.similarity) #* c(rep(20,4), rep(1,4)) # multiple voltage potentials by 20
      current.sim.metric <- sum(diff.dist^2) # squared Euclidean Distance between inflection points
      #u <- runif(1) # for MH step
      u <- 1 # no MH step
      if(total.straps == 0){
        sim.metric <- (1 / eta) * current.sim.metric # set the current similarity metric to the current metric times learning rate eta
        counter <- counter + 1
        total.straps <- total.straps + 1 # add to total straps
        z <- z + 1 # current accepted pseudo electrode counter
        sim <- 1 # break the while loop
      } else if( (sim.metric/current.sim.metric) >  u){ # Metropolis Hastings AR# test to see if the similarity metric is at least as good as before
        # if this is not the first bootstrap sample
          sim.metric <- sim.metric - eta * (sim.metric - current.sim.metric) # set the current similarity metric to the current metric
          counter <- counter + 1
          z <- z + 1 # current accepted pseudo electrode counter
          total.straps <- total.straps + 1 # add to total straps
          sim <- 1 # break the while loop
      }else{
        counter <- counter + 1
        total.straps <- total.straps + 1 # add to total straps
        }
      }
      AR.strap.vec <- c(AR.strap.vec, total.straps)
      print(paste0("Electrode ", y, " Total Straps :", total.straps, " Total Predictions: ", z))
      print(paste0("Total Strap vector :", AR.strap.vec))
      #-------------------------------------------------------------------------
      
      # only remake electrode and fit classifier if the while loop above was broken because of an acceptance
      # i.e., SKIP remaking electrode and fitting classifier if while loop was broken because it reached convergence limit
      if(counter < converge.lim && z <= straps){
      
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
        # generate pseudo - study replicate by taking indices selected above and also ut out everything but DA and voltages
        
        ########################
        
      # process dataset and add to similarity matrix and check similarity matrix
      avg.CV.matrix[z,] <- colMeans(full[indx, X.colmns])  #average CVs
      similarity.matrix[z, ] <- wiggle.fn( colMeans(full[indx, X.colmns]) ) #similarity

      end.time2 <- proc.time()
      print(paste0("study strap time, elec ", y, " strap # ", z))
      print(end.time2 - start.time2)
      #-------------------------------------------------------------------------

      ################################################################## 
      # Optimize Tunining Paramaters
      ##################################################################      
      MSE.mat <- matrix(NA, nrow = length(tune.grid), ncol = 5) # matrix of test errors associated with each tuning parameter value
      MSE.mat.metric <- matrix(NA, nrow = length(tune.grid), ncol = 2) # matrix of test errors associated with only metric metric used
      colnames(MSE.mat) <- c("Tuning Parameters", "OOB Elecs", "OOB Rows", "OOB Elecs Merged", "In Bag Elecs, OOB Rows")
      colnames(MSE.mat.metric) <- c("Tuning Parameters", "OOB Error")
      print("tuning parameter OPT")
      
      opt.mat <- matrix(NA, nrow = straps, ncol = 4) # matrix of optimal tuning parameter associated with each metric
      opt.vec <- vector(length = straps) # vector of optimal tuning parameter associated with metric used
      colnames(opt.mat) <- c("OOB Elecs", "OOB Rows", "OOB Elecs Merged", "In Bag Elecs, OOB Rows")
      
      if (tune.metric > 0 && tune.metric < 5){
      
        for ( i in 1:nrow(MSE.mat) ){
            print(paste0("tuning parameter #: ",i))
            print("fit classifier for tuning")
            
            tune.samp <- data.frame(cbind(DA = full$DA[indx.tune], t(diff(t(full[indx.tune, X.colmns])))))
            
            fit.tune <- pcr(DA ~., data = tune.samp, ncomp = tune.grid[i],
                            model = FALSE, y = FALSE) # fit classifier on tune rows with current tuning parameter value
            
            rm(tune.samp)
            fit.tune <- fatTrim(fit.tune) # save memory
            
            if (tune.metric == 1){
              # Metric 1: OOB Elecs Individually
              oob.tune <- c() # test error from OOB samples
              for (e in OOB){
                  oob.elec <- which(full$Electrode == as.numeric(e)) # make OOB electrode
                  
                  oob.pred <- predict(fit.tune, as.data.frame(t(diff(t(full[oob.elec, X.colmns])))), ncomp = tune.grid[i]) # predict on model based on current OOB electrode covariates
                  oob.tune <- c(oob.tune, sqrt(mean(( full$DA[oob.elec] - oob.pred)^2))) # add OOB test electrode test error to vector
              }
              # MSE.mat[i,2] <- mean(oob.tune) # record all OOB metrics
              MSE.mat.metric[i,2] <- mean(oob.tune) # record only one OOB metrics
              rm(oob.tune, oob.elec, oob.pred)
              
            }else if (tune.metric == 2){
                # Metric 2: All OOB Rows
                # take a sub sample of all OOB rows (length equal to the length of the study strap to save memory)
                rows.samp <- sample( seq(1, nrow(full))[-indx], replace = FALSE )  
                tune.preds <- predict(fit.tune, as.data.frame(t(diff(t(full[rows.samp, X.colmns])))), ncomp = tune.grid[i]) 
                MSE.mat.metric[i, 2] <- sqrt(mean(( full$DA[rows.samp] - tune.preds)^2)) #only one metric
                rm(tune.preds, rows.samp)
                
            }else if (tune.metric == 3){
                # Metric 3: All OOB elecs Merged
                # take a sub sample of all OOB elecs merged (length equal to the length of the study strap to save memory)
                oob.rows <- seq(1, nrow(full))[is.element(full$Electrode, OOB)]
                rows.samp <- sample( oob.rows, min(length(indx), length(oob.rows)), replace = FALSE )  
                oob.pred <- predict(fit.tune, as.data.frame(t(diff(t(full[rows.samp, X.colmns])))),
                                    ncomp = tune.grid[i]) # predict on model based on current OOB electrode covariates

                MSE.mat.metric[i, 2] <- sqrt( mean(   (full$DA[rows.samp] - oob.pred)^2  ) )
                rm(oob.pred, rows.samp, oob.rows)
                
            }else if (tune.metric == 4){
            # Metric 4: OOB - Merged all OOB observations on in bag electrodes**********
                in.bag.rows <- which(is.element(full$Electrode, bag)) # all the rows of the in bag electrodes (not necessarily rows used in Study Strap)
                # take the rows of in bag electrodes, but are out-of-sample
                in.out <- setdiff(in.bag.rows, indx) 
                rows.samp <- sample( in.out, min(length(indx), length(in.out)), replace = FALSE ) # sub sample to save memory
                oob.pred <- predict(fit.tune, as.data.frame(t(diff(t(full[rows.samp, X.colmns])))), 
                                    ncomp = tune.grid[i]) # predict on model based on current OOB electrode covariates
                MSE.mat.metric[i, 2] <- sqrt(mean(  (full$DA[rows.samp] - oob.pred)^2  )) #  
                rm(oob.pred, in.bag.rows, in.out, rows.samp) # save memory 
                }
            
            
       }
        
        
        print("opt mat")
        opt.vec[z] <- tune.grid[which.min(MSE.mat.metric[,2])] 
        rm(fit.tune)
        print("best comp selection")
        best.comp <- opt.vec[z] # choose tuning parameter associated with lowest error
        print(paste0("best comp selection elec_", y, "_best.comp: ", best.comp))
      }
      print(paste0("MSE Mat, Elec_", y))
      
      ################################################################## 
      # PCR for Ensembling
      ################################################################## 
      strap.results[r,1] <- best.comp # save the number of features used
      print(paste0(" test elec model ", y))
      strap.results[r, 2] <- 1
      
      # Full Classifier for Study Strap
      # Derivative
      ss <- data.frame(matrix(NA, ncol = 1000, nrow = length(indx)))
      ss[,1] <- full$DA[indx]
      ss[,2:1000] <- t(diff(t(full[indx, (ncol(full) - 999):ncol(full)])))
      colnames(ss) <- c("DA", paste0("V_", 1:999))
      
      pcr.fit <- pcr(DA ~., data = ss, ncomp = best.comp,
                     model = FALSE, y = FALSE) # fit classifier
      
      rm(ss)
      pcr.fit <- fatTrim(pcr.fit) # save memory
      RF.preds[,r] <- predict(pcr.fit, test.elec[,-1], ncomp = best.comp) # make predictions on test electrode
      
      # truncation
      RF.preds[,r] <- RF.preds[,r] * I(RF.preds[,r] > 0)
    
      strap.results[r, 3] <- sqrt(mean((test.elec[, 1] - RF.preds[,r])^2)) # test electrode test error on this model alone
      print("main PCR model and predictions; test error test elec against one model")
      print(strap.results[r, 3]) # print
      print(paste0("Mean Test Error Ensemble (Average all previous predictions) on elec_", y))
      if(r > 1){
        print(sqrt(mean((test.elec[, 1] - rowMeans(RF.preds[,1:r]))^2)))
      }
      ############################
      # OOB Test Error
      ############################
      print("OOB")
      oob.error <- c() # test error from OOB samples
      print(OOB)
      print(paste0("OOB1 elec_", y))
      for (e in OOB){
        oob.elec <- which(full$Electrode == as.numeric(e)) # make OOB electrode
        # took out as.data.frame
        oob.pred <- predict(pcr.fit, t(diff(t(full[oob.elec, X.colmns]))), ncomp = best.comp) # predict on model based on current OOB electrode covariates
        oob.error <- c(oob.error, sqrt(mean((full$DA[oob.elec] - oob.pred)^2))) # add OOB test electrode test error to vector
        print(oob.error)
      }
      strap.results[r, 4] <- mean(oob.error) # mean OOB test error
      rm(oob.elec, oob.pred)
      #********** OOB - include all rows not trained on **********
      # take a sub sample of all OOB rows (length equal to the length of the study strap to save memory)
      print(paste0("OOB2 elec_", y))
      rows.samp <- sample( seq(1, nrow(full))[-indx], length(indx), replace = FALSE )  
      oob.pred <- predict(pcr.fit, t(diff(t(full[rows.samp, X.colmns]))), ncomp = best.comp) 
      strap.results[r, 5] <- sqrt(mean(( as.numeric(full$DA[rows.samp]) - oob.pred)^2))
      rm(rows.samp)
      
      #********** OOB - Merged all OOB electrodes into 1 test set**********
      # take a sub sample of all OOB elecs merged (length equal to the length of the study strap to save memory)
      print(paste0("OOB3 elec_", y))
      oob.rows <- seq(1, nrow(full))[is.element(full$Electrode, OOB)]
      rows.samp <- sample( oob.rows, min(length(indx), length(oob.rows)), replace = FALSE )  
      oob.pred <- predict(pcr.fit, t(diff(t(full[rows.samp, X.colmns]))), ncomp = best.comp) # predict on model based on current OOB electrode covariates
      strap.results[r, 6] <- sqrt(mean((full$DA[rows.samp] - oob.pred)^2)) #  
      rm(rows.samp, oob.rows)
      
      #********** OOB - Merged all OOB observations on in bag electrodes**********
      print(paste0("OOB4 elec_", y))
      in.bag.rows <- which(is.element(full$Electrode, bag)) # all the rows of the in bag electrodes (not necessarily rows used in Study Strap)
      in.out <- setdiff(in.bag.rows, indx) # take the rows of in bag electrodes, but are out-of-sample
      rows.samp <- sample( in.out, min(length(indx), length(in.out)), replace = FALSE ) # sub sample to save memory
      oob.pred <- predict(pcr.fit, t(diff(t(full[rows.samp, X.colmns]))), 
                          ncomp = best.comp) # predict on model based on current OOB electrode covariates
      strap.results[r, 7] <- sqrt(mean(  (full$DA[rows.samp] - oob.pred)^2  )) #  RMSE
      rm(pcr.fit)
      rm(oob.pred, in.bag.rows, in.out, rows.samp) # save memory 
       ############################
     
      print("mean OOB error")
      print(mean(oob.error))
      
      end.time1 <- proc.time()
      print("one study strap replicate time")
      print(end.time1 - start.time1)
      
      }
         }


      if (tune.metric > 0){ # if use OOB CV to tune 
        opt.mat <- opt.mat[rowSums(!is.na(opt.mat)) > 0,] #eliminate NAs
        print(paste0("elec_", y, "_Optimal Tuning Parameters: ", (opt.mat)))
    }
    
    strap.results <- strap.results[rowSums(!is.na(strap.results)) > 0,] #eliminate NA columns
    row.total <- nrow(strap.results) # use strap results because this has the total number of pseudo electrodes fully computed and classifier fitted
    RF.preds <- RF.preds[, 1:row.total]
    similarity.matrix <- similarity.matrix[1:row.total,]
    avg.CV.matrix <- avg.CV.matrix[1:row.total,]
     # end of studystrap for loop
    #---------------------------------------------------------------------
    
    ####################################################################
    # Compile results from study strap sampling and weight predictions
    ####################################################################
    
    print(paste0("elec_", y, "_dim RF Preds: ", dim(RF.preds)))
    print(paste0("elec_", y, "_dim similarity matrix: ", dim(similarity.matrix)))
    print(paste0("elec_", y, "_dim avg.CV.matrix: ", dim(avg.CV.matrix)))
    print(paste0("elec_", y, "_dim strap.results: ", dim(strap.results)))
    #save average results from one electrode's study strap 
    print("stats: stat results")
    print(strap.results)
    strap.results.mean <- colMeans(strap.results)
    final.results[y, 4] <- strap.results.mean[5]
    final.results[y, 5] <- strap.results.mean[6]
    final.results[y, 14] <- strap.results.mean[7]
    final.results[y, 9] <- strap.results.mean[1] # mean number of within bag features/principal components used
    final.results[y, 10] <- strap.results.mean[2]
    final.results[y, 11] <- strap.results.mean[3]
    final.results[y, 12] <- strap.results.mean[4]
    final.results[y, 15] <- z # total straps
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
      diff.dist <- (similarity.matrix[x,] - test.similarity) #* c(rep(20,4), rep(1,4)) # multiple voltage potentials by 20
      print(paste0("diff.dist: ", diff.dist))
      sim.vec[x] <- mean((sum(diff.dist^2))) 
      sim.vec2[x] <- mean(sqrt(sum(diff.dist[1]^2 + diff.dist[5]^2)) + sqrt(sum(diff.dist[2]^2 + diff.dist[6]^2)) + sqrt(sum(diff.dist[3]^2 + diff.dist[7]^2)) + sqrt(sum(diff.dist[4]^2 + diff.dist[8]^2)))
      sim.vec3[x] <- sum(abs(avg.CV - avg.CV.matrix[x]))
      sim.vec4[x] <- sum((avg.CV - avg.CV.matrix[x])^2)
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
   
    # distance between each point in CV squared
    weights <- (1/sim.vec4) / sum(1/sim.vec4)
    weighted.preds <- RF.preds %*% weights
    
    # within bag test error weighted
    weights <- (1/sim.vec5) / sum(1/sim.vec5)
    weighted.preds <- RF.preds %*% weights
    final.results[y,8] <- sqrt((mean((test.elec[,1] - weighted.preds)^2)))
    
    # OOB test error weighted
    weights <- (1/sim.vec6) / sum(1/sim.vec6)
    weighted.preds <- RF.preds %*% weights
    final.results[y,13] <- sqrt((mean((test.elec[,1] - weighted.preds)^2)))
    
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

write.csv(final.results, output.filename)

end.time <- proc.time()
print("Optimized AR Study Strap Time")
print(end.time - start.time)

