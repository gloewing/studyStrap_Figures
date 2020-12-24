# Study Strap AR LR
# Gabriel Loewinger 8/1/18
# Updated: 3/2/19

timeStart <- Sys.time()
library(doParallel)
samps <- 15  #number of parallel foreach loops, Can be number of electrodes running for simultaneous studystrap sampling
num.elec <- 15 #number of electrode datasets to sample from 
straps <- 75 #study replicates to include: THIS IS MAX STUDY REPLICATES
strap.vec <- c(16,20,17,14,11,12,10,19,16,20,13,12,7,12,16) # what i did most of my experimentation on
stp.ind <- 1 # 1 if use strap vector elements for each electrode
bags.start <- 40 # number of bags to start with before between bag variability maximization step
converge.lim <- 300000
bag.size <- 2 # number of random selection with replacement samples to put into a bag 
eta <- 1 #learning rate
best.comp <- 20
tune.fraction <- 1 # proportion of rows to train classifer on to optimize tuning parameters
tune.grid <- seq(15,35, by = 2)
bag.select <- 1 # 1 is AR, 2 is Integer Program, 
sprd <- -4
tune.metric <- 6 # OOB tuning metric to use (use 1 - 4), 0 corresponds to no tuning and just using pre-specified best.comp
coords <- 27 # 0 - calculates coordinates based on distances, 1 based on electrode inflection point coordinates 2 # bagGen.eps # 5 is dist.exp
paths <- 1
prox.straps <- 0
lambda <- 0 # regularization tuning parameter
radius <- 400
samp.size <- 2500
comp.vec <- c(29, 22, 20, 29, 17, 21, 12, 30, 21, 21, 21, 21, 27, 13, 29) # official vector; changed 5/12/19
BWbag <- 0 # 0 is no between bag selection, 1 is dot product, 2 is l2 norm of difference
stack.ind <- 0
vw <- 1 #volt weights
OOB.ind <- 0
error.type <- 3 # for diversity/accuracy term
alpha <- 0.9992 #0.9992 # for diversity/accuracy term tuning parameter
lambda.vec <- c(3082, 6, 2465, 6, 170, 12, 70, 414, 108,  70, 5, 6, 6, 14, 5, 265 ) # for BagGen3
set.seed(100)
seed.type <- "y+10"
tune.type <- "MrgTune"
output.filename <- paste0("SS PCR AR LR SubSamp Opt Res_seed_", 
                          seed.type, "_samp_ ", samps, 
                          "_strap_", straps, "_conv_", converge.lim, "_tuneMet_", tune.metric, 
                          "_tuneFrac_", tune.fraction,
                          "_bag.sel_", bag.select, "_sprd_", sprd, "_coord_", 
                          coords, "_samSiz_", samp.size, "_bwBag_", BWbag, "_bagSt_", bags.start,
                          "_CF_", radius, "_lamb_", lambda, "_VW_", vw, "_alpha_", alpha, 
                          "err.typ_", error.type, "_strp.vec_", stp.ind, "_bgSz_", bag.size,"_",tune.type,
                          "_paths_", paths, "_prxStrp_", prox.straps)
print(output.filename)

source("Study Strap Functions.R") # load functions
#------------------------------------------------------
print("SS functions loaded")
#############################################################

clust <- makeCluster(samps, outfile="")
registerDoParallel(clust)

final.results <- matrix(NA, ncol = 25, nrow = samps)
colnames(final.results) <- c("Mean", "Weighted", "Weighted Distance", "OOB Error All Rows", "OOB Elecs Merged", 
                             "Bw Bag Similarity", "Within Bag RMSE", "Weighted Bag RMSE", "Within Strap Feats", 
                             "Runtime", "Indiv Strap vs Test Elec Error", "OOB test error", 
                             "OOB Weighted", "In bag elecs-OOB rows", "Total Straps",
                             "Stack: NNLS", "Stack: NNLS Int", "Stack: Ridge", "Stack: Ridge Int", "Stack: np", 
                              "Stack: np1", "Stack: pos", "Stack: pos1", "dist.var", "dist.mean")

results <- foreach(y=1:samps, .combine=list, .multicombine=TRUE) %dopar%{
  
  return.results <- function(y){
    
    best.comp <- comp.vec[y]
    bagMatrix <- matrix(NA, ncol = bag.size, nrow = straps)
    t.Start <- Sys.time()
    library(pls) #library for model fitting
    library(CVXR) # optimization
    study.seed <- y+10
    set.seed(study.seed) # set seed

    # results matrix
    final.results <- matrix(NA, ncol = 25, nrow = samps)
    colnames(final.results) <- c("Mean", "Weighted", "Weighted Distance", "OOB Error All Rows", "OOB Elecs Merged", 
                                 "Bw Bag Similarity", "Within Bag RMSE", "Weighted Bag RMSE", "Within Strap Feats", 
                                 "Runtime", "Indiv Strap vs Test Elec Error", "OOB test error", 
                                 "OOB Weighted", "In bag elecs-OOB rows", "Total Straps",
                                 "Stack: NNLS", "Stack: NNLS Int", "Stack: Ridge", "Stack: Ridge Int", "Stack: np", 
                                 "Stack: np1", "Stack: pos", "Stack: pos1", "dist.var", "dist.mean")
    
    print(paste0("Electrode ", y))
    
    # Load Data

    full <- read.csv(paste0("sub_samp", samp.size)) # subsampled data

    colmns <- (ncol(full)-1000):ncol(full) # columns of full including DA (outcome) and design matrix
    X.colmns <- (ncol(full)-999):ncol(full) # columns of full just the design matrix
    
    # make average CV matrix where first column is electrode number and each row is an average electrode
    CVs <- matrix(NA, ncol = 1002, nrow = num.elec)
    print("CVs matrix")
    for(e in seq(1,num.elec)){ # do not average the CV being tested
      CVs[e, ] <- c(e, colMeans(full[full$Electrode == e, (ncol(full) - 1000):ncol(full)]))
    }
    colnames(CVs) <- c("Electrode", "DA", paste0("V_", 1:1000))
    full.CVs <- CVs # average CVs including test elec
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
    
    if(stp.ind == 1){
        straps <- strap.vec[y]
    }
 
    # matricies for storing data
    RF.preds <- matrix(NA, nrow = nrow(test.elec), ncol = straps) # column for each study strap replicate
    similarity.matrix <- matrix(NA, ncol = 8, nrow = straps) # row for each study strap replicate
    avg.CV.matrix <- matrix(NA, ncol = 1000, nrow = straps) #matrix of average electrode CVs, study strap replicate
    r <- 0 #counter
    
    if (tune.metric == 6){
        best.comp <- comp.vec[y]
    }
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
      MSE.mat <- matrix(NA, ncol = length(tune.grid), nrow = length(train.elecs)) # matrix of test errors associated with only metric used
      colnames(MSE.mat) <- c(paste0("MSE-Parameter: ", tune.grid))
      
      # iterate through tuning parameters
      for ( i in 1:length(tune.grid) ){
        
        #iterate through electrodes
        for ( x in 1:length(train.elecs) ){
          # train classifier on one electrode in training elecs
          
          tune.fit <- pcr(DA ~., data = full[indx[[x]], colmns], ncomp = tune.grid[i],
                          model = FALSE, x = FALSE, y = FALSE) # fit classifier with current tuning parameter
          
          tune.fit <- fatTrim(tune.fit)
          
          preds <- predict(tune.fit, full[test.indx[[x]], X.colmns], 
                           ncomp = tune.grid[i])
          
          MSE.mat[x, i] <- sqrt(mean((full$DA[ test.indx[[x]] ] - preds)^2)) # save MSEs of each held out electrode for each tuning parameter
          
        }
      }
      print(MSE.mat)
      
      MSE.mean <- colMeans(MSE.mat)
      print(MSE.mean)
      best.comp <- tune.grid[which.min(MSE.mean)] # choose best tuning parameter
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
    
    
      bag.mat <- matrix(NA, ncol = bag.size, nrow = straps) # matrix of bags
      bag.recipe <- matrix(NA, ncol = bag.size, nrow = straps)
      d.vector <- c()
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
        d.vector <- c(current.sim.metric, d.vector)
      } else if( (sim.metric/current.sim.metric) >  u){ # Metropolis Hastings AR# test to see if the similarity metric is at least as good as before
        # if this is not the first bootstrap sample
          sim.metric <- sim.metric - eta * (sim.metric - current.sim.metric) # set the current similarity metric to the current metric
          counter <- counter + 1
          z <- z + 1 # current accepted pseudo electrode counter
          total.straps <- total.straps + 1 # add to total straps
          sim <- 1 # break the while loop
          d.vector <- c(current.sim.metric, d.vector)
      }else{
        counter <- counter + 1
        total.straps <- total.straps + 1 # add to total straps
        }
      }
      AR.strap.vec <- c(AR.strap.vec, total.straps)
      print(paste0("Electrode ", y, " Total Straps :", total.straps, " Total Predictions: ", z))
      print(paste0("r_",r ))
      #-------------------------------------------------------------------------
      
      # only remake electrode and fit classifier if the while loop above was broken because of an acceptance
      # i.e., SKIP remaking electrode and fitting classifier if while loop was broken because it reached convergence limit
      if(counter < converge.lim && z <= straps){
          print("bag.mat add bag")
          bag.mat[r, ] <- bag
          print("bag.recipe add 0s")
          bag.recipe[r,] <- 0 # fill in with 0s and fill the rest in below
        
        ########################################
        # remake pseudo electrode with data
        ########################################
        indx <- c() # vector of indices corresponding to rows being sub-sampled
        indx.tune <- c() # vector of indices corresponding to rows being sub-sampled for tuning
        print("generate study strap")
        for(i in 1:nrow(strap.table)){
          
            sub.elec <- as.numeric(as.character(strap.table[i,1])) # current electrode being sampled from
            elec.indx <- which(full$Electrode == sub.elec) # rows corresponding to electrode sampled in current bag
            proportion <- as.numeric(as.character(strap.table[i,2])) # proportion of rows to sample = proportion/bag.size
            num.obs <- floor(length(elec.indx) / (bag.size)) * proportion # number of rows to sample times the number of times the electrode shows up in the current study strap bag
            print("elec.indx")
            elec.indx <- elec.indx[sample.int(length(elec.indx), num.obs, replace = FALSE)] #sub-sample rows
            print("elec.tune.indx")
            elec.tune.indx <- sample(elec.indx, round(length(elec.indx) * tune.fraction), replace = FALSE) # sub sample rows to train classifier on
            indx <- c(indx, elec.indx ) # sample as many rows as indicated in elec.obs and add to indices vector
            indx.tune <- c(indx.tune, elec.tune.indx)
         }
        print(paste0(length(indx), " --number of rows in strap for classifier for elec_ ", y, ",_ strap: ", z))

        rm(sub.elec)
        rm(elec.indx)
        rm(elec.tune.indx)
        bagMatrix[z,] <- bag # add bag
       
       
        
      # process dataset and add to similarity matrix and check similarity matrix
      avg.CV.matrix[z,] <- colMeans(full[indx, X.colmns])  #average CVs
      similarity.matrix[z, ] <- wiggle.fn( colMeans(full[indx, X.colmns]) ) #similarity

      end.time2 <- proc.time()
      print(paste0("study strap time, elec ", y, " strap # ", z))
      print(end.time2 - start.time2)
      #-------------------------------------------------------------------------
      
      # wiggle processing
      print("study strap")
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
            fit.tune <- pcr(DA ~., data = full[, colmns], ncomp = tune.grid[i], subset = indx.tune,
                            model = FALSE, y = FALSE) # fit classifier on tune rows with current tuning parameter value
            fit.tune <- fatTrim(fit.tune) # save memory
            if (tune.metric == 1){
              # Metric 1: OOB Elecs Individually
              oob.tune <- c() # test error from OOB samples
              for (e in OOB){
                  oob.elec <- which(full$Electrode == as.numeric(e)) # make OOB electrode
                  # took out as.data.frame
                  oob.pred <- predict(fit.tune, full[oob.elec, X.colmns], ncomp = tune.grid[i]) # predict on model based on current OOB electrode covariates
                  oob.tune <- c(oob.tune, sqrt(mean(( full$DA[oob.elec] - oob.pred)^2))) # add OOB test electrode test error to vector
              }
        
              MSE.mat.metric[i,2] <- mean(oob.tune) # record only one OOB metrics
              rm(oob.tune, oob.elec, oob.pred)
              
            }else if (tune.metric == 2){
                # Metric 2: All OOB Rows
                
                # take a sub sample of all OOB rows (length equal to the length of the study strap to save memory)
                rows.samp <- sample( seq(1, nrow(full))[-indx], min(length(indx)), replace = FALSE )  
                tune.preds <- predict(fit.tune, full[rows.samp, X.colmns], ncomp = tune.grid[i]) 
                MSE.mat.metric[i, 2] <- sqrt(mean(( full$DA[rows.samp] - tune.preds)^2)) #only one metric
                rm(tune.preds, rows.samp)
                
            }else if (tune.metric == 3){
                # Metric 3: All OOB elecs Merged
                # take a sub sample of all OOB elecs merged (length equal to the length of the study strap to save memory)
                oob.rows <- seq(1, nrow(full))[is.element(full$Electrode, OOB)]
                rows.samp <- sample( oob.rows, min(length(indx), length(oob.rows)), replace = FALSE )  
                oob.pred <- predict(fit.tune, full[rows.samp, X.colmns],
                                    ncomp = tune.grid[i]) # predict on model based on current OOB electrode covariates
                MSE.mat.metric[i, 2] <- sqrt( mean(   (full$DA[rows.samp] - oob.pred)^2  ) )
                rm(oob.pred, rows.samp, oob.rows)
                
            }else if (tune.metric == 4){
            # Metric 4: OOB - Merged all OOB observations on in bag electrodes**********
                in.bag.rows <- which(is.element(full$Electrode, bag)) # all the rows of the in bag electrodes (not necessarily rows used in Study Strap)
                # take the rows of in bag electrodes, but are out-of-sample
                in.out <- setdiff(in.bag.rows, indx) 
                rows.samp <- sample( in.out, min(length(indx), length(in.out)), replace = FALSE ) # sub sample to save memory
                oob.pred <- predict(fit.tune, full[rows.samp, X.colmns], 
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
      
      strap.results[r, 2] <- 1
      # Full Classifier for Study Strap
      # took out as.data.frame
      pcr.fit <- pcr(DA ~., data = full[, colmns], ncomp = best.comp, subset = indx,
                     model = FALSE, y = FALSE) # fit classifier
      pcr.fit <- fatTrim(pcr.fit) # save memory
      RF.preds[,r] <- predict(pcr.fit, test.elec[,-1], ncomp = best.comp) # make predictions on test electrode
      if (stack.ind == 1){
        print("stack assign")
        assign(paste0("study_", y, "_learner_", r), pcr.fit) # save model for prediction for stacking
        assign(paste0("study_", y, "_ss_", r), indx) # save study strap for stacking
        
      }
         
      # truncation
      RF.preds[,r] <- RF.preds[,r] * I(RF.preds[,r] > 0)
      strap.results[r, 3] <- sqrt(mean((test.elec[, 1] - RF.preds[,r])^2)) # test electrode test error on this model alone
      print("main PCR model and predictions; test error test elec against one model")
      print(strap.results[r, 3]) # print
      print(paste0("Mean Test Error Ensemble (Average all previous predictions) on elec_", y, "_strap_", r))
      if(r > 1){ # was r> 1
        print(sqrt(mean((test.elec[, 1] - rowMeans(RF.preds[,1:r]))^2))) # was [,1:r]
      }else{
        print(sqrt(mean((test.elec[, 1] - RF.preds[,1])^2)))
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
        oob.pred <- predict(pcr.fit, full[oob.elec, X.colmns], ncomp = best.comp) # predict on model based on current OOB electrode covariates
        oob.error <- c(oob.error, sqrt(mean((full$DA[oob.elec] - oob.pred)^2))) # add OOB test electrode test error to vector
        print(oob.error)
      }
      strap.results[r, 4] <- mean(oob.error) # mean OOB test error
      rm(oob.elec, oob.pred)
      #********** OOB - include all rows not trained on **********
      # take a sub sample of all OOB rows (length equal to the length of the study strap to save memory)
      print(paste0("OOB2 elec_", y))
      rows.samp <- sample( seq(1, nrow(full))[-indx], length(indx), replace = FALSE )  
      oob.pred <- predict(pcr.fit, full[rows.samp, X.colmns], ncomp = best.comp) 
      strap.results[r, 5] <- sqrt(mean(( as.numeric(full$DA[rows.samp]) - oob.pred)^2))
      rm(rows.samp)
      
      #********** OOB - Merged all OOB electrodes into 1 test set**********
      # take a sub sample of all OOB elecs merged (length equal to the length of the study strap to save memory)
      print(paste0("OOB3 elec_", y))
      oob.rows <- seq(1, nrow(full))[is.element(full$Electrode, OOB)]
      rows.samp <- sample( oob.rows, min(length(indx), length(oob.rows)), replace = FALSE )  
      oob.pred <- predict(pcr.fit, full[rows.samp, X.colmns], ncomp = best.comp) # predict on model based on current OOB electrode covariates
      strap.results[r, 6] <- sqrt(mean((full$DA[rows.samp] - oob.pred)^2)) #  
      rm(rows.samp, oob.rows)
      
      #********** OOB - Merged all OOB observations on in bag electrodes**********
      print(paste0("OOB4 elec_", y))
      in.bag.rows <- which(is.element(full$Electrode, bag)) # all the rows of the in bag electrodes (not necessarily rows used in Study Strap)
      in.out <- setdiff(in.bag.rows, indx) # take the rows of in bag electrodes, but are out-of-sample
      rows.samp <- sample( in.out, min(length(indx), length(in.out)), replace = FALSE ) # sub sample to save memory
      oob.pred <- predict(pcr.fit, full[rows.samp, X.colmns], 
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
      
      print(paste0("Distances from Electrode_", y))
      print(d.vector)

      
    }


    
  print(paste0("complete model fitting elec_", y))  
  ###
  # Remove NA columns or rows since these matricies were made before the number of straps until convergence was determined
  ###
     if (tune.metric > 0 && OOB.ind != 0){ # if use OOB CV to tune 
        opt.mat <- opt.mat[rowSums(!is.na(opt.mat)) > 0,] #eliminate NAs
        print(paste0("elec_", y, "_Optimal Tuning Parameters: ", (opt.mat)))
    }
    
    print("strap results eliminate NAs")
    strap.results <- strap.results[rowSums(!is.na(strap.results)) > 0,] #eliminate 
    
    # include in case straps = 1
    if(class(strap.results) == "matrix"){
      row.total <- nrow(strap.results) # use strap results because this has the total number of pseudo electrodes fully computed and classifier fitted
    }else{
      row.total <- 1
    }
    
    print("row.total")
    print(row.total)
    
    RF.preds <- RF.preds[, 1:row.total]
    similarity.matrix <- similarity.matrix[1:row.total,]
    avg.CV.matrix <- avg.CV.matrix[1:row.total,]
    
    ####################################################################
    # Compile results from study strap sampling and weight predictions
    ####################################################################
    
    ##############
    # Stacking
    ##############
    if(stack.ind == 1){
      print("Stack Matrix and Models")
      ##############
      # Stacking
      ##############
      # this assumes everything is in global environment
      
      for (ss_mat in 1:row.total){ # iterate through Study Strap design matrices
        
        #---------------------------------------------------
        for (SSL in 1:row.total){ # iterate through Single Study Learners (models)
          # stack.mat is main matrix
          # current.ss is current study strap matrix which is appended to main stack.mat
          if ( SSL == 1 ){
            # first learner of the current study strap -- bind labels (Y - outcome) with predictions made with first SSL
            current.ss <- cbind( full$DA[get(paste0("study_", y, "_ss_", ss_mat))], 
                                 predict( get(paste0("study_", y, "_learner_", SSL)), 
                                          full[get(paste0("study_", y, "_ss_", ss_mat)),-1]) )
          }
          else{
            # subsequent learners of the current study strap -- bind predictions made with current SSL
            current.ss <- cbind( current.ss, 
                                 predict( get(paste0("study_", y, "_learner_", SSL)), 
                                          full[get(paste0("study_", y, "_ss_", ss_mat)),-1]) )            
          }
        }
        
        rm(list = paste0("study_", y, "_ss_", ss_mat)) # Delete current Study Strap to save memory
        
        
        ####
        # Create Stack Matrix based on predictions made on a single Study Strap Design Matrix based on models from all SSLs (s)
        ####
        if (ss_mat == 1){
          # if its the first study strap, then make the stack.mat matrix based on the first matrix of predictions
          stack.mat <- current.ss
          rm(current.ss) # save memory
        }
        else{
          # if its not the first study strap, then rbind the stack.mat matrix with current iteration
          stack.mat <- rbind(stack.mat, current.ss)
          rm(current.ss) # save memory
        }
        ###################
        #---------------------------------------------------
        
        # delete current study strap to save memory:
        rm(list = paste0("study_",y,"_ss_", ss_mat))
        
      }
      # end stacking matrix prep.
      
      # delete models to save memory:
      rm(list = paste0("study_", y, "_learner_", 1:straps))
      
      ###################
      # Stack Regression Weights 
      ###################
      
      # nnls
      print("Stack nnls")
      library(nnls)
      nnls.coefs <- coef(nnls(A = as.matrix(stack.mat[,-1]), b = as.matrix(stack.mat[,1])))
      nnls.coefs.int <- coef(nnls(A = as.matrix(cbind(1,stack.mat[,-1])), b = as.matrix(stack.mat[,1])))
      
      print(paste0("dim RF Preds: ", dim(RF.preds)))
      print(paste0("dim stack.mat: ", dim(stack.mat)))
      print(paste0("stack nnls coefs: ", (nnls.coefs)))
      
      library(glmnet)
      print("Stack ridge")
      
      # ridge
      lambdas <- 10^seq(3, -2, by = -.2) # -.1 gives 50 lamba values
      cv_fit <- cv.glmnet(x = (as.matrix(stack.mat[,-1])), y = as.matrix(stack.mat[,1]),
                          family = "gaussian", alpha = 0, lambda = lambdas)
      
      ridge.coefs <- as.vector(coef(cv_fit, s = "lambda.min"))
      
      
      beta.PN <- 1
      beta.PN1 <- 1
      beta.pos <- 1
      beta.pos1 <- 1
      
      rm(stack.mat, cv_fit) # save memory
      print(paste0("stack ridge coefs: ", (ridge.coefs)))
     
    }
    
        # make matrix A of wiggle.fn outputs (each column is a different elec)
        print("make A matrix")
        A <- matrix(NA, nrow = 8, ncol = num.elec)
        colnames(A) <- paste0("Elec_", 1:num.elec)
        row.names(A) <- paste0("Coordinate_", 1:8)
        
        for (cv.ind in 1:num.elec){
            A[,cv.ind] <- wiggle.fn(full.CVs[cv.ind, -c(1,2)])
        }
        
        print("dist vec")
        dist.vec <- c(length = num.elec) # vector of squared l2 norms of distance between each elec and test elec -- 0 for element corresponding to test elec (y)
        for (cv.ind in 1:num.elec){
            dist.vec[cv.ind] <- sign(sum( (A[,y]- A[,cv.ind])[5:7])) * sum(  (A[,cv.ind] - A[,y])^2 )
        }
        
        # bag distance variability
        bagMatrix <- bagMatrix[rowSums(!is.na(bagMatrix)) > 0,] #eliminate NA rows
        var.vec <- c(length = nrow(bagMatrix))
        mean.vec <- c(length = nrow(bagMatrix))
        # bag distance variability
        
        print("bag variability")
        for ( bg in 1:nrow(bagMatrix)){
            var.vec[bg] <- var( dist.vec[bagMatrix[bg,]] ) # the variance of each bag is the variance of the distances of each electrode in the bag (weighted by proportion of rows sampled from each elec)
            mean.vec[bg] <- mean( dist.vec[bagMatrix[bg,]] )
        }
        
        print("bag variability stats")
        dist.var <- mean(var.vec)
        dist.mean <- mean(mean.vec)
        final.results[y,24] <- dist.var
        final.results[y,25] <- dist.mean
    
    ###############################################################################################
    
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
    
    # run time
    t.End <- Sys.time()  
    print("time")
    final.results[y,10] <- as.numeric(difftime(t.End, t.Start, units='mins'))
    
    # between bag similarity 6 
    if(bag.select == 2){
        bag.recipe <- bag.recipe[rowSums(!is.na(bag.recipe)) > 0,]
        print(paste0("bag recipe mat diff_", y))
        final.results[y,6] <- mat_diff(bag.recipe)
    }

    
    #mean weights
    print("weights")
    weights <- (1/sim.vec) / sum(1/sim.vec)
    weighted.preds <- RF.preds %*% weights
    print(sqrt((mean((test.elec[,1] - weighted.preds)^2))))
    final.results[y,2] <- sqrt((mean((test.elec[,1] - weighted.preds)^2)))
    
    # actual distance
    print("weights 1")
    weights <- (1/sim.vec2) / sum(1/sim.vec2)
    weighted.preds <- RF.preds %*% weights
    print(sqrt((mean((test.elec[,1] - weighted.preds)^2))))
    final.results[y,3] <- sqrt((mean((test.elec[,1] - weighted.preds)^2)))
    
    # distance between each point in CV
    print("weights 2")
    weights <- (1/sim.vec3) / sum(1/sim.vec3)
    weighted.preds <- RF.preds %*% weights
    print(sqrt((mean((test.elec[,1] - weighted.preds)^2))))
    #final.results[y,4] <- sqrt((mean((test.elec[,1] - weighted.preds)^2)))
    
    # distance between each point in CV squared
    print("weights 3")
    weights <- (1/sim.vec4) / sum(1/sim.vec4)
    weighted.preds <- RF.preds %*% weights
    print(sqrt((mean((test.elec[,1] - weighted.preds)^2))))
    #final.results[y,5] <- sqrt((mean((test.elec[,1] - weighted.preds)^2)))
    
    # within bag test error weighted
    print("weights 4")
    weights <- (1/sim.vec5) / sum(1/sim.vec5)
    weighted.preds <- RF.preds %*% weights
    print(sqrt((mean((test.elec[,1] - weighted.preds)^2))))
    final.results[y,8] <- sqrt((mean((test.elec[,1] - weighted.preds)^2)))
    
    # OOB test error weighted
    print("weights 5")
    weights <- (1/sim.vec6) / sum(1/sim.vec6)
    weighted.preds <- RF.preds %*% weights
    print(sqrt((mean((test.elec[,1] - weighted.preds)^2))))
    final.results[y,13] <- sqrt((mean((test.elec[,1] - weighted.preds)^2)))
    
    if(stack.ind == 1){
        # stacked (NNLS coefs) no intercept
        print("stacked (NNLS coefs) no intercept")
        print(paste0("dim rf preds: ", dim(RF.preds)))
        print(paste0("length of nnls no int: ", length(nnls.coefs)))
        weights <- nnls.coefs
        weighted.preds <- RF.preds %*% weights
        final.results[y,16] <- sqrt((mean((test.elec[,1] - weighted.preds)^2)))
        
        # stacked (NNLS coefs) intercept
        print("stacked (NNLS coefs) intercept")
        print(paste0("dim rf preds: ", dim(RF.preds)))
        print(paste0("length of nnls with int: ", length(nnls.coefs.int)))
        weights <- nnls.coefs.int
        weighted.preds <- cbind(1, RF.preds) %*% weights
        final.results[y,17] <- sqrt((mean((test.elec[,1] - weighted.preds)^2)))
        
        # stacked (Ridge coefs) no intercept
        print("stacked (Ridge coefs) no intercept")
        print(paste0("dim rf preds: ", dim(RF.preds)))
        print(paste0("length of ridge coefs: ", length(ridge.coefs)))
        weights <- ridge.coefs[-1]
        weighted.preds <- RF.preds %*% weights
        final.results[y,18] <- sqrt((mean((test.elec[,1] - weighted.preds)^2)))
        
        # stacked (Ridge coefs) intercept
        print("stacked (Ridge coefs) intercept")
        weights <- ridge.coefs
        weighted.preds <- cbind(1, RF.preds) %*% weights
        final.results[y,19] <- sqrt((mean((test.elec[,1] - weighted.preds)^2)))
        
        # stacked (Ridge coefs) intercept
        print("stacked (Ridge coefs) intercept")
        weights <- beta.PN
        weighted.preds <- RF.preds %*% weights
        final.results[y,20] <- sqrt((mean((test.elec[,1] - weighted.preds)^2)))
        
        print("stacked (Ridge coefs) intercept")
        weights <- beta.PN1
        weighted.preds <- cbind(1, RF.preds) %*% weights
        final.results[y,21] <- sqrt((mean((test.elec[,1] - weighted.preds)^2)))

        print("stacked (Ridge coefs) intercept")
        weights <- beta.pos
        weighted.preds <- RF.preds %*% weights
        final.results[y,22] <- sqrt((mean((test.elec[,1] - weighted.preds)^2)))
        
        print("stacked (Ridge coefs) intercept")
        weights <- beta.pos1
        weighted.preds <- cbind(1, RF.preds) %*% weights
        final.results[y,23] <- sqrt((mean((test.elec[,1] - weighted.preds)^2)))
    }
        
    print(paste0("final results, electrode: ", y))
    print(final.results[y,])
    return(final.results[y,])
  }
  return.results(y)
  
  
  
}
#iterate through and collect items from parallel list results
print("Compile final results")
for(i in 1:samps){
  final.results[i,] <- results[[i]]
}

write.csv(final.results, output.filename)

timeEnd <- Sys.time()
print("Optimized AR Study Strap Time")
print(difftime(timeEnd, timeStart, units='mins'))

