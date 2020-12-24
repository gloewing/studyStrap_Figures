library(doParallel)
library(CCA)
library(MatrixCorrelation)
library(caret)
library(glmnet)
library(foreach)
library(tidyverse)
library(dplyr)
library(tibble)
library(pls)

source("Study Strap Functions.R")
source("SimFn.R")
source("ss.ar.caret.mSSL.R")
source("studyStrap.predict.caret.mSSL.R")
source("fatTrim.R")
source("SSE.caret.mSSL.R")
source("merge.caret.mSSL.R")
source("ss.caret.mSSL.featBag.R")

save.folder <- "/n/home12/gloewinger/fscv_bagSize"
comp.vec <- c(29, 22, 20, 29, 17, 21, 12, 30, 21, 21, 21, 21, 27, 13, 29) # official vector; changed 5/12/19

totalSims <- 15 # number of electrodes

# Study Strap and Model fitting parameters
bagSzs <- c(1, 2,4,6,8,10,12,14,16,18,24,28,35,56,65,75,85,100,125,250, 1000) #c(1:3, seq(4,18, by = 2), seq(20,100, by = 10), 200, 500, 1000)
# 
bagSize.vec <- bagSzs
ssBagSz <- ARbagSize <- ARBagSize <- 14
ssStrapsTune <- 150 # number of study straps for ss Tuning
ssStrapsTest <- 150 # number of study straps for Test SS
ssStrapsBagSzTest <- 150 # number of study straps for Test SS Bag Size

totalPaths <- 3 # paths for AR Test
convgLim <- 10000 # convergence limit for AR
convgLimBagSize <- 10000 #AR converge limit for AR bag size tune
convgLimBagSizeTest <- 10000 #AR converge limit for AR bag size test

# parallelize

fName <- paste0("FSCV_raw_150_straps_3_paths_convg10000_bagSize_2500")
######### CHANGE NUMBER OF THREADS BEFORE RUNNING REAL THING !!!!!!!!!!!!!!!!!!!!!!!
logfile <- paste0("outputFile_fscv.txt")
writeLines(c(""), file(logfile,'w'))

num.threads <- as.integer( ceiling( totalSims / 5 ) )# round up
threads <- makeCluster(num.threads, outfile=logfile)
registerDoParallel(threads)

# setwd(save.folder)

getDoParWorkers()
timeStart <- Sys.time()

results <- foreach(iterNum = 1:totalSims, .combine = list, .multicombine = TRUE) %dopar%{
        
        test_study <- y <- iterNum
        trainStudies <- seq(1,15)[-iterNum]
        print(paste0("start: ", iterNum))
        library(doParallel)
        library(CCA)
        library(MatrixCorrelation)
        library(caret)
        library(glmnet)
        library(foreach)
        library(tidyverse)
        library(dplyr)
        library(tibble)
        library(pls)
    
        source("Study Strap Functions.R")
        # source("SimFn.R")
        source("ss.ar.caret.mSSL.R")
        source("studyStrap.predict.caret.mSSL.R")
        source("fatTrim.R")
        source("SSE.caret.mSSL.R")
        source("merge.caret.mSSL.R")
        source("ss.caret.mSSL.featBag.R")

        resList <- vector(length = 13,  "list")
        
        # save results
        # setwd(save.folder)
        final.results <- matrix(ncol = 10, nrow = totalSims)
        colnames(final.results) <- c("mrg", "sse", "sse_stack", "sse_cps", "ss", "ss_cps", 
                                     "ss_stack", "ar", "ar_stack", "ar_cps")
        
        rmse.vec <- c() # store results
        
        # file names
        filename <- paste0("fscv_raw")
        baseName <- filename
        fileNm <- filename
        
        arTuneName <- paste0("AR Tune_", baseName, "_",iterNum, "_")
        ssTuneName <- paste0("SS Tune_", baseName, "_",iterNum, "_")
        arBagSizeTest <- paste0("AR BagSize_Test_", baseName, "_",iterNum, "_")
        ssBagSizeTest <- paste0("SS BagSize_Test_", baseName, "_",iterNum, "_")
        
        seedSet <- iterNum # ensures no repeats
        set.seed(seedSet)
        
        # simulate data
        full <- as.data.frame( read.csv("sub_samp2500")  )

        
        colmns <- (ncol(full)-1000):ncol(full) # columns of full including DA (outcome) and design matrix
        X.colmns <- (ncol(full)-999):ncol(full) # columns of full just the design matrix
        
        # set test and full sets
        ## test electrode is the held out dataset
        ## full is the rest of the data from which the study straps are samp,ed
        print("test elec matrix")
        test.elec <- full[full$Electrode == y, (ncol(full) - 1000):ncol(full)] # set test electrode
        print("full matrix")
        full <- full[full$Electrode != y, -c(1,2,4,5)] # take full dataset out of test electrode # took out dataframe
        colnames(full)[1] <- "Study"
        colnames(full)[2] <- colnames(test.elec)[1] <- "Y"

        # similarity measure
        simFn <- function(x1, x2){
            return( 1 / sum(  ( wiggle.fn( colMeans( x1) ) - wiggle.fn( colMeans(x2) ) )^2  )  )
        }
        
        testLambda <- data.frame(ncomp = comp.vec[y]) # previously tuned
        ########################################

        ##############
        # Merge
        ##############
 
        source("merge.caret.mSSL.R")
        
        print(paste0("Merge: ", iterNum))
        mrg.mod <- merged(formula = Y ~., data = full, model = FALSE,
                          ssl.method = list("pcr"),
                            ssl.tuneGrid = list( testLambda ) )
        preds <- studyStrap.predict(mrg.mod, test.elec[,-1])
        rm(mrg.mod)
        rmse <- sqrt(mean( (preds[,1] - test.elec$Y)^2   ))
        rmse.vec <- c(rmse.vec, rmse)
        ##########################################################################################

        ##############
        # SSE
        ##############
        source("SSE.caret.mSSL.R")

        print(paste0("SSE: ", iterNum))
        sse.mod <- sse(formula = Y ~., data = full, target.study = test.elec[,-1], # remove Y from test.elec
                       sim.mets = FALSE, model = FALSE,
                       ssl.method = list("pcr"),
                       ssl.tuneGrid = list( testLambda ), 
                       customFNs = list(simFn),
                       stackTune = FALSE)
                        #,
                       #tune = "cv",
                       #nfolds = 5,
                       #tune.grid = list(tune.grid)
                       #)

        preds <- studyStrap.predict(sse.mod, test.elec[,-1])
        rm(sse.mod)
        rmse.avg <- sqrt(mean( (preds[,1] - test.elec$Y)^2  ))
        rmse.stack <- sqrt(mean( (preds[,2] - test.elec$Y)^2  ))
        rmse.cps <- sqrt(mean( (preds[, ncol(preds)] - test.elec$Y)^2  ))  # last column is custom function
        rmse.vec <- c(rmse.vec, rmse.avg, rmse.stack, rmse.cps)


        # ###########
        # # SS Tune
        # ###########
        # 
        # bagSize.vec <- bagSzs
        # print(paste0("SS Tune: ", iterNum))
        # 
        ssTuneResults <- matrix(NA, ncol = length(bagSize.vec), nrow = length(trainStudies))
        ssTuneResultsCPS <- matrix(NA, ncol = length(bagSize.vec), nrow = length(trainStudies))
        ssTuneResultsStack <- matrix(NA, ncol = length(bagSize.vec), nrow = length(trainStudies))
        # # Load Data
        # original.Train <- trainStudies
        # 
        # #iterate through bag sizes
        # 
        # for(bagSize in bagSzs){
        # 
        #     # iterate through held out studies
        # 
        #     for(y in trainStudies){
        # 
        #         original.Train <- trainStudies
        #         print(paste0("bagSize: ", bagSize))
        # 
        #         original.Train <- original.Train[original.Train != y] # remove validation set from training study vector
        # 
        #         #######################
        #         # Tune Lasso parameter
        #         #######################
        # 
        #         # hold one out study out CV indices of studies
        #         indxList <- vector("list", length = length(original.Train))
        #         HOOList <- vector("list", length = length(original.Train))
        # 
        #         for(study in 1:length(original.Train)){
        #             indxList[[study]] <- which(full$Study != c(original.Train[study], y) )  # indices of studies to train on
        #             HOOList[[study]] <- which(full$Study == original.Train[study]) # indices of study to hold out
        # 
        #         }
        # 
        #         # Tune lambda
        # 
        #         control <- caret::trainControl(method="cv", index = indxList, indexOut = HOOList)
        #         # tuneModel <- caret::train(Y~., data=full[,-1], method="glmnet",
        #         #                           trControl=control, tuneGrid = tune.grid, metric="RMSE")
        #         # 
        #         # best.lambda <- tuneModel$bestTune
        #         # best.lambda <- as.data.frame(best.lambda)
        #         
        #         fullIndx <- which(full$Study != y)
        #         testIndx <- which(full$Study == y)
        # 
        # 
        #         source("ss.caret.mSSL.featBag.R")
        #         ssMod <-  ss(data = full[fullIndx,], formula = Y ~.,
        #                      target.study = full[testIndx,-c(1,2)],
        #                      bag.size = bagSize,
        #                      straps = ssStrapsTune,
        #                      stack = "standard",
        #                      ssl.method = list("lm"),
        #                      #ssl.tuneGrid = list( testLambda ),
        #                      sim.mets = FALSE,
        #                      stack.standardize = FALSE, customFNs = list(simFn)#,
        #                      #tune = "cv",
        #                      #nfolds = 2,
        #                     # tune.grid = list(tune.grid)
        #                     )
        # 
        #         preds <- studyStrap.predict(ssMod, full[testIndx,-1])
        #         preds <- preds[,c(1,2,ncol(preds))] # average, stacking, CPSW
        #         rm(ssMod) # delete model object to save memory
        #         rmse <- apply( (test.elec$Y - preds)^2, 2, FUN = function(x) sqrt(mean(x)))
        #         print(paste0("final results, study: ", y, "_bagSize: ", bagSize, "_RMSE_", rmse))
        #         indx <- which(trainStudies == y) # index for results storage
        #         bag.indx <- which(bagSize.vec == bagSize)
        # 
        #         ssTuneResults[indx, bag.indx] <- rmse[1]
        #         ssTuneResultsCPS[indx, bag.indx] <- rmse[3]
        #         ssTuneResultsStack[indx, bag.indx] <- rmse[2]
        # 
        #     }
        # }
        # 
        # 
        #     # save for output
            resList[[1]] <- colMeans(ssTuneResults)
            resList[[2]] <- colMeans(ssTuneResultsCPS)
            resList[[3]] <- colMeans(ssTuneResultsStack)
        #     
        # ssRMSE <- colMeans(ssTuneResults) # average across held out datasets
        # ssBagSize <- bagSize.vec[which.min(ssRMSE)] # optimal bag size
        # 
        # # save memory
        # rm(ssTuneResults)
        # rm(ssTuneResultsCPS)
        # rm(ssTuneResultsStack)
        ############################################

        ##############
        # Test ss
        ##############
        print(paste0("Test SS: ", iterNum))

        source("ss.caret.mSSL.featBag.R")

        ssMod <-  ss(data = full, formula = Y ~.,
                     target.study = test.elec[,-1],
                     bag.size = ssBagSz,
                     straps = ssStrapsTest,
                     stack = "standard",
                     ssl.method = list("pcr"),
                     ssl.tuneGrid = list( testLambda ),
                     sim.mets = FALSE,
                     stack.standardize = FALSE,
                     customFNs = list(simFn),
                     stackTune = FALSE #,
                     #tune = "cv",
                     #nfolds = 5,
                     #tune.grid = list(tune.grid)
                     )

        preds <- studyStrap.predict(ssMod, test.elec[,-1])
        preds <- preds[,c(1, 2, ncol(preds))] # average, stacking, CPSW
        rm(ssMod) # delete model object to save memory
        rmse.vec <- c(rmse.vec, apply( (test.elec$Y - preds)^2, 2, FUN = function(x) sqrt(mean(x))) )

        #####################################################
        # Bag Sizes Test SS
        ####################################################
        print(paste0("Test SS Bag Size: ", iterNum))
        # Load Data
        original.Train <- trainStudies

        # cycle through bagsizes
        ssBagSize <- matrix(ncol = length(bagSzs), nrow = totalSims)
        ssBagSizeCPS <- matrix(ncol = length(bagSzs), nrow = totalSims)
        ssBagSizeStack <- matrix(ncol = length(bagSzs), nrow = totalSims)

        for(bagSize in bagSzs){

            source("ss.caret.mSSL.featBag.R")

            ssMod <-  ss(data = full, formula = Y ~.,
                         target.study = test.elec[,-1],
                         bag.size = bagSize,
                         straps = ssStrapsBagSzTest,
                         stack = "standard",
                         ssl.method = list("pcr"),
                         ssl.tuneGrid = list( testLambda ),
                         sim.mets = FALSE,
                         stack.standardize = FALSE,
                         customFNs = list(simFn),
                         stackTune = FALSE #,
                         #tune = "cv",
                         #nfolds = 2,
                        # tune.grid = list(tune.grid)
                         )

            preds <- studyStrap.predict(ssMod, test.elec[,-1])
            preds <- preds[,c(1, 2, ncol(preds))] # average, stacking, CPSW
            rm(ssMod) # delete model object to save memory
            rmse <- apply( (test.elec$Y - preds)^2, 2, FUN = function(x) sqrt(mean(x)))
            print(paste0("final results, study: ", y, "_bagSize: ", bagSize, "_RMSE_", rmse))
            bagSzIndx <- which(bagSzs == bagSize)

            ssBagSize[iterNum, bagSzIndx] <- rmse[1]
            ssBagSizeCPS[iterNum, bagSzIndx] <- rmse[3]
            ssBagSizeStack[iterNum, bagSzIndx] <- rmse[2]
        }

        ssBSz <- paste0(ssBagSizeTest, "_", "Avg")
        
        # save for output
        resList[[4]] <- ssBagSize[iterNum,]
        resList[[5]] <- ssBagSizeCPS[iterNum,]
        resList[[6]] <- ssBagSizeStack[iterNum,]
        
        # save memory
        rm(ssBagSize)
        rm(ssBagSizeCPS)
        rm(ssBagSizeStack)

        ###########################################################################################################
# 
#         ###########
#         # AR Tune
#         ###########
        bagSize.vec <- bagSzs
        print(paste0("AR Tune: ", iterNum))

        arTuneResults <- matrix(NA, ncol = length(bagSize.vec), nrow = length(trainStudies))
        arTuneResultsCPS <- matrix(NA, ncol = length(bagSize.vec), nrow = length(trainStudies))
        arTuneResultsStack <- matrix(NA, ncol = length(bagSize.vec), nrow = length(trainStudies))
#         # Load Data
#         original.Train <- trainStudies
# 
#         #iterate through bag sizes
# 
#         for(bagSize in bagSzs){
# 
#             # iterate through held out studies
# 
#             for(y in trainStudies){
# 
#                 original.Train <- trainStudies
#                 print(paste0("bagSize: ", bagSize))
# 
#                 original.Train <- original.Train[original.Train != y] # remove validation set from training study vector
# 
#                 #######################
#                 # Tune Lasso parameter
#                 #######################
# 
#                 # hold one out study out CV indices of studies
#                 indxList <- vector("list", length = length(original.Train))
#                 HOOList <- vector("list", length = length(original.Train))
# 
#                 for(study in 1:length(original.Train)){
#                     indxList[[study]] <- which(full$Study != c(original.Train[study], y) )  # indices of studies to train on
#                     HOOList[[study]] <- which(full$Study == original.Train[study]) # indices of study to hold out
#                     
#                 }
# 
#                 # control <- caret::trainControl(method="cv", number=num.trainStudy, index = indxList, indexOut = HOOList)
#                 # tuneModel <- caret::train(Y~., data=full[,-1], method="glmnet",
#                 #                           trControl=control, tuneGrid = tune.grid, metric="RMSE")
#                 # 
#                 # best.lambda <- tuneModel$bestTune
#                 # best.lambda <- as.data.frame(best.lambda)
#                 
#                 fullIndx <- which(full$Study != y)
#                 testIndx <- which(full$Study == y)
# 
#                 source("ss.ar.caret.mSSL.R")
#                 
#                 arMod <-  cmss(data = full[fullIndx,], formula = Y ~.,
#                                target.study = full[testIndx, -c(1,2)],
#                                converge.lim = convgLimBagSize,#10000,
#                                bag.size = bagSize,
#                                max.straps = 150,
#                                paths = 3,
#                                stack = "standard",
#                                ssl.method = list("lm"),
#                                # ssl.tuneGrid = list( best.lambda ),
#                                sim.mets = FALSE, #model = FALSE,
#                                meanSampling = TRUE,
#                                stack.standardize = FALSE,
#                                sim.fn = simFn,
#                                customFNs = list(simFn)) #,
#                                #tune = "cv",
#                                #nfolds = 2,
#                                #tune.grid = list(tune.grid))
# 
#                 preds <- studyStrap.predict(arMod, full[testIndx, -c(1,2)])
#                 preds <- preds[,c(1, 2, ncol(preds))] # average, stacking, CPSW
#                 rm(arMod) # delete model object to save memory
#                 rmse <- apply( (test.elec$Y - preds)^2, 2, FUN = function(x) sqrt(mean(x)))
#                 print(paste0("final results, study: ", y, "_bagSize: ", bagSize, "_RMSE_", rmse))
#                 indx <- which(trainStudies == y) # index for results storage
#                 bag.indx <- which(bagSize.vec == bagSize)
# 
#                 arTuneResults[indx, bag.indx] <- rmse[1]
#                 arTuneResultsCPS[indx, bag.indx] <- rmse[3]
#                 arTuneResultsStack[indx, bag.indx] <- rmse[2]
# 
#             }
#         }
# 
#         # fill in corresponding row with colMeans
        resList[[7]] <- colMeans(arTuneResults)
        resList[[8]] <- colMeans(arTuneResultsCPS)
        resList[[9]] <-  colMeans(arTuneResultsStack)
# 
#         arRMSE <- colMeans(arTuneResults) # average across held out datasets
#         ARbagSize <- bagSize.vec[which.min(arRMSE)] # optimal bag size
#         
#         rm(arTuneResults)
#         rm(arTuneResultsCPS)
#         rm(arTuneResultsStack)
#         ############################################
#         ##############
#         # Test AR
#         ##############
#         print(paste0("Test AR: ", iterNum))
# 
#         if( ARbagSize == 1){
#             # if it is only 1 then there are very few possible combinations--no need to waste time
#             convgLim <- num.trainStudy * 2
#         }
# 
#         source("ss.ar.caret.mSSL.R")
# 
#         arMod <-  cmss(data = full, formula = Y ~.,
#                        target.study = test.elec[,-1],
#                        converge.lim = convgLim,
#                        bag.size = ARbagSize,
#                        max.straps = 200,
#                        paths = totalPaths,
#                        stack = "standard",
#                        ssl.method = list("lm"),
#                        #ssl.tuneGrid = list( testLambda ),
#                        sim.mets = FALSE, model = FALSE,
#                        meanSampling = TRUE,
#                        stack.standardize = FALSE,
#                        sim.fn = simFn,
#                        customFNs = list(simFn)) #,
#                        #tune = "cv",
#                        #nfolds = 5,
#                        #tune.grid = list(tune.grid),
#                        stackTune = FALSE)
# 
#         preds <- studyStrap.predict(arMod, test.elec[,-1])
#         preds <- preds[,c(1, 2, ncol(preds))] # average, stacking, CPSW
#         rm(arMod) # delete model object to save memory
#         rmse.vec <- c(rmse.vec, apply( (test.elec$Y - preds)^2, 2, FUN = function(x) sqrt(mean(x))) )
# 
#         #####################################################
#         # Bag Sizes Test AR
#         ####################################################
#         print(paste0("Test AR Bag Size: ", iterNum))
#         # Load Data
#         original.Train <- trainStudies
# 
#                 # cycle through bagsizes
        arBagSize <- matrix(ncol = length(bagSzs), nrow = totalSims)
        arBagSizeCPS <- matrix(ncol = length(bagSzs), nrow = totalSims)
        arBagSizeStack <- matrix(ncol = length(bagSzs), nrow = totalSims)
# 
#         for(bagSize in bagSzs){
# 
#             source("ss.ar.caret.mSSL.R")
# 
#             arMod <-  cmss(data = full, formula = Y ~.,
#                            target.study = test.elec[,-1],
#                            converge.lim = convgLimBagSizeTest,
#                            bag.size = bagSize,
#                            max.straps = 200,
#                            paths = 3,
#                            stack = "standard",
#                            ssl.method = list("lm"),
#                            #ssl.tuneGrid = list( testLambda ),
#                            sim.mets = FALSE, #model = FALSE,
#                            meanSampling = TRUE,
#                            stack.standardize = FALSE,
#                            sim.fn = simFn,
#                            customFNs = list(simFn)) #,
#                           # tune = "cv",
#                           # nfolds = 2,
#                            #tune.grid = list(tune.grid),
#                           stackTune = FALSE)
# 
#             preds <- studyStrap.predict(arMod, test.elec[,-1])
#             preds <- preds[,c(1, 2, ncol(preds))] # average, stacking, CPSW
#             rm(arMod) # delete model object to save memory
#             rmse <- apply( (test.elec$Y - preds)^2, 2, FUN = function(x) sqrt(mean(x)))
#             print(paste0("final results, study: ", y, "_bagSize: ", bagSize, "_RMSE_", rmse))
#             bagSzIndx <- which(bagSzs == bagSize)
# 
#             arBagSize[iterNum, bagSzIndx] <- rmse[1]
#             arBagSizeCPS[iterNum, bagSzIndx] <- rmse[3]
#             arBagSizeStack[iterNum, bagSzIndx] <- rmse[2]
#         }
# 
        rmse.vec <- c(rmse.vec, NA, NA, NA)
        
        resList[[10]] <- arBagSize[iterNum,]
        resList[[11]] <- arBagSizeCPS[iterNum,]
        resList[[12]] <- arBagSizeStack[iterNum,]
#         
#         rm(arBagSize)
#         rm(arBagSizeCPS)
#         rm(arBagSizeStack)

        print(paste0("ForeEAch Iteration Complete: _IterNum: ", iterNum ) )
        print(rmse.vec)

        resList[[13]] <-  rmse.vec
        
        return(resList)
        }        

print("Complete, length of results")
print(length(results))

timeEnd <- Sys.time()
print(difftime(timeEnd, timeStart, units='mins'))


### names
filename <- paste0("fscv_raw_Combined")
baseName <- paste0("fscv_raw_Combined")
fileNm <- paste0("fscv_raw_Combined")

arTuneName <- paste0("fscv_raw AR Tune_", baseName)
ssTuneName <- paste0("fscv_raw SS Tune_", baseName)
arBagSizeTest <- paste0("fscv_raw AR BagSize_Test_")
ssBagSizeTest <- paste0("fscv_raw SS BagSize_Test_")

### matrices

# ss tune
ssTuneResults <- matrix(NA, ncol = length(bagSize.vec), nrow = totalSims)
ssTuneResultsCPS <- matrix(NA, ncol = length(bagSize.vec), nrow = totalSims)
ssTuneResultsStack <- matrix(NA, ncol = length(bagSize.vec), nrow = totalSims)

# ss test
ssBagSize <- matrix(ncol = length(bagSzs), nrow = totalSims)
ssBagSizeCPS <- matrix(ncol = length(bagSzs), nrow = totalSims)
ssBagSizeStack <- matrix(ncol = length(bagSzs), nrow = totalSims)
# 
# # ar tune
# arBagSizeTune <- matrix(ncol = length(bagSzs), nrow = totalSims)
# arBagSizeCPSTune <- matrix(ncol = length(bagSzs), nrow = totalSims)
# arBagSizeStackTune <- matrix(ncol = length(bagSzs), nrow = totalSims)
# 
# # ar test
# arBagSize <- matrix(ncol = length(bagSzs), nrow = totalSims)
# arBagSizeCPS <- matrix(ncol = length(bagSzs), nrow = totalSims)
# arBagSizeStack <- matrix(ncol = length(bagSzs), nrow = totalSims)

# rmse - main results
final.results <- matrix(ncol = 10, nrow = totalSims)

### add data
for(i in 1:totalSims){
    # ss tune
    ssTuneResults[i,] <- results[[i]][[1]]
    ssTuneResultsCPS[i,] <- results[[i]][[2]]
    ssTuneResultsStack[i,] <- results[[i]][[3]]
    
    # ss test
    ssBagSize[i,] <- results[[i]][[4]]
    ssBagSizeCPS[i,] <- results[[i]][[5]]
    ssBagSizeStack[i,] <- results[[i]][[6]]
    # 
    # # ar tune
    # arBagSizeTune[i,] <- results[[i]][[7]]
    # arBagSizeCPSTune[i,] <- results[[i]][[8]]
    # arBagSizeStackTune[i,] <- results[[i]][[9]]
    # 
    # # ar test
    # arBagSize[i,] <- results[[i]][[10]]
    # arBagSizeCPS[i,] <- results[[i]][[11]]
    # arBagSizeStack[i,] <- results[[i]][[12]]
    
    # rmse - main results
    final.results[i,] <- results[[i]][[13]]
}


# write files
setwd(save.folder)
write.csv(ssTuneResults, ssTuneName, row.names = FALSE)
write.csv(ssTuneResultsCPS, paste0(ssTuneName, "_CPS"), row.names = FALSE)
write.csv(ssTuneResultsStack, paste0(ssTuneName, "_Stack"), row.names = FALSE)

write.csv(ssBagSize, paste0(ssBagSizeTest, "_", "Avg"), row.names = FALSE)
write.csv(ssBagSizeCPS, paste0(ssBagSizeTest, "_", "CPS"), row.names = FALSE)
write.csv(ssBagSizeStack, paste0(ssBagSizeTest, "_", "Stack"), row.names = FALSE)

# write.csv(arBagSizeTune, arTuneName, row.names = FALSE)
# write.csv(arBagSizeCPSTune, paste0(arTuneName, "_CPS"), row.names = FALSE)
# write.csv(arBagSizeStackTune, paste0(arTuneName, "_Stack"), row.names = FALSE)
# 
# write.csv(arBagSize, paste0(arBagSizeTest, "_", "Avg"), row.names = FALSE)
# write.csv(arBagSizeCPS, paste0(arBagSizeTest, "_", "CPS"), row.names = FALSE)
# write.csv(arBagSizeStack, paste0(arBagSizeTest, "_", "Stack"), row.names = FALSE)

# save results
colnames(final.results) <- c("mrg", "sse", "sse_stack", "sse_cps", "ss", "ss_stack", "ss_cps", 
                        "ar", "ar_stack", "ar_cps")

print("dim results: ")
print(dim(final.results))
setwd(save.folder)
write.csv(final.results, fName, row.names = FALSE)

