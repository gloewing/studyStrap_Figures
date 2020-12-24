library(doParallel)
library(CCA)
library(MatrixCorrelation)
library(caret)
library(glmnet)
library(foreach)

source("Study Strap Functions.R")
source("SimFn.R")
source("ss.ar.caret.mSSL.R")
source("studyStrap.predict.caret.mSSL.R")
source("fatTrim.R")
source("SSE.caret.mSSL.R")
source("merge.caret.mSSL.R")
source("ss.caret.mSSL.R")

# sims <- cbind(160:183, expand.grid(c(0.05, 0.25, 1, 3), c(0.05, 5, 20), c(0,4)))
# sims2 <- cbind( 184:191, expand.grid(c(0.05,0.25,1,3),c(10,15), c(0)))
# sims3 <- cbind(192:203,  expand.grid(c(0.05, 0.25, 1, 3), c(0.05, 5, 20), c(8)))
# colnames(sims) <- c("simNum", "betaVar", "XVar", "Clusts")
# colnames(sims2) <- c("simNum", "betaVar", "XVar", "Clusts")
# colnames(sims3) <- c("simNum", "betaVar", "XVar", "Clusts")
# sims <- rbind(sims,sims2, sims3)

save.folder <- "/n/home12/gloewinger/sims5"
load.folder <- "~/Desktop/Research"

totalSims <- 100

# sim parameters
trainStudies <- 1:16
num.trainStudy <- length(trainStudies)
totalStudies <- 24 # need number that is divisible by 2, 4, 8 for cluster purposes and structure of 

simNum <- 202
clusts <- 8 # totalStudies # 4
betaVar <- 1
covVar <- 20
# simulation code but below we only use 16 training studies and 1 test study

# Study Strap and Model fitting parameters
bagSzs <- c(1:3, seq(4,18, by = 2), seq(20,100, by = 10), 200, 500, 1000)
bagSize.vec<- bagSzs
test_study <- 17 # arbitrarily set to 17th study but any of the non-training studies is fine (17-24 are all random test studies)
numCovs <- 20 # number of covaraites
numBootItrs <- 500 # number of bootstrap replicates to estimate variance for similarity function
ssStrapsTune <- 150 # number of study straps for ss Tuning
ssStrapsTest <- 500 # number of study straps for Test SS
ssStrapsBagSzTest <- 250 # number of study straps for Test SS Bag Size

tune.grid <- sort( c(0.001, 0.01, 0.1, seq(0.0001, 5, length = 40) ) )
totalPaths <- 5 # paths for AR Test
convgLim <- 100000 # convergence limit for AR
convgLimBagSize <- 10000 #AR converge limit for AR bag size tune
convgLimBagSizeTest <- 10000 #AR converge limit for AR bag size test

tune.grid <- as.data.frame(tune.grid) # tuning parameters to consider
tune.grid <- cbind(1,tune.grid)
colnames(tune.grid) <- c("alpha", "lambda")
# parallelize

######### CHANGE NUMBER OF THREADS BEFORE RUNNING REAL THING !!!!!!!!!!!!!!!!!!!!!!!
logfile <- paste0("outputFile", simNum,".txt")
writeLines(c(""), file(logfile,'w'))

num.threads <- as.integer( ceiling( totalSims / 3 ) )# round up
threads <- makeCluster(num.threads, outfile=logfile)
registerDoParallel(threads)

setwd(save.folder)

getDoParWorkers()
timeStart <- Sys.time()

results <- foreach(iterNum = 1:totalSims, .combine = list, .multicombine = TRUE) %dopar%{
        
        print(paste0("start: ", iterNum))

        resList <- vector(length = 13,  "list")
        
        # save results
        setwd(save.folder)
        final.results <- matrix(ncol = 10, nrow = totalSims)
        colnames(final.results) <- c("mrg", "sse", "sse_stack", "sse_cps", "ss", "ss_cps", 
                                     "ss_stack", "ar", "ar_stack", "ar_cps")
        
        rmse.vec <- c() # store results
        
        # file names
        filename <- paste0("Sim ", simNum, "_", iterNum, "_Combined")
        baseName <- paste0("Sim ", simNum, "_Combined")
        fileNm <- paste0("LASSO_Sims5_Results_", simNum, "_", iterNum, "_Combined")
        
        arTuneName <- paste0("LASSO AR Tune_Sims5_", baseName, "_",iterNum, "_")
        ssTuneName <- paste0("LASSO SS Tune_Sims5_", baseName, "_",iterNum, "_")
        arBagSizeTest <- paste0("LASSO AR BagSize_Test_Sims5_", baseName, "_",iterNum, "_")
        ssBagSizeTest <- paste0("LASSO SS BagSize_Test_Sims5_", baseName, "_",iterNum, "_")
        
        
        seedSet <- iterNum # ensures no repeats
        set.seed(seedSet)
        
        
        # simulate data
        full <- as.data.frame( multiStudySim(sim.num = simNum,
                                  iter = iterNum,
                                  sampSize = 400,
                                  covariate.var = covVar, # scales the variance of the MVNormal that generates the true means of the covaraites
                                  beta.var = betaVar, # variance of betas
                                  clust.num.mu = clusts, # number of clusters of means (if this is set to num.studies, than each true mean of covariate is different)
                                  clust.num.beta = clusts, # number of clusters of betas (if this is set to num.studies, than each true mean of covariate is different)
                                  num.covariates = numCovs,
                                  zero_covs = 11:20, # indices of covariates which are 0
                                  num.studies = totalStudies,
                                  beta.mean.range = 10, # true means of hyperdistribution of beta are drawn from a unif(-beta.mean.range, beta.mean.range)
                                  perturb = betaVar * 0.1 / 2, # perturb = 0 means all clusters are identical. otherwise perturnance of betas are elementwise drawn from a unif(-perturb, perturb)
                                  SB = 1) )

        # delete non-training/test studies
        # --- necessary based on simulation code
        allStudies <- c(trainStudies, test_study)
        full <- as.data.frame(full[is.element(full$Study, allStudies), ])

        # test study
        test.elec <- as.data.frame(full[full$Study == test_study, -1]) # set test electrode
        full <- as.data.frame(full[is.element(full$Study, trainStudies), ]) # include only training studies in full, remove rownames

        varMat <- matrix(ncol = numCovs + 1, nrow = num.trainStudy)

        ### Fit each model tuned individually to come up with similarity metric
        betaMat <- matrix(ncol = numCovs + 1, nrow = num.trainStudy)
        for(y in 1:num.trainStudy){
            studyItr <- as.data.frame(full[full$Study == trainStudies[y], -1]) # set test electrode

            control <- caret::trainControl(method="cv", number=10)
            tuneModel <- caret::train(Y~., data = studyItr, method = "glmnet",
                                      trControl = control, tuneGrid = tune.grid, metric = "RMSE")
            betaMat[y,]<- as.numeric( coef(tuneModel$finalModel, tuneModel$finalModel$lambdaOpt) )
            optLambda <- data.frame(lambda = tuneModel$finalModel$lambdaOpt, alpha = 1)

            # variance estimation of parameter estimates
            betaBoot <- matrix(ncol = numCovs + 1, nrow = numBootItrs)
            n <- nrow(studyItr)


            for(i in 1:numBootItrs){
                # fit bootstrapped model
                bootRows <- sample.int(n, replace = TRUE) # bootstrap
                tuneModel <- caret::train(Y~., data = studyItr[bootRows,], method = "glmnet",
                                          tuneGrid = optLambda)
                # add coefficient estimates
                betaBoot[i,] <- as.numeric( coef(tuneModel$finalModel, tuneModel$finalModel$lambdaOpt) ) # add bootstrap coef estimates
            }
            varMat[y,] <- apply(betaBoot,2, var) # bootstrap variance of each coefficient estimates
            rm(studyItr)
            }

        varMat <- ifelse(varMat == 0, 10e-9, varMat) # if any elements are 0, set them to a very small number to avoid numerical instability
        wMat <- 1 / varMat[,-1] # we will make weights based on inverse of variance estimate
        wMat <- apply(wMat,2, function(x) x / sum(x)) # turn these into weights by making them sum to one

        betaMat <- betaMat[,-1] * wMat # weight estimates by the inverse of variance

        betaWeights <- abs( colMeans(betaMat) ) # remove intercept
        betaWeights <- betaWeights / sum( betaWeights ) # standardize
        rm(betaMat)
        rm(varMat)
        rm(wMat)
        print(paste0("betaWeights: ", betaWeights))
        # mod$finalModel$nulldev

        # similarity measure
        simFn <- function(x1, x2, weights = betaWeights){
            return( (1 / sum( weights * ( colMeans(x1) - colMeans(x2) )^2  )   ) )
        }

        #######################
        # Tune Lasso parameter
        #######################
        # hold one out study out CV indices of studies
        indxList <- vector("list", length = num.trainStudy)
        HOOList <- vector("list", length = num.trainStudy)

        for(study in 1:num.trainStudy){
            indxList[[study]] <- which(full$Study != trainStudies[study]) # indices of studies to train on
            HOOList[[study]] <- which(full$Study == trainStudies[study]) # indices of study to hold out

        }

        control <- caret::trainControl(method="cv", number=num.trainStudy, index = indxList, indexOut = HOOList)
        tuneModel <- caret::train(Y~., data=full[,-1], method="glmnet",
                                  trControl=control, tuneGrid = tune.grid, metric="RMSE")

        testLambda <- tuneModel$bestTune
        testLambda <- as.data.frame(testLambda)
        # use this lambda for all non-tuning training
        ########################################

        ##############
        # Merge
        ##############
        print(paste0("Merge: ", iterNum))
        mrg.mod <- merged(formula = Y ~., data = full, model = FALSE,
                          ssl.method = list("glmnet"),
                          ssl.tuneGrid = list( testLambda ))
        preds <- studyStrap.predict(mrg.mod, test.elec[,-1])
        rm(mrg.mod)
        rmse <- sqrt(mean( (preds[,1] - test.elec$Y)^2   ))
        rmse.vec <- c(rmse.vec, rmse)
        ##########################################################################################

        ##############
        # SSE
        ##############
        print(paste0("SSE: ", iterNum))
        sse.mod <- sse(formula = Y ~., data = full, target.study = test.elec[,-1], # remove Y from test.elec
                       sim.mets = FALSE, model = FALSE,
                       ssl.method = list("glmnet"),
                       ssl.tuneGrid = list( testLambda ), customFNs = list(simFn) )

        preds <- studyStrap.predict(sse.mod, test.elec[,-1])
        rm(sse.mod)
        rmse.avg <- sqrt(mean( (preds[,1] - test.elec$Y)^2  ))
        rmse.stack <- sqrt(mean( (preds[,2] - test.elec$Y)^2  ))
        rmse.cps <- sqrt(mean( (preds[, ncol(preds)] - test.elec$Y)^2  ))  # last column is custom function
        rmse.vec <- c(rmse.vec, rmse.avg, rmse.stack, rmse.cps)


        ###########
        # SS Tune
        ###########

        bagSize.vec <- bagSzs
        print(paste0("SS Tune: ", iterNum))

        ssTuneResults <- matrix(NA, ncol = length(bagSize.vec), nrow = length(trainStudies))
        ssTuneResultsCPS <- matrix(NA, ncol = length(bagSize.vec), nrow = length(trainStudies))
        ssTuneResultsStack <- matrix(NA, ncol = length(bagSize.vec), nrow = length(trainStudies))
        # Load Data
        original.Train <- trainStudies

        #iterate through bag sizes

        for(bagSize in bagSzs){

            # iterate through held out studies

            for(y in trainStudies){

                original.Train <- trainStudies
                print(paste0("bagSize: ", bagSize))

                original.Train <- original.Train[original.Train != y] # remove validation set from training study vector

                #######################
                # Tune Lasso parameter
                #######################

                # hold one out study out CV indices of studies
                indxList <- vector("list", length = length(original.Train))
                HOOList <- vector("list", length = length(original.Train))

                for(study in 1:length(original.Train)){
                    indxList[[study]] <- which(full$Study != c(original.Train[study], y) )  # indices of studies to train on
                    HOOList[[study]] <- which(full$Study == original.Train[study]) # indices of study to hold out

                }

                # Tune lambda

                control <- caret::trainControl(method="cv", index = indxList, indexOut = HOOList)
                tuneModel <- caret::train(Y~., data=full[,-1], method="glmnet",
                                          trControl=control, tuneGrid = tune.grid, metric="RMSE")

                best.lambda <- tuneModel$bestTune
                best.lambda <- as.data.frame(best.lambda)
                
                fullIndx <- which(full$Study != y)
                testIndx <- which(full$Study == y)

                ssMod <-  ss(data = full[fullIndx,], formula = Y ~.,
                             target.study = full[testIndx,-c(1,2)],
                             bag.size = bagSize,
                             straps = ssStrapsTune,
                             stack = "standard",
                             ssl.method = list("glmnet"),
                             ssl.tuneGrid = list( testLambda ),
                             sim.mets = FALSE,
                             stack.standardize = FALSE, customFNs = list(simFn))

                preds <- studyStrap.predict(ssMod, full[testIndx,-1])
                preds <- preds[,c(1,2,ncol(preds))] # average, stacking, CPSW
                rm(ssMod) # delete model object to save memory
                rmse <- apply( (test.elec$Y - preds)^2, 2, FUN = function(x) sqrt(mean(x)))
                print(paste0("final results, study: ", y, "_bagSize: ", bagSize, "_RMSE_", rmse))
                indx <- which(trainStudies == y) # index for results storage
                bag.indx <- which(bagSize.vec == bagSize)

                ssTuneResults[indx, bag.indx] <- rmse[1]
                ssTuneResultsCPS[indx, bag.indx] <- rmse[3]
                ssTuneResultsStack[indx, bag.indx] <- rmse[2]

            }
        }


            # save for output
            resList[[1]] <- colMeans(ssTuneResults)
            resList[[2]] <- colMeans(ssTuneResultsCPS)
            resList[[3]] <- colMeans(ssTuneResultsStack)
            
        ssRMSE <- colMeans(ssTuneResults) # average across held out datasets
        ssBagSize <- bagSize.vec[which.min(ssRMSE)] # optimal bag size
        
        # save memory
        rm(ssTuneResults)
        rm(ssTuneResultsCPS)
        rm(ssTuneResultsStack)
        ############################################

        ##############
        # Test ss
        ##############
        print(paste0("Test SS: ", iterNum))

        ssMod <-  ss(data = full, formula = Y ~.,
                     target.study = test.elec[,-1],
                     bag.size = ssBagSize,
                     straps = ssStrapsTest,
                     stack = "standard",
                     ssl.method = list("glmnet"),
                     ssl.tuneGrid = list( testLambda ),
                     sim.mets = FALSE,
                     stack.standardize = FALSE,
                     customFNs = list(simFn))

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

            ssMod <-  ss(data = full, formula = Y ~.,
                         target.study = test.elec[,-1],
                         bag.size = bagSize,
                         straps = ssStrapsBagSzTest,
                         stack = "standard",
                         ssl.method = list("glmnet"),
                         ssl.tuneGrid = list( testLambda ),
                         sim.mets = FALSE,
                         stack.standardize = FALSE,
                         customFNs = list(simFn))

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

        ###########
        # AR Tune
        ###########
        bagSize.vec <- bagSzs
        print(paste0("AR Tune: ", iterNum))

        arTuneResults <- matrix(NA, ncol = length(bagSize.vec), nrow = length(trainStudies))
        arTuneResultsCPS <- matrix(NA, ncol = length(bagSize.vec), nrow = length(trainStudies))
        arTuneResultsStack <- matrix(NA, ncol = length(bagSize.vec), nrow = length(trainStudies))
        # Load Data
        original.Train <- trainStudies

        #iterate through bag sizes

        for(bagSize in bagSzs){

            # iterate through held out studies

            for(y in trainStudies){

                original.Train <- trainStudies
                print(paste0("bagSize: ", bagSize))

                original.Train <- original.Train[original.Train != y] # remove validation set from training study vector

                #######################
                # Tune Lasso parameter
                #######################

                # hold one out study out CV indices of studies
                indxList <- vector("list", length = length(original.Train))
                HOOList <- vector("list", length = length(original.Train))

                for(study in 1:length(original.Train)){
                    indxList[[study]] <- which(full$Study != c(original.Train[study], y) )  # indices of studies to train on
                    HOOList[[study]] <- which(full$Study == original.Train[study]) # indices of study to hold out
                    
                }

                control <- caret::trainControl(method="cv", number=num.trainStudy, index = indxList, indexOut = HOOList)
                tuneModel <- caret::train(Y~., data=full[,-1], method="glmnet",
                                          trControl=control, tuneGrid = tune.grid, metric="RMSE")

                best.lambda <- tuneModel$bestTune
                best.lambda <- as.data.frame(best.lambda)
                
                fullIndx <- which(full$Study != y)
                testIndx <- which(full$Study == y)

                arMod <-  cmss(data = full[fullIndx,], formula = Y ~.,
                               target.study = full[testIndx, -c(1,2)],
                               converge.lim = convgLimBagSize,#10000,
                               bag.size = bagSize,
                               max.straps = 150,
                               paths = 3,
                               stack = "standard",
                               ssl.method = list("glmnet"),
                               ssl.tuneGrid = list( best.lambda ),
                               sim.mets = FALSE, #model = FALSE,
                               meanSampling = FALSE,
                               stack.standardize = FALSE,
                               sim.fn = simFn,
                               customFNs = list(simFn))

                preds <- studyStrap.predict(arMod, full[testIndx, -c(1,2)])
                preds <- preds[,c(1, 2, ncol(preds))] # average, stacking, CPSW
                rm(arMod) # delete model object to save memory
                rmse <- apply( (test.elec$Y - preds)^2, 2, FUN = function(x) sqrt(mean(x)))
                print(paste0("final results, study: ", y, "_bagSize: ", bagSize, "_RMSE_", rmse))
                indx <- which(trainStudies == y) # index for results storage
                bag.indx <- which(bagSize.vec == bagSize)

                arTuneResults[indx, bag.indx] <- rmse[1]
                arTuneResultsCPS[indx, bag.indx] <- rmse[3]
                arTuneResultsStack[indx, bag.indx] <- rmse[2]

            }
        }

        # fill in corresponding row with colMeans
        resList[[7]] <- colMeans(arTuneResults)
        resList[[8]] <- colMeans(arTuneResultsCPS)
        resList[[9]] <-  colMeans(arTuneResultsStack)

        arRMSE <- colMeans(arTuneResults) # average across held out datasets
        ARbagSize <- bagSize.vec[which.min(arRMSE)] # optimal bag size
        
        rm(arTuneResults)
        rm(arTuneResultsCPS)
        rm(arTuneResultsStack)
        ############################################
        ##############
        # Test AR
        ##############
        print(paste0("Test AR: ", iterNum))

        if( ARbagSize == 1){
            # if it is only 1 then there are very few possible combinations--no need to waste time
            convgLim <- num.trainStudy * 2
        }

        arMod <-  cmss(data = full, formula = Y ~.,
                       target.study = test.elec[,-1],
                       converge.lim = convgLim,
                       bag.size = ARbagSize,
                       max.straps = 200,
                       paths = totalPaths,
                       stack = "standard",
                       ssl.method = list("glmnet"),
                       ssl.tuneGrid = list( testLambda ),
                       sim.mets = FALSE, model = FALSE,
                       meanSampling = FALSE,
                       stack.standardize = FALSE,
                       sim.fn = simFn,
                       customFNs = list(simFn))

        preds <- studyStrap.predict(arMod, test.elec[,-1])
        preds <- preds[,c(1, 2, ncol(preds))] # average, stacking, CPSW
        rm(arMod) # delete model object to save memory
        rmse.vec <- c(rmse.vec, apply( (test.elec$Y - preds)^2, 2, FUN = function(x) sqrt(mean(x))) )

        #####################################################
        # Bag Sizes Test AR
        ####################################################
        print(paste0("Test AR Bag Size: ", iterNum))
        # Load Data
        original.Train <- trainStudies

                # cycle through bagsizes
        arBagSize <- matrix(ncol = length(bagSzs), nrow = totalSims)
        arBagSizeCPS <- matrix(ncol = length(bagSzs), nrow = totalSims)
        arBagSizeStack <- matrix(ncol = length(bagSzs), nrow = totalSims)

        for(bagSize in bagSzs){

            arMod <-  cmss(data = full, formula = Y ~.,
                           target.study = test.elec[,-1],
                           converge.lim = convgLimBagSizeTest,
                           bag.size = bagSize,
                           max.straps = 200,
                           paths = 3,
                           stack = "standard",
                           ssl.method = list("glmnet"),
                           ssl.tuneGrid = list( testLambda ),
                           sim.mets = FALSE, #model = FALSE,
                           meanSampling = FALSE,
                           stack.standardize = FALSE,
                           sim.fn = simFn,
                           customFNs = list(simFn))

            preds <- studyStrap.predict(arMod, test.elec[,-1])
            preds <- preds[,c(1, 2, ncol(preds))] # average, stacking, CPSW
            rm(arMod) # delete model object to save memory
            rmse <- apply( (test.elec$Y - preds)^2, 2, FUN = function(x) sqrt(mean(x)))
            print(paste0("final results, study: ", y, "_bagSize: ", bagSize, "_RMSE_", rmse))
            bagSzIndx <- which(bagSzs == bagSize)

            arBagSize[iterNum, bagSzIndx] <- rmse[1]
            arBagSizeCPS[iterNum, bagSzIndx] <- rmse[3]
            arBagSizeStack[iterNum, bagSzIndx] <- rmse[2]
        }

        resList[[10]] <- arBagSize[iterNum,]
        resList[[11]] <- arBagSizeCPS[iterNum,]
        resList[[12]] <- arBagSizeStack[iterNum,]
        
        rm(arBagSize)
        rm(arBagSizeCPS)
        rm(arBagSizeStack)

        print(paste0("ForeEAch Iteration Complete: ", simNum, "_IterNum: ", iterNum ) )
        print(rmse.vec)

        resList[[13]] <-  rmse.vec
        
        return(resList)
        }        

print("Complete, length of results")
print(length(results))

timeEnd <- Sys.time()
print(difftime(timeEnd, timeStart, units='mins'))


### names
filename <- paste0("Sim ", simNum,  "_Combined")
baseName <- paste0("Sim ", simNum, "_Combined")
fileNm <- paste0("LASSO_Sims5_Results_", simNum, "_Combined")

arTuneName <- paste0("LASSO AR Tune_Sims5_", baseName)
ssTuneName <- paste0("LASSO SS Tune_Sims5_", baseName)
arBagSizeTest <- paste0("LASSO AR BagSize_Test_Sims5_", simNum)
ssBagSizeTest <- paste0("LASSO SS BagSize_Test_Sims5_", simNum)

### matrices

# ss tune
ssTuneResults <- matrix(NA, ncol = length(bagSize.vec), nrow = totalSims)
ssTuneResultsCPS <- matrix(NA, ncol = length(bagSize.vec), nrow = totalSims)
ssTuneResultsStack <- matrix(NA, ncol = length(bagSize.vec), nrow = totalSims)

# ss test
ssBagSize <- matrix(ncol = length(bagSzs), nrow = totalSims)
ssBagSizeCPS <- matrix(ncol = length(bagSzs), nrow = totalSims)
ssBagSizeStack <- matrix(ncol = length(bagSzs), nrow = totalSims)

# ar tune
arBagSizeTune <- matrix(ncol = length(bagSzs), nrow = totalSims)
arBagSizeCPSTune <- matrix(ncol = length(bagSzs), nrow = totalSims)
arBagSizeStackTune <- matrix(ncol = length(bagSzs), nrow = totalSims)

# ar test
arBagSize <- matrix(ncol = length(bagSzs), nrow = totalSims)
arBagSizeCPS <- matrix(ncol = length(bagSzs), nrow = totalSims)
arBagSizeStack <- matrix(ncol = length(bagSzs), nrow = totalSims)

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
    
    # ar tune
    arBagSizeTune[i,] <- results[[i]][[7]]
    arBagSizeCPSTune[i,] <- results[[i]][[8]]
    arBagSizeStackTune[i,] <- results[[i]][[9]]
    
    # ar test
    arBagSize[i,] <- results[[i]][[10]]
    arBagSizeCPS[i,] <- results[[i]][[11]]
    arBagSizeStack[i,] <- results[[i]][[12]]
    
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

write.csv(arBagSizeTune, arTuneName, row.names = FALSE)
write.csv(arBagSizeCPSTune, paste0(arTuneName, "_CPS"), row.names = FALSE)
write.csv(arBagSizeStackTune, paste0(arTuneName, "_Stack"), row.names = FALSE)

write.csv(arBagSize, paste0(arBagSizeTest, "_", "Avg"), row.names = FALSE)
write.csv(arBagSizeCPS, paste0(arBagSizeTest, "_", "CPS"), row.names = FALSE)
write.csv(arBagSizeStack, paste0(arBagSizeTest, "_", "Stack"), row.names = FALSE)

# save results
colnames(final.results) <- c("mrg", "sse", "sse_stack", "sse_cps", "ss", "ss_stack", "ss_cps", 
                        "ar", "ar_stack", "ar_cps")

print("dim results: ")
print(dim(final.results))
setwd(save.folder)
write.csv(final.results, paste0("LASSO_Sims5_Results_", simNum), row.names = FALSE)

