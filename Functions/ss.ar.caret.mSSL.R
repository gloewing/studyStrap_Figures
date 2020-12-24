#' Covariate-Matched Study Strap for Multi-Study Learning: Fits accept/reject algorithm based on covariate similarity measure
#'
#' @param formula Model formula
#' @param data A dataframe with all the studies has the following columns in this order: "Study", "Y", "V1", ...., "Vp"
#' @param target.study Dataframe of the design matrix (just covariates) of study one aims to make predictions on
#' @param sim.fn Optional function to be used as similarity measure for accept/reject step. Default function is: |cor( bar{x}^{(r)}|,~ bar{x}_{target} ) |
#' @param converge.lim Integer indicating the number of consecutive rejected study straps to reach convergence criteria.
#' @param bag.size Integer indicating the bag size tuning parameter.
#' @param max.straps Integer indicating the maximum number of accepted straps that can be fit across all paths before the algorithm stops accepting new study straps.
#' @param paths Integer indicating the number of paths (an accept/reject path is all of the models accepted before reaching one convergence criteria).
#' @param stack String determining whether stacking matrix made on training studies "standard" or on the accepted study straps "ss." Default: "standard."
#' @param sim.covs Is a vector of names of covariates or the column numbers of the covariates to be used for the similarity measure. Default is to use all covariates.
#' @param ssl.method A list of strings indicating which modeling methods to use.
#' @param ssl.tuneGrid A list of the tuning parameters in the format of the caret package. Each element must be a dataframe (as required by caret). If no tuning parameters are required then NA is indicated.
#' @param sim.mets Boolean indicating whether to calculate default covariate profile similarity measures.
#' @param model Indicates whether to attach training data to model object.
#' @param meanSampling = FALSE Boolean determining whether to use mean covariates for similarity measure. This can be much quicker.
#' @param customFNs Optional list of functions that can be used to add custom covaraite profile similarity measures.
#' @param stack.standardize Boolean determining whether stacking weights are standardized to sum to 1. Default is FALSE
#' @return A model object of studyStrap class "ss" that can be used to make predictions.
#' @examples
#' ##########################
#' ##### Simulate Data ######
#' ##########################
#'
#' set.seed(1)
#' # create half of training dataset from 1 distribution
#' X1 <- matrix(rnorm(2000), ncol = 2) # design matrix - 2 covariates
#' B1 <- c(5, 10, 15) # true beta coefficients
#' y1 <- cbind(1, X1) %*% B1
#'
#' # create 2nd half of training dataset from another distribution
#' X2 <- matrix(rnorm(2000, 1,2), ncol = 2) # design matrix - 2 covariates
#' B2 <- c(10, 5, 0) # true beta coefficients
#' y2 <- cbind(1, X2) %*% B2
#'
#' X <- rbind(X1, X2)
#' y <- c(y1, y2)
#'
#' study <- sample.int(10, 2000, replace = TRUE) # 10 studies
#' data <- data.frame( Study = study, Y = y, V1 = X[,1], V2 = X[,2] )
#'
#' # create target study design matrix for covariate profile similarity weighting and
#' # accept/reject algorithm (covaraite-matched study strap)
#' target <- matrix(rnorm(1000, 3, 5), ncol = 2) # design matrix
#' colnames(target) <- c("V1", "V2")
#'
#' ##########################
#' ##### Model Fitting #####
#' ##########################
#'
#' # Fit model with 1 Single-Study Learner (SSL): PCA Regression
#' arMod1 <-  cmss(formula = Y ~.,
#'                data = data,
#'                target.study = target,
#'                converge.lim = 10,
#'                bag.size = length(unique(data$Study)),
#'                max.straps = 50,
#'                paths = 2,
#'                ssl.method = list("pcr"),
#'                ssl.tuneGrid = list(data.frame("ncomp" = 2))
#'                )
#'
#'# Fit model with 2 SSLs: Linear Regression and PCA Regression
#' \donttest{arMod2 <-  cmss(formula = Y ~.,
#'                data = data,
#'                target.study = target,
#'                converge.lim = 20,
#'                bag.size = length(unique(data$Study)),
#'                max.straps = 50,
#'                paths = 2,
#'                ssl.method = list("lm", "pcr"),
#'                ssl.tuneGrid = list(NA, data.frame("ncomp" = 2))
#'                )}
#'
#'
#'
#' # Fit model with custom similarity function for
#' # accept/reject step and 2 custom function for Covariate
#' # Profile Similarity weights
#'
#'# custom function for CPS
#'
#' fn1 <- function(x1,x2){
#' return( abs( cor( colMeans(x1), colMeans(x2) )) )
#' }
#'
#' fn2 <- function(x1,x2){
#' return( sum ( ( colMeans(x1) - colMeans(x2) )^2 ) )
#' }
#'
#' \donttest{arMod3 <-  cmss(formula = Y ~.,
#'                data = data,
#'                target.study = target,
#'                sim.fn = fn1,
#'                customFNs = list(fn1, fn2),
#'                converge.lim = 50,
#'                bag.size = length(unique(data$Study)),
#'                max.straps = 50,
#'                paths = 2,
#'                ssl.method = list("lm", "pcr"),
#'                ssl.tuneGrid = list(NA, data.frame("ncomp" = 2))
#'                )}
#'
#' #########################
#' #####  Predictions ######
#' #########################
#'
#' preds <- studyStrap.predict(arMod1, target)
#' @import caret
#' @import tidyverse
#' @import nnls
#' @import dplyr
#' @importFrom stats coef cor cov predict
#' @export

cmss <- function(formula = Y ~.,
                data,
                target.study,
                sim.fn = NA,
                converge.lim = 50000,
                bag.size = length(unique(data$Study)),
                max.straps = 150,
                paths = 5,
                stack = "standard",
                stackTune = TRUE,
                stackNNLS = TRUE,
                stackIntercept = TRUE,
                sim.covs = NA,
                ssl.method = list("lm"),
                ssl.tuneGrid = list(c()),
                sim.mets = TRUE,
                model = FALSE,
                meanSampling = FALSE,
                customFNs = list(),
                stack.standardize = FALSE,
                stack.lambda = sort( unique( c(0, 0.0001, 0.001, 0.01, 
                                                             exp(-seq(0,5, length = 50)),
                                                             seq(5,20, by = 5),
                                                             seq(30,100, by = 10) ) ) ),
                tune = NULL,
                nfolds = 5,
                tune.grid = list(c()) 
                ){
    # Takes data that has the following columns in this order: "Study", "Y", "V_1", ...., "V_p"
        # must be ordered in that way and must include "Study" and "Y", the names of the covariates can vary
    # data is a dataset that includes data$Study (integers indicating study number)
    # should not include rows of target study
    # testStudy is the study number of the target
    # paths is the number of random paths to sample
    # converge.lim is the number of random samples to take before ending the path
    # should include target study
    # formula is the command for the model (not including first two columns
    # -- these are automatically removed during modelling step)
    # stack is the type of stacking procedure to be used - "standard" uses the full dataframe
    # "ss" - uses all the rows from all study straps
    # sim.covs is a vector of the column numbers or a vector of the names of the covariates
        # to be used in the similarity measure
        # default is to use all covariates
    # target.study should be ONLY the design matrix of the target study (no outcome (Y), no Study labels)
    # ssl.method determines the method to use (examples: "lm", "pcr")
    # ssl.tuneGrid is a dataframe with columns that contains the tuning parameters. column names must be named appropriately
        # for more info on models see: https://topepo.github.io/caret/available-models.html
        # and for tuning parameters: http://topepo.github.io/caret/model-training-and-tuning.html#model-training-and-parameter-tuning
    # sim.metrics determines whether or not to calculate similarity measures---slower
  # model determines whether raw data is attached
  # meanSamp indicates whether to sample the covariate means for the similarity measure -- faster but less precise
  # customFNs is a list where each element is a custom function for CPS weighting (not used for accept/reject step)
  # stack.standardize indicates whether stacking coefs are standardized
    # source("fatTrim.R") # to reduce model size
    # source("similarity.metrics.R") # similarity metrics

    # library(dplyr)
    # library(caret)
    # library(nnls)
    library(glmnet)
    
    # whether to do nonnegative least squares for stacking
    if(stackNNLS){
        wLim <- 0
    }else{
        wLim <- -Inf
    }
    
    # if do not tune, then use no ridge penalty
    if(!stackTune) stackParam <- 0

    num.SSLs <- length(ssl.method) # number of SSLs per study strap

    data <- as.data.frame(data)
    original.studies <- unique(data$Study)
    study.clmn <- which(colnames(data) == "Study") # column of Study label

    # check outcome name from formula in case it is not "Y"
    Yname <- toString(formula[[2]])

    if( anyNA(sim.covs) ){
      # if not specified use all the covariates (assuming the first two columns are Study and Y)
      sim.covs <- seq(1, ncol(data) )[-c(1,2)]
    }

    target.sim.covs <- sim.covs - 2 # the target study does not have Study and outcome (Y) columns
    # so shift indices over by 2

    mean.target.covs <- colMeans(target.study[,target.sim.covs]) # mean covariates of target study

    Study.Code <- NULL # include for namespace

    # rename studies as integers in order from 1:length(unique(data$Study))
    study.mat <- tibble::tibble(Study = unique(data$Study), Study.Code = 1:length(unique(data$Study)))
    data <- data %>%
      left_join(study.mat, by = "Study") %>%
      mutate(Study = Study.Code) %>%
      dplyr::select(-Study.Code)
    rm(study.mat)

    studies <- unique(data$Study)
    row.list <- list(length = max.straps) # list of the rows that were included in study strap

    # determine number of rows in each study and produce covariate mean matrix
    sampSizes <- as.vector( table(data$Study) )
    covMeans <- matrix(nrow = length(studies), ncol = length(sim.covs) ) # matrix of study means
    for (i in studies){
      # sampSizes <- c(sampSizes, length( data$Study[data$Study == unique(data$Study)[i] ]  ) )
      covMeans[i,] <- colMeans( data[data$Study == i, sim.covs ] )
    }

    model.list <- vector("list", length = num.SSLs)

    for(mod in 1:num.SSLs){
      # each element in list is a list of models for each study strap
      model.list[[mod]] <- vector("list", length = max.straps)
    }

    ss.obj <- list(models = model.list, data = c(), sim.mat <- c(),
                   strapRows <- list(),
                   dataInfo = list(studyNames = original.studies, sampleSizes = sampSizes),
                   modelInfo = list(sampling = "ar", numStraps = max.straps, SSL = ssl.method,
                                    ssl.tuneGrid = c(), numPaths = paths,
                                     convg.vec = c(), convgCritera = converge.lim,
                                    meanSamp = meanSampling, stack.type = stack,
                                    custFNs = customFNs, bagSize = bag.size),
                   stack.coefs <- c(), simMat = c())
    class(ss.obj) <- "ss"
    convgVec <- c() # each element is the number of candidate straps before acceptance
    rm(model.list)
    # similarity matrix contains the average covariates that are specified
    similarity.matrix <- matrix(ncol = 23, nrow = max.straps)  # similarity metrics

    if(length(customFNs) > 0){
      # if there are any custom functions for CPS weighting
      custom.matrix <- matrix(ncol = length(customFNs), nrow = max.straps) # custom similarity metrics
      colnames(custom.matrix) <- paste0("customFn_", 1:length(customFNs))
      }

    if( !is.function(sim.fn) ){
      # use default similarity function if none provided
      sim.fn <- function(dat1, dat2){
        return( abs(cor(colMeans(dat1), colMeans(dat2))) )
      }
    }
#######-----------------------------------------------------

    rpts <- 1 # number of repeated paths
    z <- 1
    sim.metric <- 0 # arbitrarily small
    convgVec <- c()
    while(rpts <= paths && z <= max.straps){

      message(paste0("Study Strap: ", z, ", Path: ", rpts))

      counter <- 0 # reset counter at each new sample

      ########################
      ######## Boostraping
      ########################

      sim <- 0 # this is an indicator telling the while loop to break (if = 1)
      repeat.break <- 0
      total.straps <- 0 # this tracks total number of candidate straps
      while (sim == 0 && repeat.break == 0){ #while loop for restrapping: breaks for either acceptance (sim) or accpetance counter limit (converge.lim)

        if (counter >= converge.lim){
          counter <- 0
          sim.metric <- 0 # arbitrarily small
          rpts <- rpts + 1
          repeat.break <- 1
        }

        ########################
        ######## Study Bag
        ########################
        bag <- sample(studies, bag.size, replace = TRUE)
        strap.table <- as.data.frame(table(bag)) # put in data table format

        #####################################################################
        # Study Strap Sampling
        #####################################################################
        if( meanSampling == TRUE){
          # use means to determine similarity meausre to speed up AR sampling
          current.sim.metric <- sim.fn( target.study[,target.sim.covs],
                                        covMeans[bag,])

        }else{
          indx <- c() # vector of indices corresponding to rows being sub-sampled

          for (i in 1:nrow(strap.table)) {
            elec.indx <- which(data$Study == as.numeric(as.character(strap.table[i, 1]))) # rows corresponding to electrode sampled in current bag
            elec.obs <- round(length(elec.indx) / bag.size * as.numeric(as.character(strap.table[i, 2])))# number of rows to sample times the number of times the electrode shows up in the study strap generator. Divide rows of current electrode by the bag.size (14 electrodes)
            rows.samp <- sample(elec.indx, elec.obs, replace = FALSE)
            indx <- c(indx, rows.samp) # sample as many rows as indicated in elec.obs and add to indices vector

          }

          rm(elec.indx)
          rm(rows.samp)
          ########################
          # similarity measure
          ########################
          current.sim.metric <- sim.fn( target.study[,target.sim.covs],
                                        data[indx, sim.covs] )
        }

        ########################
        # accept/reject step
        ########################

        if(current.sim.metric > sim.metric){
          # test to see if the similarity metric is at least as good as before
          sim.metric <- current.sim.metric # set the current similarity metric to the current metric
          counter <- counter + 1
          total.straps <- total.straps + 1 # add to total straps
          sim <- 1 # break the while loop

        }else{
          # if it isn't keep testing new study straps
          counter <- counter + 1
          total.straps <- total.straps + 1 # add to total straps
        }

      }

      #-------------------------------------------------------------------------

      # only fit classifier if the while loop above was broken because of an acceptance
      # i.e., SKIP remaking study strap and fitting classifier if while loop was broken because it reached convergence limit
      if(counter < converge.lim){

        if( meanSampling == TRUE){
          # need to sample study strap if used meanSampling
          indx <- c() # vector of indices corresponding to rows being sub-sampled

          for (i in 1:nrow(strap.table)) {
            elec.indx <- which(data$Study == as.numeric(as.character(strap.table[i, 1]))) # rows corresponding to electrode sampled in current bag
            elec.obs <- round(length(elec.indx) / bag.size * as.numeric(as.character(strap.table[i, 2])))# number of rows to sample times the number of times the electrode shows up in the study strap generator. Divide rows of current electrode by the bag.size (14 electrodes)
            rows.samp <- sample(elec.indx, elec.obs, replace = FALSE)
            indx <- c(indx, rows.samp) # sample as many rows as indicated in elec.obs and add to indices vector

          }
          rm(elec.indx)
          rm(rows.samp)
        }

        row.list[[z]] <- indx # rows corresponding to current study strap
        convgVec <- c(convgVec, counter) # add number of candidate study straps before convergence

        #########################################
        # Full Classifier for Study Strap
        #########################################
        trCont <- trainControl(method = "none", trim = TRUE)

        for(mod in 1:num.SSLs){
          # iterate through SSLs

            if( anyNA(ssl.tuneGrid[[mod]]) ){

              ss.obj$models[[mod]][[z]] <- fatTrim(
                                            train( formula,
                                                   data = data[indx, -study.clmn],
                                                   method=ssl.method[[mod]],
                                                   trControl = trCont
                                                   ) )

            }else{
                
                #########################################
                # Tune
                #########################################
                if(!is.null(tune)){
                    if(tune == "cv"){
                        # ****** study balanced folds not working on cluster *******
                        # tune current strap with study balanced CV
                        
                        # trainList <- createFolds(factor(data$Study[indx]), 
                        #                          k = nfolds, 
                        #                          returnTrain = TRUE)
                        # 
                        # tOut <- lapply(trainList, function(x) indx[-x]) # testing on
                        # trainList <- lapply(trainList, function(x) indx[x]) # training on
                        trCont = caret::trainControl(method = "cv", number = nfolds) #,
                        
                        tuneModel <- caret::train(formula, 
                                                  data = data[indx, -study.clmn], # need to specify indx here because we did NOT in lapply
                                                  method = ssl.method[[mod]],
                                                  trControl= trCont, 
                                                  tuneGrid = tune.grid[[mod]], 
                                                  metric = "RMSE")
                        #                             index = trainList,
                        #                             indexOut = tOut) 
                        
                    }else if(tune =="ho"){
                        # train on study strap and validate on the rest
                        trainList <- indx
                        tOut <- seq(1, nrow(data))[-indx]
                        
                        trCont <- caret::trainControl(method = "cv",
                                                      index = trainList,
                                                      indexOut = tOut)
                        
                        tuneModel <- caret::train(formula, 
                                                  data = data[, -study.clmn], # DO NOT specify index here because we did it in lapply
                                                  method = ssl.method[[mod]],
                                                  trControl= trCont, 
                                                  tuneGrid = tune.grid[[mod]], 
                                                  metric = "RMSE")
                        
                    }
                    
                    
                    ssl.tuneGrid[[mod]] <- tuneModel$bestTune # replace with current best
                    rm(tuneModel)
                }
            
              trCont <- trainControl(method = "none", trim = TRUE)
              ss.obj$models[[mod]][[z]] <- fatTrim(
                                            train(formula,
                                                  data = data[indx, -study.clmn],
                                                  method=ssl.method[[mod]],
                                                  trControl = trCont,
                                                  tuneGrid = ssl.tuneGrid[[mod]] ) )
            }

        }

        # save memory and add model to object to be outputted

        if ( sim.mets == TRUE ){
          # if sim.metrics is provided by user, then calculate metrics to generate CPS weights
          similarity.matrix[z, ] <- sim.metrics(target.study[,target.sim.covs],
                                                data[indx, sim.covs])
        }

        if( length(customFNs) > 0){
          # add custom function CPS weights
          for(fn in 1:length(customFNs)){
            custom.matrix[z,fn] <- customFNs[[fn]](target.study[,target.sim.covs],
                                                   data[indx, sim.covs])
          }
        }

        z <- z + 1 # current accepted pseudo electrode counter
        rm(indx)
      }
    }

    # remove empty elements of model list
    for(mod in 1:num.SSLs){
      ss.obj$models[[mod]] <- ss.obj$models[[mod]][ lapply(ss.obj$models[[mod]], length) > 0 ] # remove elements of list that are empty
    }
    z <- length(ss.obj$models[[1]])

        ##############
        # Stacking
        ##############

        indx.count <- 1 # counts index

        if (stack == "ss") {
            row.list <- unlist(row.list) # vectorize rows of all study straps

            # stacking regression matrix
            stack.mat <- matrix(ncol = z * num.SSLs + 1, nrow = length(row.list))
            stack.mat[, 1] <- data[row.list, Yname] # data$Y[row.list] # add labels

            for(mod in 1:num.SSLs){

                for (SSL in 1:z ) {
                  indx.count <- indx.count + 1

                  stack.mat[, indx.count] <- predict(ss.obj$models[[mod]][[SSL]],
                                                  data[row.list,])

                }
            }

            if(stackTune){
                # tune stack parameter
                modTune <- cv.glmnet(y = as.vector(stack.mat[, 1]), 
                                     x = as.matrix(stack.mat[,-1]),
                                     alpha = 0, # no lasso penalty
                                     #foldid = data[,1],
                                     lambda = stack.lambda, # ridge penalty
                                     # standardize = FALSE,
                                     intercept = stackIntercept,
                                     lower.limits = wLim) 
                
                stackParam <- modTune$lambda.min # replace old value with tuned one
                ss.obj$stackLambda <- stackParam
                rm(modTune)
            }
            
            mod <- glmnet(y = as.vector(stack.mat[, 1]), 
                          x = as.matrix(stack.mat[,-1]),
                          alpha = 0,          # no lasso penalty
                          lambda = stackParam, # ridge penalty
                          # standardize = FALSE,
                          intercept = stackIntercept,
                          lower.limits = wLim) 
            
            ss.obj$stack.coefs <- as.vector( coef( mod ) )
            
            rm(mod)
            
            if(stack.standardize){
                ss.obj$stack.coefs <- ss.obj$stack.coefs / sum(ss.obj$stack.coefs)
            }
            
        } else if (stack == "standard") {
            
            indx.count <- 1 # counts index
            stack.mat <- matrix(ncol = z * num.SSLs + 1, nrow = nrow(data))
            stack.mat[, 1] <- as.numeric( data[,Yname] )# data$Y add labels
            
            for(mod in 1:num.SSLs){
                
                for (SSL in 1:z ) {
                    indx.count <- indx.count + 1
                    
                    stack.mat[, indx.count] <- predict(ss.obj$models[[mod]][[SSL]],
                                                       data )
                    
                }
            }
            
            if(stackTune){
                # tune stack parameter
                modTune <- cv.glmnet(y = as.vector(stack.mat[, 1]), 
                                     x = as.matrix(stack.mat[,-1]),
                                     alpha = 0, # no lasso penalty
                                     #foldid = data[,1],
                                     lambda = stack.lambda, # ridge penalty
                                     # standardize = FALSE,
                                     intercept = stackIntercept,
                                     lower.limits = wLim) 
                
                stackParam <- modTune$lambda.min # replace old value with tuned one
                ss.obj$stackLambda <- stackParam
                rm(modTune)
            }
            
            mod <- glmnet(y = as.vector(stack.mat[, 1]), 
                          x = as.matrix(stack.mat[,-1]),
                          alpha = 0,          # no lasso penalty
                          lambda = stackParam, # ridge penalty
                          # standardize = FALSE,
                          intercept = stackIntercept,
                          lower.limits = wLim) 
            
            ss.obj$stack.coefs <- as.vector( coef( mod ) )
            
            rm(mod)
            
            if(stack.standardize){
                ss.obj$stack.coefs <- ss.obj$stack.coefs / sum(ss.obj$stack.coefs)
            }
            
        }

        rm(stack.mat)
        if(model == TRUE){
          ss.obj$data <- data
        }

        ss.obj$strapRows <- row.list # rows of each observation in corresponding study strap
        ss.obj$modelInfo$convg.vec <- convgVec # vector of candidate straps before acceptance
        ss.obj$modelInfo$numStraps <- z



        # generate CPSW weights
        if( !anyNA(target.study) ){
          # if a target study is provided generate matrix of CPS weights
          similarity.matrix <- prop.table(abs(similarity.matrix[1:z,]), 2)
          #similarity.matrix <- abs(similarity.matrix[1:z,]) / colSums( abs(similarity.matrix[1:z,]) )

          # each column is a vector of weights that could be used to weight the study strap predictions
          colnames(similarity.matrix) <- c("Matcor Diag", "Matcor Sum", "Matcor Sum Abs", "|rho|",
                                           "rho sq", "UV rho sq", "UV cov sq", "UV rho", "UV cov",
                                           "diag UV rho sq", "diag UV cov", "diag UV cov sq", "Mean Corr",
                                           "SMI", "RV", "RV2", "RVadj", "PSI", "r1", "r2", "r3", "r4", "GCD")
          rownames(similarity.matrix) <- paste0("SS_", 1:z) # each row corresponds to a study strap

          # CPS weights
          if ( sim.mets == TRUE && length(customFNs) == 0){
            # default similarity metrics only
            ss.obj$simMat <- similarity.matrix # matrix of similarity measures: remove rows with NAs

          }else if ( sim.mets == TRUE && length(customFNs) > 0 ){
            # default similarity metrics and custom similarity metrics
            custom.matrix <- prop.table(abs(as.matrix(custom.matrix[1:z,])), 2)
            colnames(custom.matrix) <- paste("customFN", 1:length(customFNs))
            ss.obj$simMat <- cbind( similarity.matrix, custom.matrix ) # matrix of similarity measures: remove rows with NAs
            # colnames(ss.obj$simMat)[ncol(ss.obj$simMat)] <- "customFn"

          }else if (sim.mets == FALSE && length(customFNs) > 0 ){
            # no default similarity metrics, only custom similarity metrics
            custom.matrix <- prop.table(abs(as.matrix(custom.matrix[1:z,])), 2)
            colnames(custom.matrix) <- paste("customFN", 1:length(customFNs))
            ss.obj$simMat <- custom.matrix  # matrix of similarity measures: remove rows with NAs
            # colnames(ss.obj$simMat)[ncol(ss.obj$simMat)] <- "customFn"
            }

        }

    return(ss.obj)
}
