#' The Study Strap for Multi-Study Learning: Fits Study Strap algorithm
#'
#' @param formula Model formula
#' @param data A dataframe with all the studies has the following columns in this order: "Study", "Y", "V1", ...., "Vp"
#' @param target.study Dataframe of the design matrix (just covariates) of study one aims to make predictions on
#' @param bag.size Integer indicating the bag size tuning parameter.
#' @param straps Integer indicating the maximum number of study straps to generate and fit models with.
#' @param stack String taking values "standard" or "ss" specifying how to fit the stacking regression. "standard" option uses all studies as the "test" studies. "ss" uses all the study straps as "test" studies.
#' @param sim.covs Is a vector of names of covariates or the column numbers of the covariates to be used for the similarity measure. Default is to use all covariates.
#' @param ssl.method A list of strings indicating which modeling methods to use.
#' @param ssl.tuneGrid A list of the tuning parameters in the format of the caret package. Each element must be a dataframe (as required by caret). If no tuning parameters are required then NA is indicated.
#' @param sim.mets Boolean indicating whether to calculate default covariate profile similarity measures.
#' @param model Indicates whether to attach training data to model object.
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
#' ssMod1 <- ss(formula = Y ~.,
#'             data = data,
#'             target.study = target,
#'             bag.size = length(unique(data$Study)),
#'             straps = 5,
#'             stack = "standard",
#'             sim.covs = NA,
#'             ssl.method = list("pcr"),
#'             ssl.tuneGrid = list(data.frame("ncomp" = 1)),
#'             sim.mets = TRUE,
#'             model = TRUE,
#'             customFNs = list() )
#'
#'# Fit model with 2 SSLs: Linear Regression and PCA Regression
#' \donttest{ssMod2 <- ss(formula = Y ~.,
#'             data = data,
#'             target.study = target,
#'             bag.size = length(unique(data$Study)),
#'             straps = 10,
#'             stack = "standard",
#'             sim.covs = NA,
#'             ssl.method = list("lm","pcr"),
#'             ssl.tuneGrid = list(NA, data.frame("ncomp" = 2)),
#'             sim.mets = TRUE,
#'             model = TRUE,
#'             customFNs = list( ) )}
#'
#'
#'
#' # Fit model with custom similarity function for
#' # covariate profile similarity weighting
#'
#' fn1 <- function(x1,x2){
#' return( abs( cor( colMeans(x1), colMeans(x2) )) )
#' }
#'
#' \donttest{ssMod3<- ss(formula = Y ~.,
#'             data = data,
#'             target.study = target,
#'             bag.size = length(unique(data$Study)),
#'             straps = 10,
#'             stack = "standard",
#'             sim.covs = NA,
#'             ssl.method = list("lm","pcr"),
#'             ssl.tuneGrid = list(NA, data.frame("ncomp" = 2)),
#'             sim.mets = TRUE,
#'             model = TRUE, customFNs = list(fn1) )}
#'
#' #########################
#' #####  Predictions ######
#' #########################
#'
#' preds <- studyStrap.predict(ssMod1, target)
#' @import caret
#' @import tidyverse
#' @import nnls
#' @import dplyr
#' @importFrom stats coef cor cov predict
#' @export

ss <- function(formula = Y ~.,
               data,
               target.study = NA,
               bag.size = length(unique(data$Study)),
               straps = 150,
               mtry = NA,
               stack = "standard",
               stackTune = TRUE,
               stackNNLS = TRUE,
               stackIntercept = TRUE,
               sim.covs = NA,
               ssl.method = list("lm"),
               ssl.tuneGrid = list(c()),
               sim.mets = FALSE,
               model = FALSE,
               customFNs = list(),
               stack.standardize = FALSE,
               featureMatrix = NA,
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
  # model determines whether data is attached
  # customFNs is a list where each element is a custom function for CPS weighting (not used for accept/reject step)
  # sim.mets is a indicator of whether to use the default similarity measures. Requires this and target study design matrix to calculate
  # stack.standardize indicates whether stacking coefs are standardized
  # mtry is the number of features to be randomly included
  # featureMatrix is a matrix of features that are used for feature bagging -- number of rows equal to number of study straps
  # stackTune is an indicator of whether to do cross validation for stacking
    
    # library(dplyr)
    # library(caret)
    # library(nnls)
    
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

    if( anyNA(mtry) ){
        # if number of features not specified make it the number of features included in dataset
        mtry <- ncol(data) - 2 # exclude Study and Y when counting number of features
    }

    target.sim.covs <- sim.covs - 2 # the target study does not have Study and outcome (Y) columns
    # so shift indices over by 2

    Study.Code <- NULL # include for namespace

    # rename studies as integers in order from 1:length(unique(data$Study))
    study.mat <- tibble::tibble(Study = unique(data$Study), Study.Code = 1:length(unique(data$Study)))
    data <- data %>%
      dplyr::left_join(study.mat, by = "Study") %>%
      dplyr::mutate(Study = Study.Code) %>%
      dplyr::select(-Study.Code)
    rm(study.mat)

    studies <- unique(data$Study)
    row.list <- list(length = straps) # list of the rows that were included in study strap

    # determine number of rows in each study
    sampSizes <- as.vector( table(data$Study) )

    model.list <- vector("list", length = num.SSLs)

    for(mod in 1:num.SSLs){
      # each element in list is a list of models for each study strap
      model.list[[mod]] <- vector("list", length = straps)
    }

    ss.obj <- list(models = model.list, data = c(), sim.mat <- c(),
                   strapRows <- list(),
                   dataInfo = list(studyNames = original.studies, sampleSizes = sampSizes),
                   modelInfo = list(sampling = "ss", numStraps = straps, SSL = ssl.method,
                                    ssl.tuneGrid = c(), numPaths = NA,
                         convg.vec = c(), convgCritera = NA, meanSamp = NA, stack.type = stack,
                         custFNs = customFNs, bagSize = bag.size),
                   stack.coefs <- c(), simMat = c() )
    class(ss.obj) <- "ss"

    # similarity matrix contains the average covariates that are specified
    similarity.matrix <- matrix(ncol = 23, nrow = straps)  # similarity metrics

    if(length(customFNs) > 0){
      # if there are any custom functions for CPS weighting
      custom.matrix <- matrix(ncol = length(customFNs), nrow = straps) # custom similarity metrics
      colnames(custom.matrix) <- paste0("customFn_", 1:length(customFNs))
    }

    for (z in 1:straps) {

        message(paste0("Study Strap ", z))

        ########################
        ######## Study Bag
        ########################
        bag <- sample(studies, bag.size, replace = TRUE)
        strap.table <- as.data.frame(table(bag)) #sample with replacement and put in data table format
        
        #####################################################################
        # Study Strap Rows
        #####################################################################

        indx <- c() # vector of indices corresponding to rows being sub-sampled

        for (i in 1:nrow(strap.table)) {
            elec.indx <- which(data$Study == as.numeric(as.character(strap.table[i, 1]))) # rows corresponding to electrode sampled in current bag
            elec.obs <- round(length(elec.indx) / bag.size * as.numeric(as.character(strap.table[i, 2])))# number of rows to sample times the number of times the electrode shows up in the study strap generator. Divide rows of current electrode by the bag.size (14 electrodes)
            rows.samp <- sample(elec.indx, elec.obs, replace = FALSE)
            indx <- c(indx, rows.samp) # sample as many rows as indicated in elec.obs and add to indices vector

        }

        row.list[[z]] <- indx # rows corresponding to current study strap
        
        ########################
        ######## Feature Bag
        ########################
        if(all(is.na(featureMatrix))){
            # draw randomly if no feature matrix is given
            Featbag <- sample(3:ncol(data), mtry, replace = FALSE) # first two columns are Y and study
            Featbag <- sort( c(2, Featbag) )# add in index for Y for below
            
        }else{
            # if featureMatrix is given, use the zth row as the feature bag for the zth strap
            # each row of the featureMatrix is a matrix of covariate indices to be used
            Featbag <- as.vector( featureMatrix[z,] + 2 ) # add two because the first two columns are Study and Y
            Featbag <- sort( c(2, Featbag) )# add in index for Y for below
        }

        #########################################
        # Full Classifier for Study Strap
        #########################################

        trCont <- trainControl(method = "none", trim = TRUE)

        for(mod in 1:num.SSLs){
          # iterate through SSLs
          if( anyNA(ssl.tuneGrid[[mod]]) ){

            ss.obj$models[[mod]][[z]] <- fatTrim( caret::train(formula, data = data[indx, Featbag], # -study.clmn],
                                                        method=ssl.method[[mod]], trControl = trCont) )

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
            ss.obj$models[[mod]][[z]] <- fatTrim( caret::train(formula, data = data[indx, Featbag], # -study.clmn],
                                                        method=ssl.method[[mod]], trControl = trCont,
                                                        tuneGrid = ssl.tuneGrid[[mod]] ) )
          }

        }

        if ( sim.mets == TRUE ){
          # if a target study is provided by user, and generate CPS weights
          similarity.matrix[z, ] <- sim.metrics(target.study[,target.sim.covs], data[indx, sim.covs])
        }

        if( length(customFNs) > 0 ){
          # add custom function CPS weights
          for(fn in 1:length(customFNs)){
            custom.matrix[z,fn] <- customFNs[[fn]](target.study[,target.sim.covs],
                                                   data[indx, sim.covs])
          }


    }
        rm(indx)
    }
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

            for (SSL in 1:straps ) {
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

            for (SSL in 1:straps ) {
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
        
        if(model)   ss.obj$data <- data
        ss.obj$sim.mat <- sim.mat # matrix of average covariates
        ss.obj$strapRows <- row.list

        # generate CPSW weights
        if( !anyNA(target.study) ){
          # if a target study is provided generate matrix of CPS weights
          similarity.matrix <- prop.table(abs(similarity.matrix), 2)
          # each column is a vector of weights that could be used to weight the study strap predictions
          colnames(similarity.matrix) <- c("Matcor Diag", "Matcor Sum", "Matcor Sum Abs", "|rho|",
                                           "rho sq", "UV rho sq", "UV cov sq", "UV rho", "UV cov",
                                           "diag UV rho sq", "diag UV cov", "diag UV cov sq", "Mean Corr",
                                           "SMI", "RV", "RV2", "RVadj", "PSI", "r1", "r2", "r3", "r4", "GCD")
          rownames(similarity.matrix) <- paste0("SS_", 1:straps) # each row corresponds to a study strap

          # CPS weights
          if ( sim.mets == TRUE && length(customFNs) == 0){
            # default similarity metrics only
            ss.obj$simMat <- similarity.matrix # matrix of similarity measures: remove rows with NAs

          }else if ( sim.mets == TRUE && length(customFNs) > 0 ){
            # default similarity metrics and custom similarity metrics
            custom.matrix <- prop.table(abs(as.matrix(custom.matrix)), 2)
            ss.obj$simMat <- cbind( similarity.matrix, custom.matrix ) # matrix of similarity measures: remove rows with NAs

          }else if (sim.mets == FALSE && length(customFNs) > 0 ){
            # no default similarity metrics, only custom similarity metrics
            custom.matrix <- prop.table(abs(as.matrix(custom.matrix)), 2)
            ss.obj$simMat <- custom.matrix  # matrix of similarity measures: remove rows with NAs
          }

        }

    return(ss.obj)
}

