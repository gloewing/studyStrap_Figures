#' Merged Approach for Multi-Study Learning: fits a single model on all studies merged into a single dataframe.
#'
#' @param formula Model formula
#' @param data A dataframe with all the studies has the following columns in this order: "Study", "Y", "V1", ...., "Vp"
#' @param sim.covs Is a vector of names of covariates or the column numbers of the covariates to be used for the similarity measure. Default is to use all covariates.
#' @param ssl.method A list of strings indicating which modeling methods to use
#' @param ssl.tuneGrid A list of the tuning parameters in the format of the caret package. Each element must be a dataframe (as required by caret). If no tuning parameters are required then NA is indicated
#' @param model Indicates whether to attach training data to model object
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
#' mrgMod1 <- merged(formula = Y ~.,
#'                   data = data,
#'                  sim.covs = NA,
#'                  ssl.method = list("pcr"),
#'                  ssl.tuneGrid = list( data.frame("ncomp" = 2)),
#'                  model = FALSE )
#'
#' # 2 SSLs: Linear Regression and PCA Regression
#' mrgMod2 <- merged(formula = Y ~.,
#'                   data = data,
#'                  sim.covs = NA,
#'                  ssl.method = list("lm", "pcr"),
#'                  ssl.tuneGrid = list(NA,
#'                            data.frame("ncomp" = 2) ),
#'                  model = FALSE )
#'
#' #########################
#' #####  Predictions ######
#' #########################
#'
#' preds <- studyStrap.predict(mrgMod2, target)
#' @import caret
#' @import tidyverse
#' @import nnls
#' @import dplyr
#' @export

merged <- function(formula = Y ~.,
                   data,
                   sim.covs = NA,
                   ssl.method = list("lm"),
                   ssl.tuneGrid = list(c()),
                   model = FALSE ){
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
    #library(dplyr)
    #library(caret)

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

    Study.Code <- NULL # include for namespace

    # rename studies as integers in order from 1:length(unique(data$Study))
    study.mat <- tibble(Study = unique(data$Study), Study.Code = 1:length(unique(data$Study)))
    data <- data %>%
      left_join(study.mat, by = "Study") %>%
      mutate(Study = Study.Code) %>%
    dplyr::select(-Study.Code)
    rm(study.mat)

    studies <- unique(data$Study)
    row.list <- list(length = 1) # list of the rows that were included in study strap

    # determine number of rows in each study
    sampSizes <- c()
    for (i in studies){
      sampSizes <- c(sampSizes, length(data$Study == i) )
    }

    model.list <- vector("list", length = num.SSLs)

    for(mod in 1:num.SSLs){
      # each element in list is a list of models for each study strap
      model.list[[mod]] <- vector("list", length = 1)
    }

    ss.obj <- list(models = model.list, data = c(), sim.mat <- c(),
                   strapRows <- list(),
                   dataInfo = list(studyNames = original.studies, sampleSizes = sampSizes),
                   modelInfo = list(sampling = "merged", numStraps = 1, SSL = ssl.method,
                                    ssl.tuneGrid = c(), numPaths = NA,
                         convg.vec = c(), convgCritera = NA, meanSamp = NA, stack.type = NA,
                         custFNs = NA, bagSize = NA),
                   stack.coefs <- c(), simMat = c() )
    class(ss.obj) <- "ss"



    z <- 1

    #####################################################################
    # Merge Rows
    #####################################################################
    # all rows
    indx <- 1:nrow(data) # vector of indices corresponding to rows being sub-sampled

    row.list[[z]] <- indx # rows corresponding to current study strap

    #########################################
    # Full Classifier for Study Strap
    #########################################

    trCont <- trainControl(method = "none", trim = TRUE)

    for(mod in 1:num.SSLs){
      # iterate through SSLs

      message(paste0("SSL ", mod))
      if( anyNA(ssl.tuneGrid[[mod]]) ){

        ss.obj$models[[mod]][[z]] <- fatTrim( train(formula, data = data[indx, -study.clmn],
                                                    method=ssl.method[[mod]], trControl = trCont
                                                    ) )

      }else{
        ss.obj$models[[mod]][[z]] <- fatTrim( train(formula, data = data[indx, -study.clmn],
                                                    method=ssl.method[[mod]], trControl = trCont,
                                                    tuneGrid = ssl.tuneGrid[[mod]] ) )
      }

    }

      rm(indx)

      # fake stacking coefficients (simple average weights) for ss class format
      ss.obj$stack.coefs <- rep( NA,
                                 num.SSLs + 1)  # provide for prediction function

      if(model){
        ss.obj$data <- data
      }

      ss.obj$strapRows <- row.list


  return(ss.obj)
}

