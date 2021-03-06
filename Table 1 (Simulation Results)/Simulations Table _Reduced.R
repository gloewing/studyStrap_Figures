library(ggplot2)
library(gridExtra)
library(viridis)
library(knitr)
library(kableExtra)
library(grid)
library(latex2exp)
library(tidyverse)
# labels and tables

Mode = function(x){
    ta = table(x)
    tam = max(ta)
    if (all(ta == tam))
        mod = NA
    else
        if(is.numeric(x))
            mod = as.numeric(names(ta)[ta == tam])
    else
        mod = names(ta)[ta == tam]
    return(min(mod) ) # choose minimum when there are ties
}

# labels and tables
get_legend<-function(myggplot){
    tmp <- ggplot_gtable(ggplot_build(myggplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
}

setwd("~/Desktop/Research/sims5")

sim.name <- c(160:191) # order to correspond with figure # add 138 later
num.sims <- length(sim.name)
#cluster.list <- rep(c(0,8,4), each = 4)
#sigma.list <- rep(c(0.05,0.25,1,3), each = 3)


sims <- cbind(160:183, expand.grid(c(0.05, 0.25, 1, 3), c(0.05, 5, 20), c(0,4)))
sims2 <- cbind( 184:191, expand.grid(c(0.05,0.25,1,3),c(10,15), c(0)))
colnames(sims) <- c("simNum", "betaVar", "XVar", "Clusts")
colnames(sims2) <- c("simNum", "betaVar", "XVar", "Clusts")
sims <- rbind(sims,sims2)
sim.matrix <- sims

bagSize.vec <- c(1:3, seq(4,18, by = 2), seq(20,100, by = 10), 200, 500, 1000)
labs <- c("TOS", "TOS-Stack", "TOS-CPS", "TSS",  "TSS-Stack", "TSS-CPS", "AR", "AR-CPS", "AR-Stack")

res.mat <- matrix(ncol = 9, nrow = length(160:191))
res.matMode <- matrix(ncol = 9, nrow = length(160:191))
colnames(res.mat) <- c("TOS", "TOS-Stack",  "TOS-CPS", "TOS", "TOS-Stack", "TOS-CPS", "AR",  "AR-Stack","AR-CPS")
colnames(res.matMode) <- c("TOS", "TOS-Stack",  "TOS-CPS", "TSS", "TSS-Stack", "TSS-CPS", "AR",  "AR-Stack","AR-CPS")


ARtestmat <- matrix(ncol = length(bagSize.vec), nrow = length(160:191))
SStestmat <- matrix(ncol = length(bagSize.vec), nrow = length(160:191))

# use the mode to tune


# comparison between optimal bag size for tune and test
tuneComparison <- matrix(ncol = 12, nrow = length(160:191) )
colnames(tuneComparison) <- c("SSTune", "SSTune_stack","SSTune_cps","SSTest", "SSTest_stack","SSTest_cps", 
                              "arTune", "arTune_stack","arTune_cps", "arTest","arTest_stack","arTest_cps")

# comparison between optimal bag size for tune and test Mode
tuneComparisonMode <- matrix(ncol = 12, nrow = length(160:191) )
colnames(tuneComparisonMode) <- c("SSTune", "SSTune_stack","SSTune_cps","SSTest", "SSTest_stack","SSTest_cps", 
                              "arTune", "arTune_stack","arTune_cps", "arTest","arTest_stack","arTest_cps")

# comparison between the RMSE of the actual version vs. the tested version
testVsOpt <- matrix(ncol = 4, nrow  = length(160:191))
colnames(testVsOpt) <- c("SSTuned", "SSoptTest", "ARTuned", "ARoptTes")

# save standard error estimates of monte carlo simulations
seMat <- matrix(nrow = num.sims, ncol = 9)
colnames(seMat) <- c("TOS", "TOS-Stack",  "TOS-CPS", "TSS", "TSS-Stack", "TSS-CPS", "AR",  "AR-Stack","AR-CPS")
###########################
# Individual simulations
###########################
for( z in 1:num.sims){
    # perturbation = 0.1 * sigma_beta
    pert <- "0.1_sig"
    #####################################################
    setwd("~/Desktop/Research/sims5")
    
    num <- sim.matrix[z,1]
    sigma <- sim.matrix[z,2] # sigma beta
    clusts <- sim.matrix[z,4]
    sig.mu <- sim.matrix[z,3] # sigma X
    
    ###### optimal bag size ########

    dat <- read.csv(paste0("LASSO_Sims5_Results_",num))
    merged <- dat[,1]
    
    # scale by mean of merged
    PCR_table <- dat[,-1] / merged # standardize by merged
    colnames(PCR_table) <- c("TOS", "TOS-Stack",  "TOS-CPS", "TSS", "TSS-Stack", "TSS-CPS", "AR",  "AR-Stack","AR-CPS")
    merge.mean <- round(mean(merged),2)

    modeMat <- PCR_table
    
    res.mat[z,] <- colMeans(PCR_table) # average results
    testVsOpt[z,1] <-  mean(dat[,5]/ merged)# comparison to oracle of bag size SS
    testVsOpt[z,3] <-  mean(dat[,8]/ merged)# comparison to oracle of bag size AR
    
    ###### optimal tune SS bag size
    # avg
    SS <- read.csv(paste0("LASSO SS Tune_Sims5_Sim ",num, "_Combined"))
    meanSSbagSizeTest <- colMeans(SS)
    opt.bagSz.SS <- bagSize.vec[which.min(meanSSbagSizeTest)]
    tuneComparison[z,1] <- opt.bagSz.SS
    minApply <- apply(SS,1, which.min)
    tuneComparisonMode[z,1] <- bagSize.vec[Mode(minApply)]
    SSTuneMode <- Mode(minApply)
    
    #SS stack
    SS <- read.csv(paste0("LASSO SS Tune_Sims5_Sim ",num, "_Combined_Stack"))
    meanSSbagSizeTest <- colMeans(SS)
    tuneComparison[z,2] <- bagSize.vec[which.min(meanSSbagSizeTest)]
    minApply <- apply(SS,1, which.min)
    tuneComparisonMode[z,2] <- bagSize.vec[Mode(minApply)]
    SSStackTuneMode <- Mode(minApply)
    
    #SS cps
    SS <- read.csv(paste0("LASSO SS Tune_Sims5_Sim ",num, "_Combined_CPS"))
    meanSSbagSizeTest <- colMeans(SS)
    tuneComparison[z,3] <- bagSize.vec[which.min(meanSSbagSizeTest)]
    minApply <- apply(SS,1, which.min)
    tuneComparisonMode[z,3] <- bagSize.vec[Mode(minApply)]
    SSCPSTuneMode <- Mode(minApply)
    rm(SS)
    #----------------
    
    #### optimal tune AR bag size
    # avg
    ar <- read.csv(paste0("LASSO AR Tune_Sims5_Sim ",num, "_Combined"))
    meanARbagSizeTest <- colMeans(ar)
    opt.bagSz <- bagSize.vec[which.min(meanARbagSizeTest)]
    tuneComparison[z,7] <- opt.bagSz
    minApply <- apply(ar,1, which.min)
    tuneComparisonMode[z,7] <- bagSize.vec[Mode(minApply)]
    minApply <- apply(ar,1, which.min)
    tuneComparisonMode[z,7] <- bagSize.vec[Mode(minApply)]
    arTuneMode <- Mode(minApply)
    
    #ar stack
    ar <- read.csv(paste0("LASSO AR Tune_Sims5_Sim ",num, "_Combined_Stack"))
    meanARbagSizeTest <- colMeans(ar)
    tuneComparison[z,8] <- bagSize.vec[which.min(meanARbagSizeTest)]
    minApply <- apply(ar,1, which.min)
    tuneComparisonMode[z,8] <- bagSize.vec[Mode(minApply)]
    arStackTuneMode <- Mode(minApply)
    
    #ar cps
    ar <- read.csv(paste0("LASSO AR Tune_Sims5_Sim ",num, "_Combined_CPS"))
    meanARbagSizeTest <- colMeans(ar)
    tuneComparison[z,9] <- bagSize.vec[which.min(meanARbagSizeTest)]
    minApply <- apply(ar,1, which.min)
    tuneComparisonMode[z,9] <- bagSize.vec[Mode(minApply)]
    arCPSTuneMode <- Mode(minApply)
    rm(ar)
    #----------------
    
    # optimal test SS/AR bag size
    # avg SS
    SStestmat[z,] <- colMeans(read.csv(paste0("LASSO SS BagSize_Test_Sims5_",num,"_Avg")))
    tuneComparison[z,4] <- bagSize.vec[which.min(SStestmat[z,])]
    testVsOpt[z,4] <- colMeans( read.csv(paste0("LASSO SS BagSize_Test_Sims5_",num,"_Avg"))[which.min(SStestmat[z,])] / merged )
    SS <- read.csv(paste0("LASSO SS BagSize_Test_Sims5_",num,"_Avg"))
    minApply <- apply(SS,1, which.min)
    tuneComparisonMode[z,4] <- bagSize.vec[Mode(minApply)]
    modeMat[,4] <- read.csv(paste0("LASSO SS BagSize_Test_Sims5_",num,"_Avg"))[,SSTuneMode]/merged
    
    # stack SS
    SS <- colMeans(read.csv(paste0("LASSO SS BagSize_Test_Sims5_",num,"_Stack")))
    tuneComparison[z,5] <- bagSize.vec[which.min(SS)]
    rm(SS)
    SS <- read.csv(paste0("LASSO SS BagSize_Test_Sims5_",num,"_Stack"))
    minApply <- apply(SS,1, which.min)
    tuneComparisonMode[z,5] <- bagSize.vec[Mode(minApply)]
    modeMat[,5] <- read.csv(paste0("LASSO SS BagSize_Test_Sims5_",num,"_Stack"))[,SSStackTuneMode]/merged
    rm(SS)
    
    # cps SS
    SS <- colMeans(read.csv(paste0("LASSO SS BagSize_Test_Sims5_",num,"_CPS")))
    tuneComparison[z,6] <- bagSize.vec[which.min(SS)]
    rm(SS)
    SS <- read.csv(paste0("LASSO SS BagSize_Test_Sims5_",num,"_CPS"))
    minApply <- apply(SS,1, which.min)
    tuneComparisonMode[z,6] <- bagSize.vec[Mode(minApply)]
    modeMat[,6] <- read.csv(paste0("LASSO SS BagSize_Test_Sims5_",num,"_CPS"))[,SSCPSTuneMode]/merged
    rm(SS)
    
    # avg ar
    ARtestmat[z,] <- colMeans(read.csv(paste0("LASSO AR BagSize_Test_Sims5_",num,"_Avg")))
    tuneComparison[z,10] <- bagSize.vec[which.min(ARtestmat[z,])]
    testVsOpt[z,4] <- colMeans( read.csv(paste0("LASSO AR BagSize_Test_Sims5_",num,"_Avg"))[which.min(SStestmat[z,])] / merged )
    ar <- read.csv(paste0("LASSO AR BagSize_Test_Sims5_",num,"_Avg"))
    minApply <- apply(ar,1, which.min)
    tuneComparisonMode[z,10] <- bagSize.vec[Mode(minApply)]
    modeMat[,7] <- read.csv(paste0("LASSO AR BagSize_Test_Sims5_",num,"_Avg"))[,arTuneMode]/merged
    rm(ar)
    
    # stack ar
    ar <- colMeans(read.csv(paste0("LASSO AR BagSize_Test_Sims5_",num,"_Stack")))
    tuneComparison[z,11] <- bagSize.vec[which.min(ar)]
    rm(ar)
    ar <- read.csv(paste0("LASSO AR BagSize_Test_Sims5_",num,"_Stack"))
    minApply <- apply(ar,1, which.min)
    tuneComparisonMode[z,11] <- bagSize.vec[Mode(minApply)]
    modeMat[,8] <- read.csv(paste0("LASSO AR BagSize_Test_Sims5_",num,"_Stack"))[,arStackTuneMode]/merged
    rm(ar)
    
    # cps ar
    ar <- colMeans(read.csv(paste0("LASSO AR BagSize_Test_Sims5_",num,"_CPS")))
    tuneComparison[z,12] <- bagSize.vec[which.min(ar)]
    rm(ar)
    ar <- read.csv(paste0("LASSO AR BagSize_Test_Sims5_",num,"_CPS"))
    minApply <- apply(ar,1, which.min)
    tuneComparisonMode[z,12] <- bagSize.vec[Mode(minApply)]
    modeMat[,9] <- read.csv(paste0("LASSO AR BagSize_Test_Sims5_",num,"_CPS"))[,arCPSTuneMode] /merged
    rm(ar)
    
    res.matMode[z,] <- colMeans(modeMat) # average results using Mode as tuning
    #-------------------------
    
    PCR_tableAug <- cbind(PCR_table)
    
    # standard error of monte carlo
    seMat[z,] <- apply(PCR_table[,1:9], 2, sd) / sqrt( nrow(PCR_table)  ) # \hat{s} / sqrt{n} where n is number of monte carlo sims
    
    ## box and whisker plot with color and no legend
    Method <- rep(c("TOS", "TOS-Stack",  "TOS-CPS", "TSS", "TSS-Stack", "TSS-CPS", "AR",  "AR-Stack","AR-CPS"), 
                  each = 100)

    PCR_gg2 <- as.vector(PCR_table)
    PCR_gg2 <- data.frame(RMSE = PCR_gg2, Methods = Method ) 
    colnames(PCR_gg2) <- c("RMSE", "Methods")
    PCR_gg2[,2] <- factor(PCR_gg2[,2])

    #####################################################
}

# add on parameter values to monte carlo error matrix (seMat)

colOrds <- c(1,3,2,4,6,5,7,9,8) # order columns to match figures

param <- rep(sim.matrix[,2], each = 2)
s <- seq(1, length(param), by = 2)

reducedMat <- matrix(NA, nrow = 2 * nrow(sim.matrix),
                     ncol = length(colOrds) + 1)

# parameters
reducedMat[s, 1] <- paste0( round(sim.matrix[,2], 3), " (",
                            (sim.matrix[,3])^2, ")")
# RMSE
reducedMat[s, -1] <- round( res.mat, digits = 2 )
# standard error
reducedMat[s + 1, -1] <- paste0( "(", signif( seMat, digits = 2 ), ")" )

#######################################################
#####################
# reduced table
#####################
# No clusters table
indx <- c(1, 2, 7, 8, 17, 18, 23, 24) #1:24
kable( reducedMat[indx, c(1, colOrds + 1) ], 
      format = "latex", booktabs = T) %>% kable_styling(position = "center")

# res.mat[1:12, colOrds]


# clusters table
indx <- c( 25, 26, 31, 32, 41, 42, 47, 48 )
kable( reducedMat[indx, c(1, colOrds + 1) ], 
       format = "latex", booktabs = T) %>% kable_styling(position = "center")

#########################
# full table (supplement)
#########################
# No clusters table
indx <- 1:24
kable( reducedMat[indx, c(1, colOrds + 1) ], 
       format = "latex", booktabs = T) %>% kable_styling(position = "center")

# clusters table
indx <- 25:48
kable( reducedMat[indx, c(1, colOrds + 1) ], 
       format = "latex", booktabs = T) %>% kable_styling(position = "center")
