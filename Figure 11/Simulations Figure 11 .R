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
    
    
    
    ## box and whisker plot with color and no legend
    Method <- rep(c("TOS", "TOS-Stack",  "TOS-CPS", "TSS", "TSS-Stack", "TSS-CPS", "AR",  "AR-Stack","AR-CPS"), 
                  each = 100)

    PCR_gg2 <- as.vector(PCR_table)
    PCR_gg2 <- data.frame(RMSE = PCR_gg2, Methods = Method ) 
    colnames(PCR_gg2) <- c("RMSE", "Methods")
    PCR_gg2[,2] <- factor(PCR_gg2[,2])

    #####################################################
}

bagSizeComp <- cbind(sim.matrix,tuneComparison)

cbind(res.mat,testVsOpt)





colOrds <- c(1,3,2,4,6,5,7,9,8) # order columns to match figures
# No clusters table
indx <- 1:12


#########################################
# bag size and sigma^2_beta relationship just for clusters (c = 4)

dfBag <- as.data.frame( cbind(sim.matrix[13:24,-c(1,4)], tuneComparison[13:24,10] ) )
colnames(dfBag) <- c("betaVar", "XVar", "Opt")

p <- ggplot(data = dfBag, aes(x=XVar, y=Opt, colour = factor(betaVar))) +
    geom_line(size = 2) + 
    theme_classic(base_size = 12) + 
    #geom_jitter(height = 0.2, width = 0.75, size = 3) 
    ylab("Optimal Bag Size") + 
    # aes(fill=factor(betaVar))
    geom_point( size = 4 )+ xlab(TeX(sprintf('$\\sigma_{X}}'))) + 
    theme(
        plot.title = element_text(hjust = 0.5, color="red", size=14, face="bold")) +
    theme(axis.text.x = element_text(face="bold",  
                                     size=rel(2) ), 
          axis.title.x = element_text(face="bold", 
                                      size=rel(2) ),
          axis.title.y = element_text(face="bold", 
                                     size=rel(2),
                                     vjust = 2),
          axis.text.y = element_text(face="bold", 
                                     size=rel(2) ),
          #legend.key.size = unit(3, "line"),
          legend.text =  element_text(face="bold", 
                                        size=rel(2) ),
          legend.title = element_text(face="bold", 
                                      size=rel(2) ),
          legend.title.align = 0.7
          ) + 
    labs(color=TeX(sprintf('$\\sigma_{\\beta}^2}')),
         size = rel(2)) 


ggsave("Sims5_BagSize_vs_X.png", 
       plot =p,
       width = 7, height = 5)

