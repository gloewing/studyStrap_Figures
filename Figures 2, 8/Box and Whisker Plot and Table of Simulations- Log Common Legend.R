library(ggplot2)
library(gridExtra)
library(viridis)
library(knitr)
library(kableExtra)
library(grid)
library(latex2exp)
library(tidyverse)
dslabs::ds_theme_set(new = "theme_minimal")

# labels and tables
get_legend<-function(myggplot){
    tmp <- ggplot_gtable(ggplot_build(myggplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
}

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

#Our transformation function to ensure rounding is same decimal places for all
scaleFUN <- function(x) sprintf("%.1f", x) # sprintf("%.2f", x)


setwd("~/Desktop/Research/sims5")

sim.name <- c(160:191) # order to correspond with figure # add 138 later
num.sims <- length(sim.name)
#cluster.list <- rep(c(0,8,4), each = 4)
#sigma.list <- rep(c(0.05,0.25,1,3), each = 3)


sims <- cbind(160:183, expand.grid(c(0.05, 0.25, 1, 3), c(0.05, 5, 20)^2, c(0,4)))
sims2 <- cbind( 184:191, expand.grid(c(0.05,0.25,1,3),c(10,15)^2, c(0)))
colnames(sims) <- c("simNum", "betaVar", "XVar", "Clusts")
colnames(sims2) <- c("simNum", "betaVar", "XVar", "Clusts")
sims <- rbind(sims,sims2)
sim.matrix <- sims

bagSize.vec <- c(1:3, seq(4,18, by = 2), seq(20,100, by = 10), 200, 500, 1000)
labs <- c("OSE", "OSE-Stack",  "OSE-CPS", "SSE", "SSE-Stack",
          "SSE-CPS", "AR",  "AR-Stack","AR-CPS")

res.mat <- matrix(ncol = 9, nrow = length(160:191))
res.matMode <- matrix(ncol = 9, nrow = length(160:191))
colnames(res.mat) <- c("OSE", "OSE-Stack",  "OSE-CPS", "SSE", "SSE-Stack",
                       "SSE-CPS", "AR",  "AR-Stack","AR-CPS")
colnames(res.matMode) <- c("OSE", "OSE-Stack",  "OSE-CPS", "SSE", "SSE-Stack",
                           "SSE-CPS", "AR",  "AR-Stack","AR-CPS")


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
    colnames(PCR_table) <- c("OSE", "OSE-Stack",  "OSE-CPS", "SSE", "SSE-Stack",
                             "SSE-CPS", "AR",  "AR-Stack","AR-CPS")
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
    modeMat[,4] <- read.csv(paste0("LASSO SS BagSize_Test_Sims5_",num,"_Avg"))[,SSTuneMode] / merged
    
    # stack SS
    SS <- colMeans(read.csv(paste0("LASSO SS BagSize_Test_Sims5_",num,"_Stack")))
    tuneComparison[z,5] <- bagSize.vec[which.min(SS)]
    rm(SS)
    SS <- read.csv(paste0("LASSO SS BagSize_Test_Sims5_",num,"_Stack"))
    minApply <- apply(SS,1, which.min)
    tuneComparisonMode[z,5] <- bagSize.vec[Mode(minApply)]
    modeMat[,5] <- read.csv(paste0("LASSO SS BagSize_Test_Sims5_",num,"_Stack"))[,SSStackTuneMode] / merged
    rm(SS)
    
    # cps SS
    SS <- colMeans(read.csv(paste0("LASSO SS BagSize_Test_Sims5_",num,"_CPS")))
    tuneComparison[z,6] <- bagSize.vec[which.min(SS)]
    rm(SS)
    SS <- read.csv(paste0("LASSO SS BagSize_Test_Sims5_",num,"_CPS"))
    minApply <- apply(SS,1, which.min)
    tuneComparisonMode[z,6] <- bagSize.vec[Mode(minApply)]
    modeMat[,6] <- read.csv(paste0("LASSO SS BagSize_Test_Sims5_",num,"_CPS"))[,SSCPSTuneMode] / merged
    rm(SS)
    
    # avg ar
    ARtestmat[z,] <- colMeans(read.csv(paste0("LASSO AR BagSize_Test_Sims5_",num,"_Avg")))
    tuneComparison[z,10] <- bagSize.vec[which.min(ARtestmat[z,])]
    testVsOpt[z,4] <- colMeans( read.csv(paste0("LASSO AR BagSize_Test_Sims5_",num,"_Avg"))[which.min(SStestmat[z,])] / merged )
    ar <- read.csv(paste0("LASSO AR BagSize_Test_Sims5_",num,"_Avg"))
    minApply <- apply(ar,1, which.min)
    tuneComparisonMode[z,10] <- bagSize.vec[Mode(minApply)]
    modeMat[,7] <- read.csv(paste0("LASSO AR BagSize_Test_Sims5_",num,"_Avg"))[,arTuneMode] / merged
    rm(ar)
    
    # stack ar
    ar <- colMeans(read.csv(paste0("LASSO AR BagSize_Test_Sims5_",num,"_Stack")))
    tuneComparison[z,11] <- bagSize.vec[which.min(ar)]
    rm(ar)
    ar <- read.csv(paste0("LASSO AR BagSize_Test_Sims5_",num,"_Stack"))
    minApply <- apply(ar,1, which.min)
    tuneComparisonMode[z,11] <- bagSize.vec[Mode(minApply)]
    modeMat[,8] <- read.csv(paste0("LASSO AR BagSize_Test_Sims5_",num,"_Stack"))[,arStackTuneMode] / merged
    rm(ar)
    
    # cps ar
    ar <- colMeans(read.csv(paste0("LASSO AR BagSize_Test_Sims5_",num,"_CPS")))
    tuneComparison[z,12] <- bagSize.vec[which.min(ar)]
    rm(ar)
    ar <- read.csv(paste0("LASSO AR BagSize_Test_Sims5_",num,"_CPS"))
    minApply <- apply(ar,1, which.min)
    tuneComparisonMode[z,12] <- bagSize.vec[Mode(minApply)]
    modeMat[,9] <- read.csv(paste0("LASSO AR BagSize_Test_Sims5_",num,"_CPS"))[,arCPSTuneMode] / merged
    rm(ar)
    
    res.matMode[z,] <- colMeans(modeMat) # average results using Mode as tuning
    #-------------------------
    
    PCR_tableAug <- cbind(PCR_table)
    
    ## box and whisker plot with color and no legend
    Method <- rep(c("OSE", "OSE-Stack",  "OSE-CPS", "SSE", "SSE-Stack",
                    "SSE-CPS", "AR",  "AR-Stack","AR-CPS"),
                  each = 100)
    # Method <- rep(c("TOS", "TOS-CPS", "TOS-Stack", "TSS", "TSS-CPS", "TSS-Stack", 
    #                 "AR","AR-CPS", "AR-Stack"), 
    #               each = 100)
    
    PCR_gg2 <- as.vector(PCR_table)
    PCR_gg2 <- data.frame(RMSE = PCR_gg2, Methods = Method ) 
    colnames(PCR_gg2) <- c("RMSE", "Methods")
    PCR_gg2[,2] <- factor(PCR_gg2[,2])
    
    ensemble <- rep( c("OSE", "SSE", "AR"), each = 100 )
    weights <- rep( c("AVG", "CPS", "Stack"), each = 100 )
    #PCR_gg2 <- cbind(PCR_gg2, ensemble, weights)
    
    pf <- PCR_table %>%
        as_tibble() %>%
        gather(method, rmse) %>%
        mutate(method = factor(method, levels =
                                   c("OSE", "OSE-Stack",  "OSE-CPS", "SSE", "SSE-Stack",
                                     "SSE-CPS", "AR",  "AR-Stack","AR-CPS"))) 
    
    pf$ensemble <- rep( c("Observed-Studies", "Study Strap", "Accept/Reject"), each = 300 )
    pf$Weights  <- rep( 
        rep( c("Average", "Stacking","Covariate Profile Similarity"), 3), 
        each = 100 )
    
    pf$ensemble <- factor(pf$ensemble, levels = c("Observed-Studies", "Study Strap", "Accept/Reject"))
    pf$Weights <- factor(pf$Weights, levels = c("Average", "Stacking","Covariate Profile Similarity"))
    
    
    p <- pf %>%
        ggplot(aes(x = ensemble, y = log(rmse), fill = Weights)) + # , color = Weights
        geom_boxplot(
            lwd = 1.5, 
            fatten = 0.5, 
            alpha = 0.5 #, 0.5
        ) +
        theme_void(base_size = 12) + 
        scale_y_continuous(labels=scaleFUN) +
        guides(col = guide_legend(ncol = 9, title="Weights: ")) + 
        ################################
    guides(colour = guide_legend("Weights:"),
           fill = guide_legend("Weights:")) + 
        ################################
    geom_hline(yintercept=0, 
               linetype="dashed", 
               color = "black", 
               size = rel(0.5),
               alpha = 0.7) + # 
        ggtitle( TeX(sprintf(paste0('$\\mathbf{\\sigma_{\\beta}^2 =$', sigma, 
                                    ', $\\sigma_{\\mathbf{X}}^2 =$', sig.mu, "}$"
        )))) + 
        xlab("Method") + ylab(TeX('$\\mathbf{log(RMSE/RMSE_{TOM})}$')) + 
        scale_color_manual(values = c("#ca0020", "#0868ac", "#E69F00")) + # "#525252",
        scale_fill_manual(values = c("#ca0020", "#0868ac", "#E69F00")) + # "#525252",
        theme(
            axis.text.y=element_text(face="bold",color="black", size=rel(3.25)), #3.25
            axis.text.x=element_text(face="bold",color="black", size=rel(3.4)), # 3.4 GCL
            axis.title = element_text(face="bold", color="black", size=rel(1.5)),
            axis.title.y = element_blank(),
            axis.title.x = element_blank(),
            legend.position = "bottom", 
            legend.key.size = unit(7, "line"), # added in to increase size
            legend.text = element_text(face="bold", color="black", size = rel(3.75)), # 3.4
            legend.title = element_text(face="bold", color="black", size = rel(4)), #3.75
            plot.title = element_text(hjust = 0.5, color="black", size=44, face="bold"),
            legend.spacing.x = unit(5, "mm"),
            legend.background = element_blank(),
            legend.box.background = element_rect(colour = "black", size = 2),
            legend.box.margin = margin(l = 20, r = 20)
            ) +
        scale_x_discrete(labels = c("Observed-Studies", "   Study Strap", "Accept/Reject") ) # x-axis text 3.25

    ##### remove x-axis labels (ensembling scheme)
    if(num == 160 | num == 163 | num == 172 | num == 175){
        print(num)
        p <- p + theme( axis.text.x = element_blank() )
    }
    
    if(num == 171 | num == 163 | num == 173 | num == 175 | num == 183){
        print(num)
        p <- p + theme( axis.text.y = element_blank() )
    }
    
    
    if(num == 172 | num == 180 | num == 160 | num == 168){
        p <- p + theme( plot.margin = margin(l = 30) ) # add margin on left so y axis has room
    }
    
    if(num == 175 | num == 183| num == 163 | num == 171){
        p <- p + theme( plot.margin = margin(r = 100) )  # add margin on right to compensate for y axis and margin added on left above
    }
    
    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    # supplemental figures
    
    assign(paste0("psup", num), 
           p  + theme(axis.text.x=element_text(face="bold",color="black", size=rel(3.5)),
                      axis.text.y=element_text(face="bold",color="black", size=rel(3.25)),
                      plot.margin = margin(l = 15, r = 15, t = 15, b = 15),
                      plot.title = element_text(hjust = 0.5, color="black", size=43, face="bold"),
                      legend.key.size = unit(5, "line"), # added in to increase size
                      legend.text = element_text(face="bold", color="black", size = rel(3.75)), # 3.4
                      legend.title = element_text(face="bold", color="black", size = rel(4))) +
               scale_y_continuous(breaks = pretty(log(pf$rmse), min.n = 3, n = 3), labels=scaleFUN)+   
               scale_x_discrete(labels = c("  Observed- \n Studies", "  Study \n Strap", "  Accept/ \n Reject") ) + # x-axis
               geom_hline(yintercept=0, 
                          linetype="dashed", 
                          color = "black", 
                          size = rel(1),
                          alpha = 0.7) 
               ) #,
                                  # limits= c(pretty(log(pf$rmse), min.n = 3, n = 3)[1],
                                  #           pretty(log(pf$rmse), min.n = 3, n = 3)[length(pretty(log(pf$rmse), min.n = 3, n = 3))])
            

    
    if(num == 168 | num == 169){
        assign(paste0("psup", num), 
               p  + theme(axis.text.x=element_text(face="bold",color="black", size=rel(3.5)),
                          axis.text.y=element_text(face="bold",color="black", size=rel(3.25)),
                          plot.margin = margin(l = 15, r = 15, t = 15, b = 15),
                          plot.title = element_text(hjust = 0.5, color="black", size=43, face="bold"),
                          legend.key.size = unit(5, "line"), # added in to increase size
                          legend.text = element_text(face="bold", color="black", size = rel(3.75)), # 3.4
                          legend.title = element_text(face="bold", color="black", size = rel(4))) +
                   scale_y_continuous(breaks = c(-4,-2,0,2 ), labels=scaleFUN ) +   
                   scale_x_discrete(labels = c("  Observed- \n Studies", "  Study \n Strap", "  Accept/ \n Reject") ) # x-axis
        )
        
        
    }
    
    # change y axis ticks
    if(num == 173)  assign(paste0("psup", num), 
                           get(paste0("psup", num)) + 
                               scale_y_continuous(breaks = c(0,-1, -2 ), 
                                                  labels=scaleFUN,
                                                  limits = c(-2.1, 0.6)) )  


    # remove x axis since these arent bottom row
if(is.element(num, c(160,161,164,165,168,169,
                     172,173,176,177,180,181) ) ){
    assign(paste0("psup", num), 
           get(paste0("psup", num))  + theme(axis.text.x=element_blank())
    )
                # x-axis
    
    
    
}
    
    legend2 <- get_legend(get(paste0("psup", num)))
    assign(paste0("psup", num), get(paste0("psup", num)) + theme(legend.position="none") )
    
    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    # 2. Save the legend
    #+++++++++++++++++++++++
    legend <- get_legend(p)
    # 3. Remove the legend from the box plot
    #+++++++++++++++++++++++
    p <- p + theme(legend.position="none")
    
    assign(paste0("p", num), p) # rename
    #-------------------------------------------
    #+++++++++++++++++++++++
    assign(paste0("pa", num), p +
               scale_y_continuous(breaks = c(-0.2, 0, 0.2), limits = c(-.35, 0.35) ) )
    assign(paste0("pl", num), p + 
               scale_y_continuous(breaks = round(c(-4.0, 0, 4.0),2), limits = c(-5.0, 5.0), labels=scaleFUN ) )
    
    
    assign(paste0("ps", num), p + 
               scale_y_continuous(breaks = round(c(-0.5, 0, 0.5),2), limits = c(-1, 0.5), labels=scaleFUN ) )
    
    assign(paste0("pq", num), p + 
               scale_y_continuous(breaks = round(c(-3.0, 0, 3.0),2), limits = c(-3.5, 4), labels=scaleFUN ) ) #+ coord_cartesian(ylim = c(-3.5, 4)) )
    #####################################################
        
}

bagSizeComp <- cbind(sim.matrix,tuneComparison)

cbind(res.mat,testVsOpt)

# save file

# main paper figures
# clusters
ggsave("Sims_Box_sims5_ClustersLog_Legend_reduced2.pdf", 
       plot = grid.arrange( 
           grobs = list(ps172, pq180, 
                        ps175, pq183, legend),  
           layout_matrix = rbind(matrix(seq(1,9), ncol = 2, nrow = 2), rep(10, 3) ),
           nrow = 3, ncol = 2,
           heights = c(2.3, 2.3, 1),
           left = textGrob(TeX('$\\mathbf{log(RMSE/RMSE_{TOM})}$'), rot = 90,
                           vjust = 0.43,
                           hjust = 0.3,
                           gp = gpar(fontface = "bold", cex = 4.25))),
       width = 32, height = 15)

# no clusters
ggsave("Sims_Box_sims5_HybridShiftLog_Legend_reduced2.pdf", 
       plot = grid.arrange(grobs = list(pa160, pl168,
                                        pa163,   pl171,
                                        legend),
                           layout_matrix = rbind(matrix(seq(1,5), ncol = 2, nrow = 2), rep(10, 3) ),
                           nrow = 3, ncol = 2,
                           heights = c(2.3, 2.3, 1), #2.3, 2.3, 1
                           left = textGrob(TeX('$\\mathbf{log(RMSE/RMSE_{TOM})}$'), rot = 90,
                                           vjust = 0.43,
                                           hjust = 0.3,
                                           gp = gpar(fontface = "bold", cex = 4.25))), # 3.5
       width = 32, height = 15)


#+++++++++++++++++++++++++++++++++++++++++++
# simulation figures

ggsave("Sims_Box_sims5_Clusters_Legend2.png", 
       plot = grid.arrange( 
           grobs = list(psup172, psup173,  psup175,# pz174,
                        psup176, psup177,  psup179, # po178,
                        psup180, psup181, psup183, #po182, 
                        legend2), #po182, 
           layout_matrix = rbind(matrix(seq(1,9), ncol = 3, nrow = 3), rep(10, 3) ),
           nrow = 4, ncol = 3,
           heights = c(2.3, 2.3, 2.3, 0.75), #2.3, 2.3, 1
           left = textGrob(TeX('$\\mathbf{log(RMSE/RMSE_{TOM})}$'), rot = 90,
                           vjust = 0.43,
                           hjust = 0.3,
                           gp = gpar(fontface = "bold", cex = 4.5))),
       width = 32, height = 20)


ggsave("Sims_Box_sims5_HybridShift_Legend2.png", 
       plot = grid.arrange(grobs = list(psup160, psup161,  psup163, #p162,
                                        psup164, psup165, psup167, #p166,
                                        psup168, psup169,  psup171,
                                        legend2), #po182, 
                           layout_matrix = rbind(matrix(seq(1,9), ncol = 3, nrow = 3), rep(10, 3) ),
                           nrow = 4, ncol = 3,
                           heights = c(2.3, 2.3, 2.3, 0.75), #2.3, 2.3, 1
                           left = textGrob(TeX('$\\mathbf{log(RMSE/RMSE_{TOM})}$'), rot = 90,
                                           vjust = 0.43,
                                           hjust = 0.3,
                                           gp = gpar(fontface = "bold", cex = 4.5))),
       width = 32, height = 20)

