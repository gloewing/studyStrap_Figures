library(ggplot2)
library(latex2exp)
library(grid)
library(gridExtra)
library(viridis)

################
# Raw Data
################
setwd("~/Desktop/Research")
merged <- read.csv("Merged PCR Expedited Tuning 0 Trunc_seed_y_samps_ 1_15_ncomp_28_tune.met_5_tune.frac_1_samp_size_full")[,4]

# load data
setwd("~/Desktop/Research/AR_multiseed")
a1 <- "SS PCR AR Repeat Opt Res_seed_y+"
a2<- "_samp_ 15_ncomp_20_LR_1_strap_250_conv_3e+05_tuneMet_6_tnFrc_1_seqLen_11_bag.sel_1_sprd_-1_coord_0_svBagMat_0_ssTune_10_Rpts"

l <- seq(10,90, by = 10)

dat <- matrix(nrow = length(l), ncol = 2)
colnames(dat) <- c("Average", "CPS")

for (i in 1:length(l)){
    a <- read.csv(paste0(a1,l[i],a2))[,2:3]/merged
    
    dat[i,] <- c(colMeans(a))

}

a <- read.csv("SS PCR AR Repeat Opt Res_seed_y_samp_ 15_ncomp_20_LR_1_strap_250_conv_3e+05_tuneMet_6_tnFrc_1_seqLen_11_bag.sel_1_sprd_-1_coord_0_svBagMat_0_ssTune_10_Rpts")[,2:3]
dat <- rbind(colMeans(a/merged), dat)
setwd("~/Desktop/Research")
merged <- read.csv("Merged PCR Expedited Tuning 0 Trunc_seed_y_samps_ 1_15_ncomp_28_tune.met_5_tune.frac_1_samp_size_full")

dat <- data.frame(dat)
dat.raw <- c(dat[,1], dat[,2])
dat.raw <- data.frame( Mean = dat.raw, Method = c(rep("Average", 10), rep("CPS", 10))    )
dat.raw$Method <- factor(dat.raw$Method, levels = c("Average", "CPS"))

pp.raw <- ggplot(dat.raw, aes(y = Mean, x = Method, legend = NA, color = factor(Method))) +
    geom_boxplot(lwd = 1, fatten = 0.5) + 
    geom_hline(yintercept= 1, linetype="dashed", color = "red") + 
    # scale_x_discrete("Average", "CPS") +
    ylab(TeX('$\\mathbf{RMSE/RMSE_{TOM}}$')) + 
    xlab("Weights") +
    scale_color_manual(values=c("darkred", "darkblue")) +
    ylim(.86,1) + 
    theme(axis.text.x=element_text(face="bold", color="black", size = rel(2)),
          axis.text.y=element_text(face="bold",color="black", size = rel(2)),
          strip.text = element_text(size = 10,face="bold",color="black"), 
          axis.title.x = element_text(face="bold", color="black", 
                                      size = rel(2),
                                      vjust = -2),
          axis.title.y = element_text(face="bold", color="black", 
                                      size = rel(2),
                                      vjust = 2.5),
          legend.position = "none", 
          plot.margin = margin(r = 20, b = 20, l = 15, t = 5) ) 


################
# Derivative 
################

setwd("~/Desktop/Research")
merged <- read.csv("Merged PCR Deriv_Expedited Tuning 0 Trunc_seed_y_samps_ 1_15_ncomp_28_tune.met_5_tune.frac_1_samp_size_full")[,4]

# load data
setwd("~/Desktop/Research/AR_multiseed")
a1 <- "SS PCR AR LR Derivative Opt Repeated Res_seed_y+"
a2<- "_samp_ 15_ncomp_20_LR_1_strap_250_conv_3e+05_tuneMet_5_tuneFrac_1_seqLen_18_smpSz_full_bgSz_75_Rpts_10"

# individual Box plots
a <- read.csv("SS PCR AR LR Derivative Opt Repeated Res_seed_y_samp_ 15_ncomp_20_LR_1_strap_250_conv_3e+05_tuneMet_5_tuneFrac_1_seqLen_18_smpSz_full_bgSz_75_Rpts_10")[,2:3]

l <- seq(10, 90, by = 10)#c(60:63, 65:71)#60:71

dat <- matrix(nrow = length(l), ncol = 2)
colnames(dat) <- c("Average", "CPS")
for (i in 1:length(l)){
    a <- read.csv(paste0(a1,l[i],a2))[,2:3]/merged
    dat[i,] <- c(colMeans(a))

}
a <- read.csv("SS PCR AR LR Derivative Opt Repeated Res_seed_y_samp_ 15_ncomp_20_LR_1_strap_250_conv_3e+05_tuneMet_5_tuneFrac_1_seqLen_18_smpSz_full_bgSz_75_Rpts_10")[,2:3]
dat <- rbind(colMeans(a/merged), dat)
setwd("~/Desktop/Research")
merged <- read.csv("Merged PCR Deriv_Expedited Tuning 0 Trunc_seed_y_samps_ 1_15_ncomp_28_tune.met_5_tune.frac_1_samp_size_full")[,4]


dat <- data.frame(dat)
dat.deriv <- c(dat[,1], dat[,2])
dat.deriv <- data.frame( Mean = dat.deriv, Method = c(rep("Average", 10), rep("CPS", 10))    )
dat.deriv$Method <- factor(dat.deriv$Method, levels = c("Average", "CPS"))


pp.deriv <- ggplot(dat.deriv, aes(y = Mean, x = Method, legend = NA, color = factor(Method))) +
    geom_boxplot(lwd = 1, fatten = 0.5) + 
    geom_hline(yintercept= 1, linetype="dashed", color = "red") + 
    ylab(TeX('$\\mathbf{RMSE/RMSE_{TOM}}$')) + 
    xlab("Weights") +
    scale_color_manual(values=c("darkred", "darkblue")) +
    ylim(.86,1) + 
    theme(axis.text.x=element_text(face="bold", color="black", size = rel(2)),
          axis.text.y=element_text(face="bold",color="black", size = rel(2)),
          strip.text = element_text(size = 10,face="bold",color="black"), 
          axis.title.x = element_text(face="bold", color="black", 
                                      size = rel(2),
                                      vjust = -2),
          axis.title.y = element_text(face="bold", color="black", 
                                      size = rel(2),
                                      vjust = 2.5),
          legend.position = "none", 
          plot.margin = margin(r = 20, b = 20, l = 15, t = 5) ) 

setwd("~/Desktop/Research Final/Methods Paper/Draft/Figures") # setwd for saving file
ggsave("Bw_seed_raw.png", 
       plot = pp.raw,
       width = 5, height = 4)
ggsave("Bw_seed_deriv.png",
       plot = pp.deriv,
       width = 5, height = 4)

#######################################################################
