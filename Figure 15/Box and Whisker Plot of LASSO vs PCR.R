# figure dimensions
# Width: 750, Height 550

library(ggplot2)
library(gridExtra)
library(viridis)

setwd("~/Desktop/Research Final/PCR Merged")
PCR_merged <- read.csv("Merged PCR Expedited Tuning 0 Trunc_seed_y_samps_ 1_15_ncomp_28_tune.met_5_tune.frac_1_samp_size_full")[,4]

setwd("~/Desktop/Research Final/PCR Diff Merged")
pcr.diff <- read.csv("Merged PCR Deriv_Expedited Tuning 0 Trunc_seed_y_samps_ 1_15_ncomp_28_tune.met_5_tune.frac_1_samp_size_full")[,4]

setwd("~/Desktop/Harvard Biostatistics/FSCV with Ken/Final Scripts/Merged Full/LASSO DIff/Tuned 4-15-19")
lasso.diff <- read.csv("Merged LASSO Deriv 0 Trunc_seed_y_samps_ 1_15_ncomp_28_tune.met_5_tune.frac_1_sampSize_full")[,4]

setwd("~/Desktop/Harvard Biostatistics/FSCV with Ken/Final Scripts/Merged Full/LASSO")
lasso <- read.csv("Merged LASSO 0 Trunc_seed_y_samps_ 1_15_ncomp_28_tune.met_5_tune.frac_1_sampSize_full")
lasso <- lasso$RMSE

PCR_table <- cbind(PCR_merged, lasso, pcr.diff, lasso.diff)
colnames(PCR_table) <- c("PCR", "LASSO", "PCR Derivative", "LASSO Derivative")

Method <- rep(c("PCR", "LASSO", "PCR Derivative", "LASSO Derivative"), each = 15)
PCR_gg2 <- c()
for (i in 1:ncol(PCR_table)){
  PCR_gg2 <- c(PCR_gg2, PCR_table[,i])
}
PCR_gg2 <- as.data.frame(cbind(PCR_gg2, Method) )
colnames(PCR_gg2) <- c("RMSE", "Methods")
PCR_gg2$Methods <- factor(PCR_gg2$Methods, levels = c("LASSO", "PCR", "LASSO Derivative", "PCR Derivative"))

reg <- rep( rep( c("PCR", "LASSO"), each = 15), 2)
cov <- rep( c("Raw", "Derivative"), each = 30)

PCR_gg2 <- data.frame( RMSE = as.numeric(as.vector(PCR_gg2$RMSE)), 
                       Regression = as.factor(reg), 
                       Covariates = as.factor(cov) )

PCR_gg2$Covariates <- factor(PCR_gg2$Covariates, levels = c("Raw", "Derivative"), ordered = TRUE)

# ggplot
p <- ggplot(PCR_gg2, aes(x = Covariates, y = RMSE, 
                         fill = Regression)) + 
    geom_boxplot(lwd = 1, fatten = 0.5)  +
    scale_color_manual(values=c("#ca0020", "#0868ac")) +
    scale_fill_manual(values=c("#ca0020", "#0868ac")) +
    labs(x = "Covariates")  +
    theme(
        legend.text =  element_text(face="bold", color="black", size=rel(1.5)),
        legend.title = element_text(face="bold", color="black", size=rel(2)),
        axis.title.x = element_text(face="bold", color="black", 
                            size = rel(2),
                            vjust = -2),
        axis.title.y = element_text(face="bold", color="black", 
                            size = rel(2),
                            vjust = 2.5),
        axis.text = element_text(face="bold", color="black", 
                                                 size = rel(1.5)),
        plot.margin = margin(r = 20, b = 20, l = 15, t = 5) ) 

setwd("~/Desktop/Research Final/Methods Paper/Draft/Figures") # setwd for saving file
ggsave("PCR_vs_lasso.png", 
       plot = p,
       width = 7, height = 4)
