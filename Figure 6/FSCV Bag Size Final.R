##########################################################
# Raw and Derivative FSCV Bag Size Tuning Curve
##########################################################

##################
# Derivative
##################
library(ggplot2)
library(latex2exp)
library(grid)
library(gridExtra)
setwd("~/Desktop/Research Final/Study Strap AR/Diff/SubSample Bag Size Curve")

nam <- "SS PCR AR LR Derivative Opt Res_seed_y_samp_ 15_ncomp_20_LR_1_strap_75_conv_3e+05_tuneMet_5_tuneFrac_1_seqLen_11_smpSz_2500_bgSz_"

bagSzs <- c(1, 2,4,6,8,10,12,14,16,18,24,28,35,56,65,75,85,100,125,250, 1000)
bags <- c()
res <- c()
merged <- read.csv("Merged PCR Deriv_Expedited Tuning 0 Trunc_seed_y_samps_ 1_15_ncomp_28_tune.met_6_tune.frac_1_samp_size_2500")
merged.mean <- merged$RMSE
for(sz in bagSzs){
    a <- read.csv(paste0(nam, sz))
    res <- c(res, mean(a$Mean / merged.mean))
    bags <- c(bags, (sz))
    print(mean(a$Mean))
}
naive <- read.csv("Stacked PCR Diff Results 0 Truncate (8-17)y_tuneMet_6_smpSz_2500_subSampInd_1")

bags <- factor(bags, levels = c(bagSzs))
df <- data.frame(RMSE = res, Bag.Size = bags, Derivative = "Derivative")

# minDeriv <- c(1, bagSzs)[which.min(res)]
minDeriv <- c(bagSzs)[which.min(res)]

##################
# Raw
##################
setwd("~/Desktop/Research Final/Study Strap AR/Raw/Subsample 2500 Bag Size")
merged <- read.csv("Merged PCR Expedited Tuning 0 Trunc_seed_y_samps_ 1_15_ncomp_28_tune.met_6_tune.frac_1_samp_size_2500")
merged.mean <- merged$RMSE
nam <- "SS PCR AR LR SubSamp Opt Res_seed_y_samp_ 15_strap_75_conv_3e+05_tuneMet_6_tuneFrac_1_seqLen_11_bag.sel_1_sprd_-2.5_coord_27_samSiz_2500_bwBag_0_bagSt_75_CF_400_lamb_0_VW_1_alpha_0.9992err.typ_2_strp.vec_0_bgSz_"
bagSzs <- c(1,2,4,6,8,10,12,14,16,18,24,28,35,56,65,75,85,100,125,250, 1000)
bags <- c()
res <- c()
for(sz in bagSzs){
    a <- read.csv(paste0(nam, sz))
    res <- c(res, mean(a$Mean / merged.mean))
    bags <- c(bags, (sz))
    print(mean(a$Mean))
}
naive <- read.csv("Stacked PCR Results 0 Truncate (8-17)y_tuneMet_6_smpSz_2500_subSampInd_1")

bags <- factor(bags, levels = c(bagSzs))#levels = c(1, bagSzs))

df2 <- data.frame(RMSE = res, Bag.Size = bags, Derivative = "Raw")
df <- rbind(df, df2)
minRaw <- c( bagSzs)[which.min(res)]

df$Derivative <- factor(df$Derivative, levels = c("Raw", "Derivative"))
##################
# plot
##################
p <- ggplot(data = df, aes(x=factor(Bag.Size), y=log(RMSE), color = factor(Derivative),
                           fill = factor(Derivative), group = factor(Derivative)) ) + 
    geom_line( size = 2) + 
    geom_point(size = 4, aes(shape = factor(Derivative ))) + 
    theme_classic(base_size = 12) + 
    xlab("Bag Size") +
    scale_x_discrete(breaks = df$Bag.Size[seq(1,21, 2)] )  +# c(1, 4, 8, 14, 24,35, 75, 100, 1000))  +
    theme(legend.position = "none", 
          plot.title = element_text(hjust = 0.5, size=14, face="bold")) +
    ylab(TeX('$\\mathbf{log(RMSE/RMSE_{TOM})}$')) +
    scale_color_manual(values = c("#ca0020", "#0868ac")) +
    theme(axis.text.x = element_text(face="bold", # color="#993333", 
                                     angle=45, size=rel(2),
                                     vjust = 0.6 ),  
          axis.text.y = element_text(face="bold", size=rel(2) ), 
          axis.title.y = element_text(face="bold", size=rel(2),
                                      vjust = 3), 
          axis.title.x = element_text(face="bold", size=rel(2) )
          ) +
        annotate("text", x=16, y=1.5,  
             label='bold("Derivative")', 
             size=10, color="#0868ac", parse=TRUE) + 
    annotate("text", x=16, y=1.85,  
             label='bold("Raw")', 
             size=10, color="#ca0020", parse=TRUE) + 
    scale_fill_manual(values=c("#ca0020", "#0868ac")) +
    geom_hline(yintercept=0, 
               linetype="dashed", 
               color = "black", 
               size = rel(0.5),
               alpha = 0.7)  
    

p <- p 
    
##################
# zoom
##################
p2 <- p + coord_cartesian(ylim = c(-0.15, 0.2)) + 
        geom_vline(xintercept = which(df$Bag.Size == minRaw), color = "#ca0020" ) +
            geom_vline(xintercept = which(df$Bag.Size == minDeriv), color = "#0868ac") +
                theme( axis.title.y = element_blank() ) + 
    theme( plot.margin = margin(r = 30, b = 15, l = 15) )  # add margins to counterbalance other plot having space on left for y axis title

p <- p + theme( plot.margin = margin(l = 15, r = 25, b = 15) )  # add margins for y-axis title space


assign(paste0("bag.size.plot"), p ) # unscaled
assign(paste0("bag.size.plot.zoom"), p2  ) # zoom

##################
# save
##################

setwd("~/Desktop/Research Final/Methods Paper/Draft/Figures") # setwd for saving file
ggsave("bagSize_plotLog.pdf", 
       plot = bag.size.plot,
       width = 7, height = 4.5)
ggsave("bagSize_plot_zoomLog.pdf",
       plot = bag.size.plot.zoom,
       width = 7, height = 4.5)
