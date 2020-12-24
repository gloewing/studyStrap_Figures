library(ggplot2)
library(gridExtra)
library(viridis)

####################################
# Original (Raw) Covariates
####################################
setwd("~/Desktop/Research Final/PCR Merged")

#Our transformation function to ensure rounding is same decimal places for all
scaleFUN <- function(x) sprintf("%.1f", x) # sprintf("%.2f", x)


merged <- read.csv("Merged PCR Expedited Tuning 0 Trunc_seed_y_samps_ 1_15_ncomp_28_tune.met_5_tune.frac_1_samp_size_full")
merged <- merged$RMSE

setwd("~/Desktop/Research")
ss <- read.csv("Study Strap PCR Opt Res_seed_y_samp_ 15_strap_500_bgSz_14sampSzfull")
cpsw <- ss$Weighted
ssStack <- ss$Stacking
ss <- ss$Mean

setwd("~/Desktop/Research Final/Naive Ensemble PCR")
ens <- read.csv("Directed Stacked PCR Results 0 Truncate (7-27-19)y_tuneMet_6_smpSz_2500_subSampInd_0")
stack <- ens$Stacked.Int
ens.cpsw <- ens$Wiggle.Distance
ens<- ens$Mean

########################################
setwd("~/Desktop/Research/AR_multiseed")

# labels and tables
get_legend<-function(myggplot){
    tmp <- ggplot_gtable(ggplot_build(myggplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
}


a1 <- "SS PCR AR Repeat Opt Res_seed_y+"
a2<- "_samp_ 15_ncomp_20_LR_1_strap_200_conv_3e+05_tuneMet_6_tnFrc_1_seqLen_11_bag.sel_1_sprd_-1_coord_0_svBagMat_0_ssTune_10_Rpts"

a <- read.csv("SS PCR AR Repeat Opt Res_seed_y_samp_ 15_ncomp_20_LR_1_strap_250_conv_3e+05_tuneMet_6_tnFrc_1_seqLen_11_bag.sel_1_sprd_-1_coord_0_svBagMat_0_ssTune_10_Rpts")[,2:3]/merged

l <- seq(10,90, by = 10) # change to 100
num.elecs <- 15
num.seeds <- length(l)
dat <- matrix(nrow = num.elecs, ncol = num.seeds)
dat.weight <- matrix(nrow = num.elecs, ncol = num.seeds)
dat.stack <- matrix(nrow = num.elecs, ncol = num.seeds)

# iterate through the rest of the seeds
for (i in 1:length(l)){
    a <- read.csv(paste0(a1,l[i],a2))/merged # read in and scale by merged
    dat[,i ] <- a$Mean # add for simple average across straps
    dat.weight[,i ] <- a$Weighted # CPSW weights
    dat.stack[,i] <- a$Stacking
    
}

#avg across seeds
ar <- rowMeans(dat) #avg across seeds
ar.cpsw <- rowMeans(dat.weight)
ar.stack <- rowMeans(dat.stack)

PCR_table <- cbind(ens / merged,  stack / merged, ens.cpsw / merged, 
                   ss / merged, ssStack / merged, cpsw / merged, ar, ar.stack,  ar.cpsw ) #AR already scaled by merged
colnames(PCR_table) <- c("OSE", "OSE-Stack",  "OSE-CPS", "SSE", "SSE-Stack",
                         "SSE-CPS", "AR", "AR-Stack","AR-CPS")
merge.mean <- mean(merged)
PCR_table.raw <- PCR_table
# scale by mean of merged
## box and whisker plot with color and no legend
Method <- rep(c("OSE", "OSE-Stack",  "OSE-CPS", "SSE", "SSE_Stack",
                "SSE-CPS", "AR", "AR_Stack","AR-CPS"), each = 15)
Method = factor(Method, levels = 
                    c("OSE", "OSE-Stack",  "OSE-CPS", "SSE", "SSE-Stack",
                      "SSE-CPS", "AR", "AR-Stack","AR-CPS") )

pf <- PCR_table %>%
    as_tibble() %>%
    gather(method, rmse) %>%
    mutate(method = factor(method, levels =
                               c("OSE", "OSE-Stack",  "OSE-CPS", "SSE", "SSE-Stack",
                                 "SSE-CPS", "AR",  "AR-Stack","AR-CPS"))) 


pf$ensemble <- rep( c("Observed-Studies", "Study Strap", "Accept/Reject"), each = 45 )
pf$Weights  <- rep( 
    rep( c("Average", "Stacking","Covariate Profile Similarity"), 3), 
    each = 15 )

pf$ensemble <- factor(pf$ensemble, levels = c("Observed-Studies", "Study Strap", "Accept/Reject"))
pf$Weights <- factor(pf$Weights, levels = c("Average", "Stacking","Covariate Profile Similarity"))

p <- ggplot(pf, aes(x = ensemble, y = log(rmse), fill = Weights)) + # , color = Weights
    geom_boxplot(
        lwd = 1.5, 
        fatten = 0.5, 
        alpha = 0.5 
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
    ggtitle( TeX(sprintf(paste0('$\\mathbf{Raw}$')))) + 
    xlab("Method") + ylab(TeX('$\\mathbf{log(RMSE/RMSE_{TOM})}$')) + 
    scale_color_manual(values = c("#ca0020", "#0868ac", "#E69F00")) + # "#525252",
    scale_fill_manual(values = c("#ca0020", "#0868ac", "#E69F00")) + # "#525252",
    theme(
        axis.text.y=element_text(face="bold",color="black", size=rel(1.75)), #3.25
        axis.text.x=element_text(face="bold",color="black", size=rel(1.75)), # 3.25 GCL
        axis.title = element_text(face="bold", color="black", size=rel(1.5)),
        axis.title.y = element_blank(), #element_text(face="bold", color="black", size=rel(1.25)),
        axis.title.x = element_blank(),
        legend.position = "bottom", 
        legend.key.size = unit(3, "line"), # added in to increase size
        legend.text = element_text(face="bold", color="black", size = rel(2)), # 3 GCL
        legend.title = element_text(face="bold", color="black", size = rel(2)),
        plot.title = element_text(hjust = 0.5, color="black", size=30, face="bold"),
        legend.spacing.x = unit(5, "mm"),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black", size = 1.5),
        legend.box.margin = margin(l = 20, r = 20)
    ) +
    scale_x_discrete(labels = c("Observed-Studies", "   Study Strap", "Accept/Reject") ) # x-axis text 3.25
# scale_x_discrete(labels = c("  Observed- \n Studies", "  Study \n Strap", "  Accept/ \n Reject") ) # x-axis text 3.5


# 2. Save the legend
#+++++++++++++++++++++++
legend <- get_legend(p)
# 3. Remove the legend from the box plot
#+++++++++++++++++++++++
p <- p + theme(legend.position="none")
    
assign(paste0("raw"), p + 
           theme( plot.margin = margin(l = 20) ) +
           scale_y_continuous(labels=scaleFUN, 
                              breaks = c(-0.5, 0.0, 0.5, 1.0, 1.5, 2.0), 
                              limits = c(-0.5, 2.02))  
        ) # unscaled
#-------------------------------------------------------------------------------

####################################
# Derivative Covariates
####################################
setwd("~/Desktop/Research Final/PCR Diff Merged")

#Our transformation function to ensure rounding is same decimal places for all
scaleFUN <- function(x) sprintf("%.1f", x) # sprintf("%.2f", x)

merged <- read.csv("Merged PCR Deriv_Expedited Tuning 0 Trunc_seed_y_samps_ 1_15_ncomp_28_tune.met_5_tune.frac_1_samp_size_full")
merged <- merged$RMSE
merge.mean <- mean(merged)

setwd("~/Desktop/Research")
ss <- read.csv("Study Strap PCR Derivative Opt Res_seed_y_samp_ 15_strap_500_bgSz_75sampSzfull")
cpsw <- ss$Weighted
ssStack <- ss$Stacking
ss <- ss$Mean

setwd("~/Desktop/Research Final/Naive Ensemble PCR Diff")
ens <- read.csv("Stacked PCR Diff Results 0 Truncate (8-17)y_tuneMet_6_smpSz_2500_subSampInd_0")
stack <- ens$Stacked.Int
ens.cpsw <- ens$Mode
ens<- ens$Mean


########################################
setwd("~/Desktop/Research/AR_multiseed")

# labels and tables
get_legend<-function(myggplot){
    tmp <- ggplot_gtable(ggplot_build(myggplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
}


a1 <- "SS PCR AR LR Derivative Opt Repeated Res_seed_y+"
a2<- "_samp_ 15_ncomp_20_LR_1_strap_200_conv_3e+05_tuneMet_0_tuneFrac_1_seqLen_18_smpSz_full_bgSz_75_Rpts_10"

l <- seq(10, 100, by = 10) 
num.elecs <- 15
num.seeds <- length(l)
dat <- matrix(nrow = num.elecs, ncol = num.seeds)
dat.weight <- matrix(nrow = num.elecs, ncol = num.seeds)
dat.stack <- matrix(nrow = num.elecs, ncol = num.seeds)

# iterate through the rest of the seeds
for (i in 1:length(l)){
    a <- read.csv(paste0(a1,l[i],a2))/merged # read in and scale by merged
    dat[,i ] <- a$Mean # add for simple average across straps
    dat.weight[,i ] <- a$Weighted # CPSW weights
    dat.stack[,i] <- a$Stacking
    
}

#avg across seeds
ar <- rowMeans(dat) #avg across seeds
ar.cpsw <- rowMeans(dat.weight)
ar.stack <- rowMeans(dat.stack)

PCR_table <- cbind(ens / merged,  stack / merged, ens.cpsw / merged, 
                   ss / merged, ssStack / merged, cpsw / merged, ar, ar.stack,  ar.cpsw ) #AR already scaled by merged
colnames(PCR_table) <- c("OSE", "OSE-Stack",  "OSE-CPS", "SSE", "SSE-Stack",
                         "SSE-CPS", "AR", "AR-Stack","AR-CPS")
merge.mean <- mean(merged)
PCR_table.deriv <- PCR_table
# scale by mean of merged
## box and whisker plot with color and no legend
Method <- rep(c("OSE", "OSE-Stack",  "OSE-CPS", "SSE", "SSE_Stack",
                "SSE-CPS", "AR", "AR_Stack","AR-CPS"), each = 15)
Method = factor(Method, levels = 
                    c("OSE", "OSE-Stack",  "OSE-CPS", "SSE", "SSE-Stack",
                      "SSE-CPS", "AR", "AR-Stack","AR-CPS") )

pf <- PCR_table %>%
    as_tibble() %>%
    gather(method, rmse) %>%
    mutate(method = factor(method, levels =
                               c("OSE", "OSE-Stack",  "OSE-CPS", "SSE", "SSE-Stack",
                                 "SSE-CPS", "AR",  "AR-Stack","AR-CPS"))) 


pf$ensemble <- rep( c("Observed-Studies", "Study Strap", "Accept/Reject"), each = 45 )
pf$Weights  <- rep( 
    rep( c("Average", "Stacking","Covariate Profile Similarity"), 3), 
    each = 15 )

pf$ensemble <- factor(pf$ensemble, levels = c("Observed-Studies", "Study Strap", "Accept/Reject"))
pf$Weights <- factor(pf$Weights, levels = c("Average", "Stacking","Covariate Profile Similarity"))

p <- ggplot(pf, aes(x = ensemble, y = log(rmse), fill = Weights)) + # , color = Weights
    geom_boxplot(
        lwd = 1.5, 
        fatten = 0.5, 
        alpha = 0.5 
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
    ggtitle( TeX(sprintf(paste0('$\\mathbf{Derivative}$')))) + 
    xlab("Method") + ylab(TeX('$\\mathbf{log(RMSE/RMSE_{TOM})}$')) + 
    scale_color_manual(values = c("#ca0020", "#0868ac", "#E69F00")) + # "#525252",
    scale_fill_manual(values = c("#ca0020", "#0868ac", "#E69F00")) + # "#525252",
    theme(
        axis.text.y=element_text(face="bold",color="black", size=rel(1.75)), #3.25
        axis.text.x=element_text(face="bold",color="black", size=rel(1.75)), # 3.25 GCL
        axis.title = element_text(face="bold", color="black", size=rel(1.5)),
        axis.title.y = element_blank(), #element_text(face="bold", color="black", size=rel(1.25)),
        axis.title.x = element_blank(),
        legend.position = "bottom", 
        legend.key.size = unit(3, "line"), # added in to increase size
        legend.text = element_text(face="bold", color="black", size = rel(1.75)), # 3 GCL
        legend.title = element_text(face="bold", color="black", size = rel(2)),
        plot.title = element_text(hjust = 0.5, color="black", size=30, face="bold"),
        legend.spacing.x = unit(5, "mm"),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black", size = 1),
        legend.box.margin = margin(l = 10, r = 10, b = 0)
    ) +
    scale_x_discrete(labels = c("Observed-Studies", "   Study Strap", "Accept/Reject") ) # x-axis text 3.25
# scale_x_discrete(labels = c("  Observed- \n Studies", "  Study \n Strap", "  Accept/ \n Reject") ) # x-axis text 3.5

# change for derivative only because it is on the right side
p <- p + theme( 
                axis.text.y = element_blank() )   # add margin on right to compensate for y axis and margin added on left above


# 2. Save the legend
#+++++++++++++++++++++++
legend <- get_legend(p)
# 3. Remove the legend from the box plot
#+++++++++++++++++++++++
p <- p + theme(legend.position="none")

assign(paste0("deriv"), p + 
           theme( plot.margin = margin(r = 60) ) +
           scale_y_continuous(labels=scaleFUN, 
                              breaks = c(-0.5, 0.0, 0.5, 1.0, 1.5), 
                              limits = c(-0.5, 2.02))   #-0.5, 2.1
) # unscaled

#######################
# Figures
#######################
setwd("~/Desktop/Research Final/Methods Paper/Draft/Figures") # setwd for saving file

ggsave("Both_Methods_Final.pdf", 
       plot = grid.arrange(grobs = list(raw, deriv, legend),
                           layout_matrix = rbind(c(1,2), c(3,3) ),
                           nrow = 2, ncol = 2,
                           heights = c(2.3, 0.75), # 4.6, 1 or 2.3, 0.5
       left = textGrob(TeX('$\\mathbf{log(RMSE/RMSE_{TOM})}$'), rot = 90,
                       hjust = 0.3,
                       vjust = 0.75,
                       gp = gpar(fontface = "bold", cex = 1.85))),
       width = 16, height = 5)

############################################################################
# different sizing
ggsave("Both_Methods_Final.pdf", 
       plot = grid.arrange(grobs = list(raw, deriv, legend),
                           layout_matrix = rbind(c(1,2), c(3,3) ),
                           nrow = 2, ncol = 2,
                           heights = c(4.5, 0.75), # 4.6, 1 or 2.3, 0.5; 7, 1.75
                           left = textGrob(TeX('$\\mathbf{log(RMSE/RMSE_{TOM})}$'), rot = 90,
                                           hjust = 0.45,
                                           vjust = 0.75,
                                           gp = gpar(fontface = "bold", cex = 2.25))),
       width = 16, height = 8)


# #####################################
# # Tables
# #####################################
library(kableExtra)
PCR_table <- round( rbind(colMeans(PCR_table.raw), colMeans(PCR_table.deriv) ), 3)
row.names(PCR_table) <- c("Raw", "Derivative")
kable((PCR_table),  format = "latex", booktabs = T) %>% kable_styling(position = "center")
