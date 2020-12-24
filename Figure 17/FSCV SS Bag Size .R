##########################################################
# Raw and Derivative FSCV SS Bag Size Curve
##########################################################

##################
# functions
##################
library(ggplot2)
library(latex2exp)
library(grid)
library(gridExtra)
library(tidyverse)
setwd("~/Desktop/Research/fscv_bagSize")

#Our transformation function to ensure rounding is same decimal places for all
scaleFUN <- function(x) sprintf("%.2f", x)

##################
# Raw
##################
merged <- read.csv("FSCV_raw_150_straps_3_paths_convg10000_bagSize_2500")
merged.mean <- merged$mrg
bagSzs <- c(1,2,4,6,8,10,12,14,16,18,24,28,35,56,65,75,85,100,125,250, 1000)
bags <- c()
res <- c()

avg <- colMeans( read.csv("fscv_raw SS BagSize_Test__Avg") / merged.mean )
stack <- colMeans( read.csv("fscv_raw SS BagSize_Test__Stack") / merged.mean )
cps <- colMeans( read.csv("fscv_raw SS BagSize_Test__CPS") / merged.mean )

names(avg) <- names(stack) <- names(cps) <- bagSzs

df <- rbind(
    cbind(bagSzs, (avg), weights = "Average"),
    cbind(bagSzs, (stack), weights = "Stacking"),
    cbind(bagSzs, (cps), weights = "CPS")
)
df <- data.frame(df)
df <- cbind(df, "raw")
colnames(df) <- c("bag", "rmse", "weights", "raw")


df2 <- df

##################
# derivative
##################
merged <- read.csv("FSCV_deriv_150_straps_3_paths_convg10000_bagSize_2500")
merged.mean <- merged$mrg
bagSzs <- c(1,2,4,6,8,10,12,14,16,18,24,28,35,56,65,75,85,100,125,250, 1000)

avg <- colMeans( read.csv("fscv_deriv SS BagSize_Test__Avg") / merged.mean )
stack <- colMeans( read.csv("fscv_deriv SS BagSize_Test__Stack") / merged.mean )
cps <- colMeans( read.csv("fscv_deriv SS BagSize_Test__CPS") / merged.mean )

names(avg) <- names(stack) <- names(cps) <- bagSzs

df <- rbind(
    cbind(bagSzs, (avg), weights = "Average"),
    cbind(bagSzs, (stack), weights = "Stacking"),
    cbind(bagSzs, (cps), weights = "CPS")
)
df <- data.frame(df)
df <- cbind(df, "deriv")
colnames(df) <- c("bag", "rmse", "weights", "raw")


df <- data.frame( rbind(df2, df) )
df$bag <- factor(df$bag, levels = c(bagSzs))

df$rmse <- as.numeric( as.vector(df$rmse) )

##################
# plot
##################
for(pro in c("raw", "deriv")){
    p <- ggplot(data = df[df$raw == pro,], 
                        aes( x=factor(bag), y=log(rmse), color = factor(weights),
                               fill = factor(weights), group = factor(weights) ) 
                ) + 
        geom_line( size = 2) + 
        geom_point(size = 4, aes(shape = factor(weights ))) + 
        theme_classic(base_size = 12) + 
        xlab("Bag Size") +
        scale_x_discrete(breaks = df$bag[seq(1,21, 2)] )  +
        scale_y_continuous(labels=scaleFUN) +
        scale_color_manual(values = c("#ca0020", "#E69F00", "#0868ac")) + 
        scale_fill_manual(values = c("#ca0020", "#E69F00", "#0868ac")) +
        theme(legend.position = "none", 
              plot.title = element_text(hjust = 0.5, size=14, face="bold")) +
        ylab(TeX('$\\mathbf{log(RMSE/RMSE_{TOM})}$')) +
        theme(axis.text.x = element_text(face="bold", 
                                         angle=45, size=rel(2),
                                         vjust = 0.6 ),  
              axis.text.y = element_text(face="bold", size=rel(2) ), 
              axis.title.y = element_text(face="bold", size=rel(2),
                                          vjust = 3), 
              axis.title.x = element_text(face="bold", size=rel(2) )
              ) +
        geom_hline(yintercept=0, 
                   linetype="dashed", 
                   color = "black", 
                   size = rel(0.5),
                   alpha = 0.7)  +
        theme( plot.margin = margin(r = 20, b = 15) )  # add margins
    
        # add annotation if
    if(pro == "raw"){
        p <- p + annotate("text", x=16, y=1.6,  
                 label='bold("Stacking")',
                 size=10, color="#0868ac", parse=TRUE) + #  "#619CFF",
            annotate("text", x=16, y=1.35,  
                     label='bold("CPS")', 
                     size=10, color="#E69F00", parse=TRUE) +   # "#00BA38"  "#ca0020", , 
            annotate("text", x=16, y=1.85,  
                     label='bold("Average")', 
                     size=10, color="#ca0020", parse=TRUE)
    }else{
        p <- p + annotate("text", x=16, y=1.3,  
                          label='bold("Stacking")',
                          size=10, color="#0868ac", parse=TRUE) + #  "#619CFF",
            annotate("text", x=16, y=1.1,  
                     label='bold("CPS")', 
                     size=10, color="#E69F00", parse=TRUE) +   # "#00BA38"  "#ca0020", , 
            annotate("text", x=16, y=1.5,  
                     label='bold("Average")', 
                     size=10, color="#ca0020", parse=TRUE)
    }
   
        
    ##################
    # zoom
    ##################
    p2 <- p + coord_cartesian(ylim = c(-0.035, 0.05)) + 
                    theme( axis.title.y = element_blank() ) + 
        theme( plot.margin = margin(r = 30, b = 15, l = 15) )  # add margins to counterbalance other plot having space on left for y axis title
    
    p <- p + theme( plot.margin = margin(l = 15, r = 25, b = 15) )  # add margins for y-axis title space
    
    
    assign(paste0("bag.size.plot.", pro), p ) # unscaled
    assign(paste0("bag.size.plot.zoom.", pro), p2  ) # zoom
    
}

##################
# save
##################

setwd("~/Desktop/Research Final/Methods Paper/Draft/Figures") # setwd for saving file
# raw
ggsave("bagSize_SS_plotLogRaw.png", 
       plot = bag.size.plot.raw,
       width = 7, height = 4.5)
ggsave("bagSize_plot_zoomLogRaw.png",
       plot = bag.size.plot.zoom.raw,
       width = 7, height = 4.5)

# raw
ggsave("bagSize_SS_plotLogDeriv.png", 
       plot = bag.size.plot.deriv,
       width = 7, height = 4.5)
ggsave("bagSize_plot_zoomLogDeriv.png",
       plot = bag.size.plot.zoom.deriv,
       width = 7, height = 4.5)
