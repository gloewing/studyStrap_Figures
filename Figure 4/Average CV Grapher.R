#########################
# Average CVs -- color by height for cluster
library(viridis)
library(ggplot2)

title = "Average Electrodes (Cyclic Voltammograms)"
setwd("~/Desktop/Research Final/Methods Paper/Draft/Figures/Electrode Graphs")
cvs <- read.csv("Electrode Means")[,-1]
height <- rowSums(abs(cvs))
height <- height/max(height) - 0.25 # tune the color scale 
cv.means <- as.vector(t(cvs)) #transforms matrix into a vector
elec <- c(rep(1:15, each = 1000))
cv <- matrix(cbind(cv.means, factor(rep(1:1000, 15)), elec, factor(rep(1,15000)),
                   rep(height, each = 1000)), ncol = 5)
colnames(cv) <- c("current", "sample", "electrode", "line", "height")

row.names(cv) <- c()
volt.df <- cv
volt.df[,c(3,4)] <- as.factor(volt.df[,c(3,4)])
volt.df <- as.data.frame(volt.df)  

fscv.plot <- ggplot(volt.df, aes(x = sample, y = current, 
                                 color = factor(height))) +
    geom_line(aes(linetype = factor(line))) + 
    xlab("Covariate (Voltage Potential) Index") +
    ylab("Current (nA)") +
    scale_x_continuous(breaks = c(1, seq(200,1000,200))) +
    scale_y_continuous(breaks = seq(-1000,2000,500)) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), legend.position = "none",
          axis.title.y =  element_text(face="bold", 
                                       size=rel(2.5),
                                       vjust = 1.25),
          axis.title.x = element_text(face="bold", 
                                       size=rel(2.5),
                                      vjust = -1.75),
          axis.text.x = element_text(face="bold", 
                                     size=rel(2.5) ), 
          axis.text.y = element_text(face="bold", 
                                     size=rel(2.5) )) +
    labs(color = "Electrode") +
    scale_color_viridis(option = "D", discrete = TRUE) +
    theme( plot.margin = margin(l = 10, b = 20, r = 10) )  # add margins
    
fscv.plot
#######################
# 1000 600

ggsave("avg_cvs.pdf", 
       plot = fscv.plot,
       width = 9, height = 5)
